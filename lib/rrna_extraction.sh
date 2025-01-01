#!/bin/bash

#######################################################################
# predict_rrna
#  barrnapを使用して単一ゲノムからrRNA配列を抽出する
#  barrnapの出力fastaにはpartialの情報が含まれないので、自前で抽出D
#######################################################################
function predict_rrna () {
    if [ $# -eq 0 ]; then echo "Usage: predict_rrna <input.fasta> <output.gff> [arc|bac]" ; return 1; fi
    local input=$1
    local output=$2
    local mode=${3:-'bac'}
    local prefix

    # 入出力チェック
    if [[ ! -f "$input" ]]; then echo "ERROR: Input file $input not found." >&2 ; return 1 ; fi
    if [[ -z "$output" ]]; then prefix=$(basename "${input%.*}"); output="${prefix}_barrnap.gff" ; fi
    if [[ -f "$output" ]]; then echo "ERROR: Output file $output aleady exist." >&2 ; return 1 ; fi
    if [[ $mode != "arc" && $mode != "bac" ]]; then echo "ERROR: Invalid mode. Use 'arc' or 'bac'." >&2 ; return 1 ; fi

    # barrnapの実行
    cmd="barrnap --kingdom $mode $input > $output"
    echo "[CMD] $cmd" >&2
    eval "$cmd" || return 1
    echo -e "[INFO] Barrnap completed: $output" >&2
    return 0

}

#######################################################################
# predict_rrna_parallel
#  barrnapを使用して複数ゲノムを並列で処理
#######################################################################
function predict_rrna_parallel () {
    if [ $# -eq 0 ]; then echo "Usage: predict_rrna_parallel <input_dir> [suffix] [arc|bac] [output_dir] [jobs] " ; return 1; fi
    local input_dir=$1
    local suffix=${2:-'fna'}
    local mode=${3:-'bac'}
    local output_dir=${4:-./out_rrna}
    local njobs=${5:-4}
    local genomes

    # 入力ディレクトリと出力ディレクトリの確認
    if [[ ! -d "$input_dir" ]]; then echo "Error: Input directory $input_dir not found." >&2 ; return 1 ; fi
    if [[ ! -d "$output_dir" ]]; then mkdir -p "$output_dir" ; fi

    # ファイルリストの取得
    mapfile -t genomes < <(find "$input_dir" -type f -name "*.${suffix}")
    if [[ ${#genomes[@]} -eq 0 ]]; then echo "Error: No files with suffix .$suffix found in $input_dir." >&2 ; return 1 ; fi
    if [[ $mode != "arc" && $mode != "bac" ]]; then echo "Error: Invalid mode. Use 'arc' or 'bac'." >&2 ; return 1 ; fi

    # 並列処理の実行
    export suffix mode output_dir
    export -f predict_rrna
    parallel --no-notice --jobs "$njobs" predict_rrna {} "$suffix" "$mode" "$output_dir" ::: "${genomes[@]}"
    
}

#######################################################################
# extract_gff: gffから指定文字列を検索してその領域の配列を抽出
#  デフォルトでは文字列"product=16S ribosomal RNA"を検索する
#  部分配列の場合はヘッダに追記される >Contig_98_38_1144_-_partial
#######################################################################
function extract_gff () {
  if [ $# -eq 0 ]; then echo "Usage: extract_gff <input.gff> <input.fna> [search_string]" ; return 1; fi
  local IN_GFF=$1 
  local IN_FAS=$2
  local STR=${3:-"product=16S ribosomal RNA"}

  # 入力ファイルの確認
  if [[ ! -f "$IN_GFF" || ! -f "$IN_FAS" ]]; then
    echo "Error: Missing input files." >&2
    return 1
  fi

  # reverse complement
  function revcomp() {
    echo "$1" | tr ATGCYRKMWBVDHNatgcyrkmwbvdhn TACGRYMKWVBHDNtacgrymkwvbhdn | rev
  }

  # GFFファイルの処理
  grep -v "^#" "$IN_GFF" | awk -F"\t" -v str="$STR" '
    $3 ~ /rRNA/ && $9 ~ str { 
      partial = ($9 ~ /partial/) ? "_partial" : "";
      print $1, $4, $5, $7, partial
    }
  ' | while read -r CNTG ST ED STRND PARTIAL; do
    if [[ -z "$CNTG" || -z "$ST" || -z "$ED" || -z "$STRND" ]]; then
      echo "Error: Invalid GFF entry for CNTG=$CNTG, ST=$ST, ED=$ED, STRND=$STRND" >&2
      continue
    fi

    # 対応するFASTA配列の抽出
    SEQ=$(awk -v c="$CNTG" -v st="$ST" -v ed="$ED" '
      BEGIN { seq_found=0; seq="" }
      /^>/ { 
        if (seq_found) exit; 
        if ($0 ~ ">"c) seq_found=1; next 
      }
      seq_found { seq=seq $0 }
      END { print substr(seq, st, ed-st+1) }
    ' "$IN_FAS")

    if [[ -z "$SEQ" ]]; then
      echo "Error: Sequence not found for $CNTG:$ST-$ED" >&2
      continue
    fi

    # ストランド方向に応じて配列を修正
    if [[ "$STRND" == "-" ]]; then
      SEQ=$(revcomp "$SEQ")
    fi

    # 出力（部分配列の場合は"_partial"を付与）
    echo -e ">${CNTG}_${ST}_${ED}_${STRND}${PARTIAL}\n${SEQ}"
  done
  return 0
}

#######################################################################
# filter_16s : マルチコピーとして予測された遺伝子配列をフィルタする
#  1.全て全長配列の場合: そのまま出力
#  2.全長配列と部分配列が混在の場合: 全長配列のみ出力 
#    除外されたfastaのヘッダー行を標準エラー出力に
#  3.全て部分配列の場合: 最も長い配列を出力
#  filter_16s <STDIN>
#######################################################################
function filter_16s() {
    local in_fa
    local temp_fa

    # 標準入力を処理する場合
    temp_fa=$(mktemp)
    cat > "$temp_fa"
    in_fa="$temp_fa"

    # 内部関数: 最も長い配列のみを返す
    function longest_fa() {
        awk '
        BEGIN { max_length = 0 }
        /^>/ {
            if (seq) {
                if (length(seq) > max_length) {
                    if (longest_header) {
                        print "[INFO] Reject: " longest_header > "/dev/stderr"
                    }
                    max_length = length(seq)
                    longest_header = current_header
                    longest_seq = seq
                } else {
                    print "[INFO] Reject: " current_header > "/dev/stderr"
                }
            }
            current_header = $0
            seq = ""
            next
        }
        { seq = seq $0 }
        END {
            if (length(seq) > max_length) {
                if (longest_header) {
                    print "[INFO] Reject: " longest_header > "/dev/stderr"
                }
                longest_header = current_header
                longest_seq = seq
            } else {
                print "[INFO] Reject: " current_header > "/dev/stderr"
            }
            print longest_header
            print longest_seq
        }'
    }

    # 処理: barrnapの結果から抽出された16SrRNAのfastaのヘッダ行を解析
    mapfile -t header_all < <(grep -e ">" "$in_fa")
    mapfile -t header_part < <(grep -e ">" "$in_fa" | grep -e "_partial")
    mapfile -t header_full < <(grep -e ">" "$in_fa" | grep -v "_partial")

    if [[ ${#header_full[@]} == "${#header_all[@]}" ]]; then
        # 全て全長配列の場合: そのまま出力
        cat "$in_fa"
    elif [[ ${#header_full[@]} -gt 0 && ${#header_part[@]} -gt 0 && ${#header_full[@]} != "${#header_all[@]}" ]]; then
        # 全長配列と部分配列が混在の場合: 全長配列のみ出力
        for i in "${header_full[@]}"; do 
            awk -v search="$i" '
            BEGIN { found = 0; }
            /^>/ { if (index($0, search) > 0) { found = 1; print; } else { found = 0; } }
            !/^>/ { if (found) print; }
            ' "$in_fa"
        done
        # 除外されたfastaのヘッダー行を標準エラー出力に
        for rj in "${header_part[@]}"; do 
            echo "Rejects: $rj" >&2
        done
    elif [[ ${#header_full[@]} == 0 && ${#header_part[@]} -gt 0 ]]; then
        # 全て部分配列の場合: 最も長い配列を出力
        longest_fa < "$in_fa"
    else
        echo "Error: Unexpected condition encountered." >&2
        return 2
    fi

    # 一時ファイルを削除
    [[ -n "$temp_fa" ]] && rm -f "$temp_fa"
}

#######################################################################
# unique_fa
# - genome データから抽出したrRNA配列を入力とする
# - 16S rRNAのみ抽出し、vsearchでクラスタリング
# - fastaのヘッダ行を変更する(例. >GCA_XXXXX_cp1)
# - IDの変換履歴を標準エラー出力
#######################################################################
function unique_fa () {
    if [ $# -eq 0 ]; then echo "Usage: unique_fa <input.fasta> [<output_dir>] [<thresh_of_identity>]" ; return 1; fi
    local FA=$1
    local OUT_DIR=${2:-./out_16s}
    local THRESH=${3:-'0.97'}
    local PFX OUTPUT OUTUC VSLOGS
    

    # 入力ファイルの存在確認
    if [[ -z "$FA" || ! -f "$FA" ]]; then
        echo "Error: Input file '$FA' not found or not specified." >&2
        return 1
    fi

    # 出力ディレクトリの作成
    if [[ ! -d "$OUT_DIR" ]]; then mkdir -p "$OUT_DIR"; fi

    # 出力ファイル名の設定
    PFX=$(basename "${FA%.*}")
    OUTPUT="${OUT_DIR}/${PFX}_ucopy.fa"
    OUTUC="${OUT_DIR}/${PFX}_clusters.uc"
    VSLOGS="${OUT_DIR}/combined_vsearch.log"

    # クラスタリング処理
    vsearch --cluster_fast "$FA" --id "$THRESH" --centroids - --uc "$OUTUC" 2>>"$VSLOGS" \
    | awk -v pfx="$PFX" '
        BEGIN { n = 0; }
        /^>/ {
            n++;
            new_id = ">" pfx "_cp" n;
            print pfx "\t" new_id "\t" $0 > "/dev/stderr";
            print new_id;
        }
        !/^>/ { print; }
    ' > "$OUTPUT"
    
    return 0
}

#######################################################################
# extract_unique16s <input_dir> 
#######################################################################
function extract_unique16s() {
    if [ $# -eq 0 ]; then 
        echo "Usage: extract_unique16s <input_dir> [<output_dir>] [<suffix>] [<bac|arc>] [<identity>] [<threads>]" >&2
        return 1
    fi
    local input_dir=$1
    local output_dir=${2:-"./out_rrna"}
    local suffix=${3:-'fna'}
    local mode=${4:-'bac'}
    local identity=${5:-'0.97'}
    local njobs=${6:-4}

    # 入力ディレクトリの存在確認
    if [[ -z "$input_dir" || ! -d "$input_dir" ]]; then
        echo "[ERROR] Either $input_dir is missing, not a directory, or not defined!" >&2
        return 1
    fi

    # 出力ディレクトリの作成
    if [[ ! -d "$output_dir" ]]; then mkdir -p "$output_dir"; fi

    # ファイルリストの取得
    mapfile -t genomes < <(find "$input_dir" -type f -name "*.${suffix}")
    if [[ ${#genomes[@]} -eq 0 ]]; then
        echo "[ERROR] No files with suffix .$suffix found in $input_dir." >&2
        return 1
    fi

    # モードのバリデーション
    if [[ $mode != "arc" && $mode != "bac" ]]; then
        echo "[ERROR] Invalid mode. Use 'arc' or 'bac'." >&2
        return 1
    fi

    # 関数内のログメッセージを処理
    # shellcheck disable=SC2317    
    function process_genome() {
        local genome="$1"
        local output_dir="$2"
        local suffix="$3"
        local identity="$4"

        local pfx ; pfx=$(basename "$genome" ."$suffix")
        local out_gff="${output_dir}/${pfx}_barrnap.gff"
        local out_16s="${output_dir}/${pfx}_filt_16s.fa"
        local log_file="${output_dir}/${pfx}_processing.log"

        {
            echo "[INFO] Processing genome: $genome"
            predict_rrna "$genome" "$out_gff" &&
            echo "[INFO] Barrnap completed: $out_gff" &&
            extract_gff "$out_gff" "$genome" | filter_16s > "$out_16s" &&
            echo "[INFO] Filtering completed: $out_16s" &&
            unique_fa "$out_16s" "$output_dir" "$identity" &&
            echo "[INFO] Unique 16S sequences saved for $genome"
        } >> "$log_file" 2>&1 || echo "[ERROR] Processing failed for $genome" >> "$log_file"
    }
    export -f process_genome predict_rrna extract_gff filter_16s unique_fa

    # 並列処理
    parallel --jobs "$njobs" --halt soon,fail=1 --line-buffer process_genome {} "$output_dir" "$suffix" "$identity" ::: "${genomes[@]}"

    echo "[INFO] All processing completed. Check logs in $output_dir"
}

#######################################################################
# merge_ufasta <input_dir> <output_fasta>
# - 抽出された16s-rrnaのfastaをマージしてヘッダ行を変更する
# - 複数コピーが存在する場合は警告を出す 
#######################################################################
function merge_fasta() {
    # 引数のチェック
    if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
        echo "Usage: merge_fasta <input_dir> <output_file> [file_extension]" >&2
        return 1
    fi

    local input_dir="$1"
    local output_file="$2"
    local file_extension="${3:-_filt_16s_ucopy.fa}"  # デフォルト値を設定

    # 出力ファイルを初期化
    : > "$output_file"

    # 入力ディレクトリ内の指定された拡張子のFASTAファイルを処理
    while IFS= read -r -d '' file; do
        # シーケンス数を確認
        local num_sequences
        num_sequences=$(grep -c "^>" "$file")

        if [ "$num_sequences" -eq 1 ]; then
            # ヘッダーを変更（ファイル名から識別子を抽出）
            awk '/^>/{sub(".[0-9]_filt_16s_.*","",$1)}{print}' "$file" >> "$output_file"
        else
            # 警告を出力し、ヘッダーは変更せずそのまま出力
            echo "Warning: File '$file' contains $num_sequences sequences. Including original headers." >&2
            cat "$file" >> "$output_file"
        fi
    done < <(find "$input_dir" -type f -name "*${file_extension}" -print0)
}

