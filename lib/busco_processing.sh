#!/bin/bash

###########################################################################################
# buscoを実行
#  - ゲノムデータについて実行することが前提 -m genome 
#  - buscoはconda環境なので、関数内でbuscoのコマンド確認を行う
#  - run_busco <REF> <INPUT_FASTA> 
###########################################################################################
function run_busco () {
    # Usage: run_busco <reference_db> <genome.fasta> [<output_dir>] [<num_threads>]
    local REF=$1
    local IN_FA=$2
    local OUT_DIR=${3:-'out_busco'}
    local THRD=${4:-4}

    # コマンド
    command -v busco >/dev/null 2>&1 || { echo "[ERROR] Command busco not found!" ; return 1 ; }

    # サンプルプレフィクスを生成
    local ID
    ID=$(basename "$IN_FA" | sed -e "s/\..*//")

    # BUSCO出力ディレクトリの確認（既に存在する場合エラー）
    local RESULT_DIR="$OUT_DIR/$ID"
    if [[ -d "$RESULT_DIR" ]]; then
        echo "[ERROR] BUSCO result directory already exists: $RESULT_DIR" >&2
        return 1
    fi

    # 出力ディレクトリの作成
    if [[ ! -d "$OUT_DIR" ]]; then mkdir -p "$OUT_DIR"; fi

    # buscoコマンドを構築
    local CMD="busco -m genome -i \"$IN_FA\" --offline --out_path \"$OUT_DIR\" -o \"$ID\" -l \"$REF\" -c \"$THRD\""
    echo "[CMD] $CMD" >&2

    # コマンドの実行
    eval "$CMD" || { echo "[ERROR] Error in busco command: $CMD" >&2; return 1; }
}

###########################################################################################
# GNU parallelでrun_buscoを実行
#  - ゲノムデータについて実行することが前提 -m genome 
#  - run_busco_parallel <REF> <IN_DIR> <suffix_fasta> <OUT_DIR> <num_jobs>
#  - parallelのジョブ数は4で固定
###########################################################################################
function run_busco_parallel () {
    # Usage : run_busco_parallel <REF> <IN_DIR> <suffix_fasta> <OUT_DIR> <num_jobs>
    local REF=$1 
    local IN_DIR=$2
    local SFX=${3:-'fna'}
    local OUT_DIR=${4:-'out_busco'}
    local THREADS=${5:-4}
    local LOG_FILE
    local JOBS=4
    local PTHREADS

    # 入力ファイル取得
    local FAS
    mapfile -t FAS < <(find "$IN_DIR" -type f -name "*.$SFX")

    if [[ ${#FAS[@]} -eq 0 ]]; then 
        echo "[ERROR] No input files found in ${IN_DIR} with suffix ${SFX}" >&2
        return 1
    fi

    # 出力ディレクトリの確認と作成（既存の場合エラー）
    if [[ -d "$OUT_DIR" ]]; then
        echo "[ERROR] Output directory $OUT_DIR already exists. Please specify a new directory." >&2
        return 1
    else
        mkdir -p "$OUT_DIR"
    fi

    # プロセスのスレッド数
    PTHREADS=$(( THREADS / JOBS ))
    if [[ $PTHREADS -lt 1 ]]; then
        echo "[WARNING] Threads per job too low, setting to 1." >&2
        PTHREADS=1
    fi

    # ログファイルの準備
    LOG_FILE=$(date +"%Y%m%dT%H%M")_busco.log

    # run_busco 関数を export
    export -f run_busco
    export REF OUT_DIR JOBS PTHREADS LOG_FILE

    # 並列処理の実行
    parallel -j "$JOBS" --halt soon,fail=1 --no-run-if-empty run_busco "$REF" {} "$OUT_DIR" "$PTHREADS" 1>>"$LOG_FILE"  ::: "${FAS[@]}"

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Some BUSCO jobs failed. Check $LOG_FILE for details." >&2
        return 1
    fi

    echo "[INFO] BUSCO analysis completed. Results saved in $OUT_DIR." >&2
    echo "[INFO] Logs available in $LOG_FILE." >&2
}

###########################################################################################
# buscoのサマリーデータを整形
# - buscoの出力ディレクトリを入力として、buscoの結果を要約する
# - 完全性および汚染度の情報を追記する
###########################################################################################
function parse_busco_summary () {
    # Usage: parse_busco_summary <busco_dir>
    local IN_BSC=$1
    
    # IDの取得
    local ID
    ID=$(basename "$IN_BSC")
    
    # short_summary.txt のパスを確認
    local SSUM
    SSUM=$(find "$IN_BSC" -name "short_summary.txt" -print -quit)
    
    if [[ -z "$SSUM" ]]; then 
        echo "[WARNING] No short_summary.txt found in $IN_BSC" >&2
        return 1
    fi
    
    # short_summary.txt をパース
    awk -v id="$ID" '
        $0 ~ /Complete and single-copy BUSCOs/ {CBS=$1}
        $0 ~ /Complete and duplicated BUSCOs/  {CBD=$1}
        $0 ~ /Fragmented BUSCOs/              {FB=$1}
        $0 ~ /Missing BUSCOs/                 {MB=$1}
        END {
            SUM=CBS+CBD+FB+MB
            if (SUM > 0) {
                CMPL=(CBS/SUM)*100
                CNTM=(CBD/SUM)*100
            } else {
                CMPL=0
                CNTM=0
            }
            printf "%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n", id, CBS, CBD, FB, MB, CMPL, CNTM
        }
    ' "$SSUM"
    return 0
}

###########################################################################################
# buscoのサマリーデータを整形する
# - 複数のbusco結果をまとめて要約する
###########################################################################################
function summarize_busco () {
    # Usage: summarize_busco <busco_results_dir> [output_file]
    local BUSCO_DIR=$1
    local OUTPUT_FILE=${2:-}

    # 入力ディレクトリ確認
    if [[ ! -d "$BUSCO_DIR" ]]; then 
        echo "[ERROR] Directory $BUSCO_DIR does not exist" >&2
        return 1
    fi

    # ヘッダー作成
    local HEADER="ID\tComplete_Single\tComplete_Duplicated\tFragmented\tMissing\tCompleteness(%)\tContamination(%)"

    # 出力ファイル確認
    if [[ -n "$OUTPUT_FILE" ]]; then
        if [[ -e "$OUTPUT_FILE" ]]; then
            echo "[ERROR] Output file $OUTPUT_FILE already exists" >&2
            return 1
        fi
        echo -e "$HEADER" > "$OUTPUT_FILE"
    else
        echo -e "$HEADER"
    fi

    # サブディレクトリを並列で処理
    export -f parse_busco_summary
    export OUTPUT_FILE
    find "$BUSCO_DIR" -mindepth 1 -maxdepth 1 -type d | parallel --no-notice parse_busco_summary {} >> "${OUTPUT_FILE:-/dev/stdout}"

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Failed to process summaries of busco." >&2
        return 2
    fi

    return 0
}

###########################################################################################
# Core-geneの抽出
# - buscoの結果からComplete_Singleの配列を取得
# - summaryファイルを入力としゲノムのクォリティでフィルタ
# - ゲノムごとのbusco結果が格納されたディレクトリを入力として
# - 複数ゲノムに共通するComplete buscoを取得する
# - デフォルトで完全性50%以上、汚染度10%未満でゲノムをフィルタ
# - この時に採用されたゲノム数および除外されたゲノム数を標準エラー出力
###########################################################################################
function core_extraction () {
    # Usage: core_extraction <busco_summary.tsv> <busco_dir> [<completeness>] [<contamination>] [<output_dir>]
    local BUSCO_SUM=$1
    local BUSCO_DIR=$2
    local COMP=${3:-50}
    local CNTM=${4:-10}
    local OUT_CORE=${5:-out_core}

    # 入力チェック
    if [[ ! -f "$BUSCO_SUM" ]]; then
        echo "[ERROR] Summary file $BUSCO_SUM does not exist." >&2
        return 1
    fi
    if [[ ! -d "$BUSCO_DIR" ]]; then
        echo "[ERROR] Directory $BUSCO_DIR does not exist." >&2
        return 1
    fi
    if [[ -d "$OUT_CORE" ]]; then
        echo "[WARNING] Directory $OUT_CORE already exists. Remove it to rerun." >&2
        return 1
    fi

    mkdir -p "$OUT_CORE"

    # 条件を満たすゲノムを選択
    mapfile -t MAGS < <(awk -F"\t" -v comp="$COMP" -v cntm="$CNTM" '$6 > comp && $7 < cntm {print $1}' "$BUSCO_SUM")
    mapfile -t RMG < <(tail -n +2 "$BUSCO_SUM" | awk -F"\t" -v comp="$COMP" -v cntm="$CNTM" '$6 <= comp || $7 >= cntm {print $1}')
    
    local NMAGS=${#MAGS[@]}
    if [[ $NMAGS -eq 0 ]]; then
        echo "[ERROR] No genomes meet the criteria (Completeness >${COMP}%, Contamination <${CNTM}%)." >&2
        return 1
    else
        echo "[INFO] Accepted $NMAGS genomes and removed ${#RMG[@]} genomes." >&2
    fi

    # 全てのゲノムで共通するComplete BUSCOのIDを取得
    mapfile -t CBSC_ID < <(for MAG in "${MAGS[@]}"; do
        IN_DIR=$(find "${BUSCO_DIR}/${MAG}" -type d -name "run_*" -print -quit)
        if [[ -z "$IN_DIR" ]]; then
            echo "[WARNING] Missing BUSCO run directory for $MAG." >&2
            continue
        fi
        FTAB=$(find "$IN_DIR" -name "full_table.tsv" -print -quit)
        if [[ ! -f "$FTAB" ]]; then
            echo "[WARNING] Missing full_table.tsv in $IN_DIR." >&2
            continue
        fi
        grep -v "^#" "$FTAB" | awk -F"\t" '$2=="Complete" {print $1}'
    done | sort | uniq -c | awk -v nmags="$NMAGS" '$1 == nmags {print $2}')

    if [[ ${#CBSC_ID[@]} -eq 0 ]]; then
        echo "[ERROR] No common Complete BUSCOs found across all genomes." >&2
        return 1
    else 
        echo "[INFO] ${#CBSC_ID[@]} common Complete BUSCOs found across all genomes." >&2
    fi

    # 各ゲノムのComplete BUSCO配列を取得しマージ
    for CBSC in "${CBSC_ID[@]}"; do
        for MAG in "${MAGS[@]}"; do
            SEQ_DIR=$(find "${BUSCO_DIR}/${MAG}" -type d -path "*/busco_sequences/single_copy_busco_sequences" -print -quit)
            if [[ -z "$SEQ_DIR" ]]; then
                echo "[WARNING] Missing sequence directory for $MAG." >&2
                continue
            fi
            FNA_BUSCO="${SEQ_DIR}/${CBSC}.fna"
            FAA_BUSCO="${SEQ_DIR}/${CBSC}.faa"
            if [[ -f "$FNA_BUSCO" ]]; then
                sed -e "s/^>.*/>${MAG}/" "$FNA_BUSCO" >> "${OUT_CORE}/${CBSC}.fna"
            fi
            if [[ -f "$FAA_BUSCO" ]]; then
                sed -e "s/^>.*/>${MAG}/" "$FAA_BUSCO" >> "${OUT_CORE}/${CBSC}.faa"
            fi
        done
    done

    echo "[INFO] Complete BUSCO sequences have been extracted to $OUT_CORE."
    return 0
}
