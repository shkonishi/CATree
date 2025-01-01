#!/bin/bash

###########################################################################################
# fastaのヘッダーを変換
#  - IDルックアップテーブル(タブ区切りテキスト)を予め作成しておく 
###########################################################################################
function convert_fasta_headers() {
    if [ $# -eq 0 ]; then echo "Usage: convert_fasta_headers <input.fasta> <lookup_table.tsv> <output.fasta>" ; return 1; fi
    local input_fasta=$1
    local lookup_table=$2
    local output_fasta=$3

    # 入力ファイルの存在確認
    if [[ ! -f "$input_fasta" ]]; then
        echo "[ERROR] Input fasta file not found: $input_fasta" >&2
        return 1
    fi
    if [[ ! -f "$lookup_table" ]]; then
        echo "[ERROR] Lookup table file not found: $lookup_table" >&2
        return 1
    fi

    # 一時ファイルの生成 (同名ファイルの上書きを防止)
    local temp_file="${output_fasta}.tmp"

    awk 'BEGIN {
        FS = "\t"; OFS = "\n";
    }
    NR == FNR {
        lookup[$1] = $2; 
        next;
    }
    /^>/ {
        header = substr($0, 2); 
        if (header in lookup) {
            print ">" lookup[header];
        } else {
            print ">" header;
        }
        next;
    }
    {
        print;
    }' "$lookup_table" "$input_fasta" > "$temp_file"

    # awk 実行結果の確認
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Failed to process fasta headers" >&2
        rm -f "$temp_file"
        return 1
    fi

    # 出力ファイルを移動
    mv "$temp_file" "$output_fasta" || {
        echo "[ERROR] Failed to move temporary file to $output_fasta" >&2
        rm -f "$temp_file"
        return 1
    }

    echo "[INFO] Header conversion completed successfully: $output_fasta"
    return 0

}

###########################################################################################
# FastTreeを実行する
# # https://microbesonline.org/fasttree/#BranchLen
# # https://purduercac-applications.readthedocs.io/en/latest/Biocontainers/fasttree/fasttree.html
# run_fasttree input_nuc.aln "-nt -gtr"
# run_fasttree input_aa.aln "-lg" # -lg: LG+CAT（amino acid）
###########################################################################################
function run_fasttree () {
    if [ $# -eq 0 ]; then echo "Usage: run_fasttree <ALIGNMENT_FASTA> [<FASTTREE_OPTS>]" ; return 1; fi
    local IN_ALN=$1
    local FASTTREE_OPTS=${2:-'-nt'}
    local PFX OUT_DIR OUT_TRE LOG_FT 

    # 出力ファイル名の定義
    PFX=$(basename "${IN_ALN%.*}")
    OUT_DIR=$(dirname "$IN_ALN")
    OUT_TRE="${OUT_DIR}/${PFX}.nwk"


    # 入力ファイルの存在チェック
    if [[ ! -f "$IN_ALN" ]]; then
        echo "[ERROR] Input alignment $IN_ALN does not exist" >&2
        return 1
    fi

    # 出力ファイルの上書き防止チェック
    if [[ -f "$OUT_TRE" ]]; then
        echo "[ERROR] Output tree $OUT_TRE already exists" >&2
        return 1
    fi

    # ログファイルの準備
    LOG_FT="${OUT_DIR}/$(date +"%Y%m%dT%H%M")_fasttree.log"

    echo "$LOG_FT" "$OUT_TRE" 

    # # FastTreeコマンドの実行
    cmd="FastTree ${FASTTREE_OPTS} -quiet -log ${LOG_FT} ${IN_ALN} > ${OUT_TRE}"
    echo "[CMD] $cmd" >&2
    eval "$cmd" || {
        echo "[ERROR] FastTree execution failed for $IN_ALN" >&2
        return 1
    }

    return 0
}
