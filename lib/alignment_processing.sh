#!/bin/bash

###########################################################################################
# Core-gene alignment and trimming
# - コアジーンセットをmafftでマルチプルアラインメント
# - アライメントfastaをtrimalでトリミング
# - core_align <input_fasta> <output_directory> [threads]
###########################################################################################
function core_align() {
    # Usage: core_align <input_fasta> <output_directory> [mafft_opts] [trimal_opts]
    # Example: core_align input.fna ./out --auto -automated1
    local IN_FA=$1
    local OUT_DIR=$2
    local MAFFT_OPTS=${3:-"--auto"}
    local TRIMAL_OPTS=${4:-"-automated1"}

    # 出力ディレクトリの作成
    if [[ ! -d "$OUT_DIR" ]]; then mkdir -p "$OUT_DIR"; fi


    # ファイル名と出力パス設定
    local ID ; ID=$(basename "$IN_FA" | sed -e "s/\..*//")
    local ALIGN_FILE="$OUT_DIR/${ID}.aln"
    local TRIM_FILE="$OUT_DIR/${ID}_trim.aln"    

    # MAFFT実行
    local CMD1="mafft $MAFFT_OPTS --quiet \"$IN_FA\" > \"$ALIGN_FILE\""
    echo "[CMD] $CMD1" >&2
    eval "$CMD1" || { echo "[ERROR] Error in mafft command: $CMD1" >&2; return 1; }

    # TRIMAL実行
    local CMD2="trimal -in \"$ALIGN_FILE\" $TRIMAL_OPTS > \"$TRIM_FILE\""
    echo "[CMD] $CMD2" >&2
    eval "$CMD2" || { echo "[ERROR] Error in trimal command: $CMD2" >&2; return 1; }


    echo "[INFO] Processed: $IN_FA -> $TRIM_FILE"
    return 0

}


###########################################################################################
# Core-gene alignment and trimming with GNU parallel
# - 塩基配列 (nuc) or アミノ酸配列 (aa) を選択し並列実行
# core_align_parallel --input input_dir --output output_dir --threads 8 --type aa --mafft-opts "--globalpair --maxiterate 1000" --trimal-opts "-gappyout"
###########################################################################################
function core_align_parallel() {
    # Usage: core_align_parallel <input_directory> <output_directory> [threads] [--type nuc|aa] \
    #                           [--mafft-opts "<options>"] [--trimal-opts "<options>"]

    local IN_DIR OUT_DIR="out_core_align" THRD=4 TYPE="nuc" MAFFT_OPTS="--auto" TRIMAL_OPTS="-automated1"

    # オプション解析
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                IN_DIR="$2"; shift 2 ;;
            -o|--output)
                OUT_DIR="$2"; shift 2 ;;
            -t|--threads)
                THRD="$2"; shift 2 ;;
            --type)
                TYPE="$2"; shift 2 ;;
            --mafft-opts)
                MAFFT_OPTS="$2"; shift 2 ;;
            --trimal-opts)
                TRIMAL_OPTS="$2"; shift 2 ;;
            *)
                echo "[ERROR] Unknown option: $1" >&2
                return 1 ;;
        esac
    done

    # 入力チェック
    if [[ -z "$IN_DIR" ]]; then echo "[ERROR] Input directory not specified. Use -i or --input" >&2 ; return 1 ; fi

    # 出力ディレクトリの作成
    [[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"

    # 配列ファイルのリスト取得
    local EXT="fna"  # デフォルトは塩基配列
    [[ "$TYPE" == "aa" ]] && EXT="faa"
    mapfile -t FASTA_FILES < <(find "$IN_DIR" -type f -name "*.$EXT")

    if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
        echo "[ERROR] No $TYPE FASTA files (.$EXT) found in $IN_DIR" >&2
        return 1
    fi

    # 並列処理
    export -f core_align
    export OUT_DIR
    export MAFFT_OPTS
    export TRIMAL_OPTS

    parallel -j "$THRD" --keep-order --halt soon,fail=1 \
        "core_align {} \"$OUT_DIR\" \"$MAFFT_OPTS\" \"$TRIMAL_OPTS\"" ::: "${FASTA_FILES[@]}"

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Parallel processing failed." >&2
        return 1
    fi

    echo "[INFO] All tasks completed successfully. Results are in $OUT_DIR"
    return 0
}



###########################################################################################
# Concatenate all core-gene alignments
# - core_concatenate <input_directory> <output_fasta> 
###########################################################################################
function core_concatenate () {
    # Usage: core_concatenate <input_directory> <output_fasta> 
    # Example: core_concatenate ./out_core_align concatenated.aln 

    local IN_DIR=$1
    local OUT_FASTA=${2:-'result_cores.aln'}

    # 入力ディレクトリの確認
    if [[ ! -d "$IN_DIR" ]]; then
        echo "[ERROR] Input directory $IN_DIR does not exist" >&2
        return 1
    fi

    # 各FASTAファイルを処理
    declare -A CONCATENATED
    for FASTA_FILE in "$IN_DIR"/*_trim.aln; do
        if [[ ! -f "$FASTA_FILE" ]]; then
            echo "[WARNING] No _trim.aln files found in $IN_DIR" >&2
            continue
        fi

        echo "[INFO] Processing $FASTA_FILE..." >&2
        local ID=""
        while IFS= read -r LINE; do
            if [[ "$LINE" == \>* ]]; then
                # 新しい配列のID行
                ID=${LINE#>}
                ID=${ID%% *} # IDの余計な空白を削除
            else
                # シーケンスを連結
                CONCATENATED["$ID"]+="$LINE"
            fi
        done < "$FASTA_FILE"
    done

    # 結果をFASTA形式で出力
    echo "[INFO] Writing concatenated sequences to $OUT_FASTA..." >&2
    : > "$OUT_FASTA" # ファイルを初期化
    for GENOME in "${!CONCATENATED[@]}"; do
        echo ">$GENOME" >> "$OUT_FASTA"
        echo "${CONCATENATED[$GENOME]}" >> "$OUT_FASTA"
    done

    echo "[INFO] Concatenation of all core-genes alignments completed: $OUT_FASTA" >&2

}
