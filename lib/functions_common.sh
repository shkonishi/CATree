#!/bin/bash
###############################################################################
# ここには共通関数が記述されています。
# スクリプト内で使用する場合は以下のようにしてインクルードする必要があります。
# source "$(dirname "$0")/../lib/common.sh"
###############################################################################

# conda 環境を実行する
function act_conda() {
    # Usage: act_conda <conda_env_name> [<conda_env_file>]
    # Example: 
    #   act_conda qiime2-2022.2 ${HOME}/miniconda3/etc/profile.d/conda.sh
    #   if [[ $? -ne 0 ]]; then echo "アクティベートに失敗しました" ; exit 1; fi

    # 引数の取得
    local env_name="$1"
    local conda_env_file="${2:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

    # conda環境名が指定されているか確認
    if [[ -z "$env_name" ]]; then
        echo "[ERROR] Please specify the conda environment name." >&2
        return 1
    fi

    # conda環境ファイルが存在するか確認
    if [[ ! -f "$conda_env_file" ]]; then
        echo "[ERROR] Conda environment file not found: $conda_env_file" >&2
        return 1
    fi

    # condaを初期化
    # shellcheck disable=SC1090
    source "$conda_env_file"
    if ! command -v conda &> /dev/null; then
        echo "[ERROR] Failed to initialize conda. Please check the environment file: $conda_env_file" >&2
        return 1
    fi

    # すでに指定された環境がアクティベートされている場合
    if [[ "$CONDA_DEFAULT_ENV" == "$env_name" ]]; then
        echo "[INFO] Conda environment is already activated: $env_name"
        return 0
    fi

    # conda環境をアクティベート
    conda activate "$env_name"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Failed to activate conda environment: $env_name" >&2
        return 1
    fi

    echo "[INFO] Successfully activated conda environment: $env_name"
    return 0
}
