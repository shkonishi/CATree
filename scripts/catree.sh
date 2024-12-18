#!/bin/bash
set -eC

# Constants
CMDNAME=$(basename "$0")
VERSION="0.1.0"
LOG_FILE=$(date +"%Y%m%dT%H%M")_${CMDNAME%%.*}.log

# Source common functions
source "$(dirname "$0")/../lib/functions_common.sh"
source "$(dirname "$0")/../lib/busco_processing.sh"
source "$(dirname "$0")/../lib/alignment_processing.sh"
source "$(dirname "$0")/../lib/tree_building.sh"

# Default parameters
CONDA_ENV_FILE="${HOME}/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_NAME="busco5"
BUSCO_REF="${HOME}/db/busco_downloads/lineages/bacteria_odb10"
SUFFIX_FASTA="fna"
OUTPUT_PREFIX="results"
COMPLETENESS=50
CONTAMINATION=10
THREADS=4
CONFIG_FILE="$(dirname "$0")/../config/id_lookup.tsv"

# Print help
function print_help() {
  cat << EOS
Usage: $CMDNAME [OPTIONS] <input_directory> <output_directory>
Options:
  -e, --conda_env <FILE>         Conda environment file (default: $CONDA_ENV_FILE)
  -b, --conda_name <NAME>        Conda environment name (default: $CONDA_ENV_NAME)
  -r, --reference <PATH>         BUSCO reference path (default: $BUSCO_REF)
  -s, --suffix <STRING>          Suffix of input files (default: $SUFFIX_FASTA)
  -m, --completeness <NUM>       Genome completeness threshold (default: $COMPLETENESS)
  -n, --contamination <NUM>      Genomic contamination threshold (default: $CONTAMINATION)
  -o, --output_prefix <STRING>   Output file prefix (default: $OUTPUT_PREFIX)
  -t, --threads <INT>            Number of threads (default: $THREADS)
  -c, --config <FILE>            Configuration file (default: $CONFIG_FILE)
  -h, --help                     Display this help and exit
  -v, --version                  Show script version
EOS
}

# Handle errors
function handle_error() {
  echo "[ERROR] $1" | tee -a "$LOG_FILE"
  exit 1
}

# Check dependencies
for cmd in conda mafft trimal FastTree; do
  command -v $cmd >/dev/null 2>&1 || handle_error "Command '$cmd' not found!"
done

# Parse options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -e|--conda_env) CONDA_ENV_FILE="$2"; shift 2 ;;
    -b|--conda_name) CONDA_ENV_NAME="$2"; shift 2 ;;
    -r|--reference) BUSCO_REF="$2"; shift 2 ;;
    -s|--suffix) SUFFIX_FASTA="$2"; shift 2 ;;
    -m|--completeness) COMPLETENESS="$2"; shift 2 ;;
    -n|--contamination) CONTAMINATION="$2"; shift 2 ;;
    -o|--output_prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -c|--config) CONFIG_FILE="$2"; shift 2 ;;
    -h|--help) print_help; exit 0 ;;
    -v|--version) echo "$VERSION"; exit 0 ;;
    *) # Positional arguments
      if [[ -z "$INPUT_DIR" ]]; then
        INPUT_DIR="$1"
      elif [[ -z "$OUTPUT_DIR" ]]; then
        OUTPUT_DIR="$1"
      else
        handle_error "Unexpected argument: $1"
      fi
      shift ;;
  esac
done

# Validate input
[[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]] && handle_error "Input and output directories are required."
[[ ! -d "$INPUT_DIR" ]] && handle_error "Input directory not found: $INPUT_DIR"
[[ -d "$OUTPUT_DIR" ]] && handle_error "Output directory already exists: $OUTPUT_DIR"

# Log setup
echo "$CMDNAME version $VERSION" | tee -a "$LOG_FILE"
echo "[CMD] $CMDNAME ${ARGS[*]}" | tee -a "$LOG_FILE"

# Arguments
cat << EOS | tee -a "$LOG_FILE" >&2
### Arguments ###
| Option                  | Value                      
|-------------------------|----------------------------
| Conda env name          | $CONDA_ENV_NAME           
| Conda env file          | $CONDA_ENV_FILE           
| Reference               | $BUSCO_REF                
| Suffix of input         | $SUFFIX_FASTA             
| Thresh of completeness  | $COMPLETENESS             
| Thresh of contamination | $CONTAMINATION            
| Input directory         | $INPUT_DIR                
| Output directory        | $OUTPUT_DIR               
| Prefix of output        | $OUTPUT_PREFIX            
| Num of threads          | $THREADS                  
| Config file             | $CONFIG_FILE              
EOS


# Main process
echo "[INFO] Activating conda environment $CONDA_ENV_NAME..." | tee -a "$LOG_FILE"
act_conda "$CONDA_ENV_NAME" "$CONDA_ENV_FILE" || handle_error "Failed to activate conda environment."

echo "[INFO] Running BUSCO..." | tee -a "$LOG_FILE"
BUSCO_DIR="$OUTPUT_DIR/busco"
run_busco_parallel "$BUSCO_REF" "$INPUT_DIR" "$SUFFIX_FASTA" "$BUSCO_DIR" "$THREADS" || handle_error "BUSCO failed."

echo "[INFO] Summarizing BUSCO results..." | tee -a "$LOG_FILE"
OUT_BUSCOSUM="$OUTPUT_DIR/busco_summary.tsv"
summarize_busco "$BUSCO_DIR" > "$OUT_BUSCOSUM" || handle_error "Failed to summarize BUSCO results."

echo "[INFO] Extracting core genes..." | tee -a "$LOG_FILE"
CORE_DIR="$OUTPUT_DIR/core_genes"
core_extraction "$OUT_BUSCOSUM" "$BUSCO_DIR" "$COMPLETENESS" "$CONTAMINATION" "$CORE_DIR" 2>&1 | tee -a "$LOG_FILE"
if [ "${PIPESTATUS[0]}" -ne 0 ]; then
    handle_error "Core gene extraction failed."
fi

echo "[INFO] Aligning core genes..." | tee -a "$LOG_FILE"
CORE_ALN_DIR="$OUTPUT_DIR/core_alignments"
core_align_parallel -i "$CORE_DIR" -o "$CORE_ALN_DIR" -t "$THREADS" || handle_error "Core gene alignment failed."

echo "[INFO] Concatenating alignments..." | tee -a "$LOG_FILE"
OUT_ALN="$OUTPUT_DIR/${OUTPUT_PREFIX}_core.aln"
core_concatenate "$CORE_ALN_DIR" "$OUT_ALN" || handle_error "Alignment concatenation failed."

echo "[INFO] Convert headers of alignment fasta..." | tee -a "$LOG_FILE"
convert_fasta_headers "$OUT_ALN" "$CONFIG_FILE" "$OUT_ALN" || echo "[ERROR] Conversion process failed."

echo "[INFO] Building phylogenetic tree..." | tee -a "$LOG_FILE"
run_fasttree "$OUT_ALN" || handle_error "Tree building failed."

echo "[INFO] Pipeline completed successfully!" | tee -a "$LOG_FILE"
