#!/bin/bash
set -eC

# Constants
CMDNAME=$(basename "$0")
VERSION="0.1.0"
LOG_FILE=$(date +"%Y%m%dT%H%M")_${CMDNAME%%.*}.log

# Source common functions
source "$(dirname "$0")/../lib/functions_common.sh"
source "$(dirname "$0")/../lib/rrna_extraction.sh"
source "$(dirname "$0")/../lib/alignment_processing.sh"
source "$(dirname "$0")/../lib/tree_building.sh"

# Default parameters
SUFFIX_FASTA="fna"
OUTPUT_PREFIX="results_16s"
MODE='bac'
THREADS=4
IDENTITY='0.97'
MAFFT_OPTS="--auto" # "--globalpair --maxiterate 1000"
TRIMAL_OPTS="-automated1" # "-gappyout"
CONFIG_FILE="$(dirname "$0")/../config/id_lookup.tsv"
FTREE_OPTS="-nt -gtr"

# Print help
function print_help() {
  cat << EOS
Usage: $CMDNAME [OPTIONS] <input_directory> <output_directory>
Options:
  -s, --suffix <STRING>          Suffix of input files (default: $SUFFIX_FASTA)
  -o, --output_prefix <STRING>   Output file prefix (default: $OUTPUT_PREFIX) 
  -i, --identity <NUMERIC>       reject if identity lower (default: 0.97)       
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
for cmd in barrnap vsearch mafft trimal FastTree; do
  command -v $cmd >/dev/null 2>&1 || handle_error "Command '$cmd' not found!"
done

# Parse options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--suffix) SUFFIX_FASTA="$2"; shift 2 ;;  
    -o|--output_prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
    -i|--identity) IDENTITY="$2"; shift 2 ;;
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
if [[ -z "$INPUT_DIR" && -z "$OUTPUT_DIR" ]] ; then print_help ; exit 1 ; fi
[[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]] && handle_error "Input and output directories are required."
[[ ! -d "$INPUT_DIR" ]] && handle_error "Input directory not found: $INPUT_DIR"
[[ -d "$OUTPUT_DIR" ]] && handle_error "Output directory already exists: $OUTPUT_DIR"

# Log setup
echo "$CMDNAME version $VERSION" | tee -a "$LOG_FILE"
echo "[CMD] $CMDNAME $*" | tee -a "$LOG_FILE"

# Arguments
cat << EOS | tee -a "$LOG_FILE" >&2
### Arguments ###
| Option                  | Value                      
|-------------------------|----------------------------     
| Suffix of input         | $SUFFIX_FASTA             
| Input directory         | $INPUT_DIR                
| Output directory        | $OUTPUT_DIR               
| Prefix of output        | $OUTPUT_PREFIX            
| Identity of clustering  | $IDENTITY                  
| Options for mafft       | $MAFFT_OPTS		      
| Options for TrimAL      | $TRIMAL_OPTS	      
| Options for FastTree    | $FTREE_OPTS    
| Num of threads          | $THREADS                  
| Config file             | $CONFIG_FILE              
EOS

# Main process
# 16S rRNA specific workflow
echo "[INFO] Starting 16S rRNA extraction..." | tee -a "$LOG_FILE"
OUTRRNA_DIR="${OUTPUT_DIR}/rrna"
batch_extract_unique16s "$INPUT_DIR" "$OUTRRNA_DIR" "$SUFFIX_FASTA" "$MODE" "$IDENTITY" "$THREADS" || handle_error "Failed to extract 16S rRNA sequences."

# Merge all 16S rRNA copies
echo "[INFO] Performing merged fasta ..." | tee -a "$LOG_FILE"
U16S_FA="${OUTPUT_DIR}/merged_16s.fa"
merge_fasta "$OUTPUT_DIR" "$U16S_FA" || handle_error "Failed to merge sequences."

# Alignment all 16S rRNA copies
echo "[INFO] Alignment & trimming ..." | tee -a "$LOG_FILE"
core_align "$U16S_FA" "$OUTPUT_DIR" || handle_error "Failed to align sequences."

# Header exchange
U16S_ALN=$(find "$OUTPUT_DIR" -type f -name "*_trim.aln")
if [[ ! -f "$U16S_ALN" ]] ; then handle_error "Cannot find fasta after alignment and trimming" ; fi
OUT_ALN="$OUTPUT_DIR/${OUTPUT_PREFIX}.aln"

echo "[INFO] Convert headers of alignment fasta..." | tee -a "$LOG_FILE"
convert_fasta_headers "$U16S_ALN" "$CONFIG_FILE" "$OUT_ALN" || echo "[ERROR] Conversion process failed."

# Tree building
echo "[INFO] Building phylogenetic tree of 16s rRNA ..." | tee -a "$LOG_FILE"
run_fasttree "$OUT_ALN" "$FTREE_OPTS" || handle_error "Tree building failed."

echo "[INFO] 16S rRNA analysis completed successfully." | tee -a "$LOG_FILE"

