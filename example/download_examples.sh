#!/bin/bash

# Exit script on error
set -e
set -o pipefail

# Define input directory
SCRIPT_DIR="$(realpath "$(dirname "$0")")"
IN_DIR="${SCRIPT_DIR}/genomes"

# Create directory for demo-data
mkdir -p "$IN_DIR"
echo "Downloading example genome data into ${IN_DIR}..."

# Define URLs and filenames
declare -A GENOMES=(
    ["GCA_000205025.1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/205/025/GCA_000205025.1_ASM20502v1/GCA_000205025.1_ASM20502v1_genomic.fna.gz"
    ["GCA_000250875.1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/250/875/GCA_000250875.1_ASM25087v1/GCA_000250875.1_ASM25087v1_genomic.fna.gz"
    ["GCA_001411495.1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/411/495/GCA_001411495.1_ASM141149v1/GCA_001411495.1_ASM141149v1_genomic.fna.gz"
    ["GCA_003315195.1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/315/195/GCA_003315195.1_ASM331519v1/GCA_003315195.1_ASM331519v1_genomic.fna.gz"
    ["GCA_003402575.1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/402/575/GCA_003402575.1_ASM340257v1/GCA_003402575.1_ASM340257v1_genomic.fna.gz"
    ["GCA_003609995.1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/609/995/GCA_003609995.1_ASM360999v1/GCA_003609995.1_ASM360999v1_genomic.fna.gz"
    ["GCA_006337085.1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/337/085/GCA_006337085.1_ASM633708v1/GCA_006337085.1_ASM633708v1_genomic.fna.gz"
    ["GCA_016696765.1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/696/765/GCA_016696765.1_ASM1669676v1/GCA_016696765.1_ASM1669676v1_genomic.fna.gz"
)

# Download genomes
for ID in "${!GENOMES[@]}"; do
    FILE="${IN_DIR}/${ID}.fna.gz"
    if [[ -f "$FILE" ]]; then
        echo "File ${FILE} already exists. Skipping download."
    else
        echo "Downloading ${ID}..."
        wget -q -O "$FILE" "${GENOMES[$ID]}"
    fi
done

# Decompress downloaded files
echo "Decompressing downloaded files..."
gunzip -f "${IN_DIR}"/*.fna.gz
echo "All genomes downloaded and decompressed successfully!"