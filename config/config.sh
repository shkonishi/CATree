#!/bin/bash

# Default parameters
# shellcheck disable=SC2034
CONDA_ENV_FILE="${HOME}/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_NAME="busco5"
BUSCO_REF="${HOME}/db/busco_downloads/lineages/bacteria_odb10"
SUFFIX_FASTA="fna"
OUTPUT_PREFIX="results"
COMPLETENESS=50
CONTAMINATION=10
THREADS=4
TYPE='nuc'
MAFFT_OPTS="--auto" # "--globalpair --maxiterate 1000"
TRIMAL_OPTS="-automated1" # "-gappyout"
CONFIG_FILE="$(dirname "$0")/../config/id_lookup.tsv"