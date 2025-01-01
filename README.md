# CATree: Core-gene Alignment and Tree Construction
CATree is a bash-based pipeline to construct phylogenetic trees from core genes of bacterial genomes. It automates core gene detection using BUSCO, followed by alignment, trimming, and tree building.

## Requirement
Ensure the following dependencies are installed and accessible in your PATH:

- BUSCO >= 5.0 (conda installation required)
- MAFFT >= 7.490
- TrimAl >= 1.4
- FastTree >= 2.1
- (Optional, future support planned) CheckM

## Installation
Clone this repository into your local machine
```bash
git clone https://github.com/shkonishi/CATree
cd CATree
```

## Running with Example Data
### Download the example data 
Download the example genome data using the provided script:
```bash
cd examples
bash ./download_examples.sh
```

### Basic Usage : 
If you want to use the default values for all options, simply run:
```bash
bash ../scripts/catree.sh ./genomes ./output
```

### Example with All Options:
```bash
bash ./scripts/catree.sh \
  -e ${HOME}/miniconda3/etc/profile.d/conda.sh \  # Path to conda environment activation script
  -b busco5 \                                    # BUSCO environment name (default: busco5)
  -r ${HOME}/db/busco_downloads/lineages/bacteria_odb10 \  # BUSCO lineage dataset (default: bacteria_odb10)
  -m 50 \                                        # Completeness threshold (default: 50)
  -n 10 \                                        # Contamination threshold (default: 10)
  -y nuc \                                       # Type of coregene sequence nuc or aa (default: nuc)
  -o result \                                    # Prefix for output files (default: "result")
  -t 4 \                                         # Number of threads (default: 4)
  -c ./config/id_lookup.tsv \                    # Path to ID lookup table (default: ./config/id_lookup.tsv)
  ./examples/genomes ./output_dir                # Input genomes directory and output directory

```
### Notes:
- **-e and -r option** : 
If the default paths for -e and -r do not exist, create a symbolic link to the appropriate location.

- **id_lookup.tsv** : 
The file `config/id_lookup.tsv` is a conversion table for rewriting header lines in FASTA files after multiple sequence alignment.  
Important: Periods (.) are not allowed in either the genome ID or the leaf label.
```txt
<Genome_ID>    <Leaf_Label>
GCA_000205025    Parasutterella_excrementihominis_GCA_000205025
```
- **Applicability to Eukaryotes** : 
While this pipeline is designed for prokaryotic genomes, it is also likely applicable to eukaryotic genomes.


## Future Features
I plan to implement the following features in future updates:
- Support for transcriptome data to core-gene extraction with BUSCO.
- Option to use CheckM for core gene extraction instead of BUSCO.

## Citation
If you use CATree in your research, please cite the following paper: (under preparation)
- Author Name et al. (Year) Title of the paper.

## Contact
If you have any questions, suggestions, or issues related to this project, please feel free to open a GitHub Issue. We will do our best to respond promptly.