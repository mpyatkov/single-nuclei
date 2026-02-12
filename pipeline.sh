#!/bin/bash

# Run the single-nuclei RNA-seq pipeline
# Usage: ./pipeline.sh

set -eu

# Convert samples.xlsx to samples.csv
echo "Converting samples.xlsx to samples.csv..."
module load miniconda
conda activate /projectnb2/wax-es/routines/condaenv/rlang4
Rscript bin/convert_samples_xlsx.R samples.xlsx --output samples.csv
conda deactivate
module unload miniconda

# Run Nextflow pipeline
module load nextflow/21.10.6
NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf -c nextflow.config -resume
module unload nextflow
