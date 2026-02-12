#!/bin/bash

# Run the single-nuclei RNA-seq pipeline
# Usage: ./pipeline.sh [options]

set -eu

module load nextflow/21.10.6

NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf -c nextflow.config -resume

module unload nextflow
