#!/bin/bash

#module load nextflow
module load nextflow/21.10.6

NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf -c configuration/nextflow.config -entry module2_preconfiguration -resume 
