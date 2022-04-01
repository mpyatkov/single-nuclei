#!/bin/bash

#module load nextflow
module load nextflow/21.04.2

NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf \
	-c configuration/nextflow.config \
	-entry module_1 -resume 
