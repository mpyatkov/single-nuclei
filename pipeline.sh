#!/bin/bash

# example1: ./pipeline.sh -s module2 # show module 2 parameters
# example2: ./pipeline.sh -r module2 # run module 2

set -eu

module load nextflow/21.10.6

VALID_ARGS=$(getopt -o s:r:h --long show:,run:,help -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

COMMAND=""

eval set -- "$VALID_ARGS"
while [ : ]; do
    case "$1" in
	-s | --show)
            # echo "Processing 'show' option. Input argument is '$2'"
	    COMMAND="NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf -c configuration/nextflow.config --current_module $2 -entry show_parameters -resume"
            shift 2
            ;;
	-r | --run)
            # echo "Processing 'run' option. Input argument is '$2'"
	    COMMAND="NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf -c configuration/nextflow.config -entry $2 -resume "
            shift 2
            ;;
	-h | --help)
	    echo "Usage: $(basename $0) [-r|run] [-s|show] <module#>"
	    exit 0
	    ;;

	:)
	    echo -e "option requires an argument.\nUsage: $(basename $0) [-r|run] [-s|show] <module#>"
	    exit 1
	    ;;
	?)
   	    echo -e "option requires an argument.\nUsage: $(basename $0) [-r|run] [-s|show] <module#>"
	    exit 1
	    ;;
	
	--)
	    shift; 
            break 
            ;;
    esac
done

eval $COMMAND

module unload nextflow
