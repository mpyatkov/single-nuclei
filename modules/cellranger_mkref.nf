#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow CHECK_DB {
    
    take:

    indexes
    
    main:

    indexes_ch=channel.of(indexes.split(','))
	.map{it ->
	    def tmp = it.split('_')
	    def fastaix = tmp[1]
	    [it,
	     file("${params.preprocessing.fasta_dir}/genome_${fastaix}.fa", checkIfExists: true),
	     file("${params.preprocessing.gtfs_dir}/${it}.gtf", checkIfExists: true)
	    ]
	} | build_cellranger_index

    emit:
    build_cellranger_index.out.db_index
}

process build_cellranger_index {
    
    tag "${index_id}"

    cpus 16
    memory '64 GB'
    time '24h'
    beforeScript 'source $HOME/.bashrc; module load bcl2fastq/2.20; module load cellranger/3.0.1'
    
    storeDir "${params.preprocessing.indexes_output_dir}"

    input:
    tuple val(index_id), path(fasta), path(gtf)
    
    output:
    path("${index_id}"), emit: db_index
    
    script:

    """
    cellranger mkref \
	   --genome=${index_id} \
	   --fasta=${fasta} \
	   --genes=${gtf} \
	   --nthreads=${task.cpus} \
	   --memgb=64
    """
}
