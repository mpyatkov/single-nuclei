#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/**
 * Cell Ranger Reference Index Builder
 * 
 * This module checks if required reference indexes exist and builds them
 * if they haven't been created yet. Indexes are stored persistently.
 */

/**
 * Workflow: Check and build Cell Ranger reference indexes
 * 
 * @param indexes Comma-separated list of index names (e.g., "exonic_mm10,genebody_mm10")
 * @return Channel of built index directories
 * 
 * For each index name, expects:
 * - FASTA file: ${params.preprocessing.fasta_dir}/genome_<genome>.fa
 * - GTF file: ${params.preprocessing.gtfs_dir}/<index_name>.gtf
 */
workflow CHECK_DB {
    
    take:

    indexes
    
    main:

    // Parse comma-separated index list and create channel
    // Map each index name to its corresponding FASTA and GTF files
    indexes_ch=channel.of(indexes.split(','))
	.map{index_name ->
	    // Parse genome name from index (e.g., "exonic_mm10" -> "mm10")
	    def parts = index_name.split('_')
	    def genome = parts[1]
	    
	    // Return tuple: [index_name, fasta_path, gtf_path]
	    [index_name,
	     file("${params.fasta_dir}/genome_${genome}.fa", checkIfExists: true),
	     file("${params.gtfs_dir}/${index_name}.gtf", checkIfExists: true)
	    ]
	} | build_cellranger_index

    emit:
    build_cellranger_index.out.db_index
}

/**
 * Process: Build Cell Ranger reference index
 * 
 * Creates a Cell Ranger compatible reference genome from FASTA and GTF files.
 * Uses storeDir to persist indexes across pipeline runs.
 */
process build_cellranger_index {
    
    tag "${index_id}"

    cpus 16
    memory '64 GB'
    time '24h'
    
    // Load required modules - Cell Ranger for building reference indexes
    beforeScript "source \$HOME/.bashrc; module load ${params.modules.bcl2fastq}; module load ${params.modules.cellranger}"
    
    // Persist built indexes to shared location (not work directory)
    storeDir "${params.indexes_output_dir}"

    input:
    // Tuple: index_name, fasta_file, gtf_file
    tuple val(index_id), path(fasta), path(gtf)
    
    output:
    // Built index directory
    path("${index_id}"), emit: db_index
    
    script:

    """
    # Build Cell Ranger reference index
    cellranger mkref \
	   --genome=${index_id} \
	   --fasta=${fasta} \
	   --genes=${gtf} \
	   --nthreads=${task.cpus} \
	   --memgb=64
    """
}
