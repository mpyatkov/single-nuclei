#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

params {
    scc_project = 'wax-dk'
    output_dir = 'output'

    preprocessing {
	indexes = 'exonic_mm10,intronic-without-mono_mm10,genebody_mm10,intronic-with-mono_mm10'
	main_db_path = '/projectnb/wax-es/routines'
	fasta_dir = "${params.preprocessing.main_db_path}/FASTA"
	gtfs_dir = "${params.preprocessing.main_db_path}/SC_GTFS"
	indexes_output_dir = "${params.preprocessing.main_db_path}/SC_INDEXES"
    }

    module1 {
	samples_general = 'samples.csv'
    }
}

def vget_prefix(path_to_r1) {
    
    prefixes = path_to_r1.collect{it -> {
	fname = new File(it.toString()).getName()
	fname[0..fname.indexOf("_S")-1]
    }}

    fastqdirs = path_to_r1.collect{it -> {
	new File(it.toString()).parent
    }}
    return [prefixes, fastqdirs]
}

include { CHECK_DB } from './modules/cellranger_mkref.nf'

samples_ch = channel.value(file(params.module1.samples_general))

samples_reads_ch = Channel
    .fromPath(params.module1.samples_general)
    .ifEmpty{exit 1, "Cannot find ${params.module1.samples_general} configuration file"}
    .splitCsv(skip:1)
    .groupTuple(by:0)
    .map{it -> [it[0],              //sample_id
		it[1][0],           //chemistry
		vget_prefix(it[-1])[0].join(","), // prefixes
		vget_prefix(it[-1])[1].join(",") // fastqdirs
	]}

workflow {
    
    db_indexes = CHECK_DB(params.preprocessing.indexes)

    samples_reads_ch.combine(db_indexes) | cellranger_count

    samples_dirs_ch = cellranger_count.out.h5
	.map{it -> [it[0],
		    it[1].toString().split("_").last().toString(),
		    it[1].toString(),
		    it[2].parent]}
	.groupTuple(by:0)
	.map{sample_id,
	     genome_ixs,
	     gtf_names,
	     sample_paths -> [sample_id,
			      genome_ixs[0],
			      sample_paths[0]]}

    rds_and_h5_processing(samples_dirs_ch)

    summary_path = channel.value("${projectDir}/${params.output_dir}/raw_h5_and_cloupe_files/")
    calc_metrics(summary_path)
}

process cellranger_count {
    
    tag "${sample_index}"
    cpus 16
    memory '112 GB'
    time '24h'
    beforeScript 'source $HOME/.bashrc; module load bcl2fastq/2.20; module load cellranger/6.0.1'
    storeDir "${params.output_dir}/raw_h5_and_cloupe_files/${sample_id}"

    input:
    tuple val(sample_id), val(chemistry), val(fastq_prefix), val(fastqdir), val(index_id)
        
    output:
    tuple val(sample_id), val(genome_ix), path("${sample_index}"), emit: h5
    
    script:
    
    genome_ix = index_id.getName()
    sample_index="${sample_id}_${genome_ix}"
    
    """
    cellranger count \
           --id=${sample_index} \
           --sample=${fastq_prefix} \
           --fastqs=${fastqdir} \
           --transcriptome=${index_id} \
           --localcores=${task.cpus} \
           --chemistry=${chemistry} \
           --disable-ui

    pushd \$PWD/${sample_index}/outs
    rm -rf *.bam *.bai raw_feature_bc_matrix filtered_feature_bc_matrix analysis
    popd
    
    pushd ${sample_index}
    find . -mindepth 1 -name "*" | grep -v outs | xargs rm -rf
    popd    
      
    mv \$PWD/${sample_index}/outs/* \$PWD/${sample_index}/
    rm -rf \$PWD/${sample_index}/outs
    """
}

process calc_metrics {
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'
    publishDir path: "${params.output_dir}/module_1_outputs/", mode: "copy", pattern: "*.csv", overwrite: true
    
    executor 'local'

    input:
    val(summary_files_location)

    output:
    path("summary_by_samples.csv")
    

    script:

    """
    00_calculate_metrics.R --summary_path ${summary_files_location} --output summary_by_samples.csv
    """
}

process rds_and_h5_processing {
    tag "${sample_id}"
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'
    
    cpus 8
    memory '32 GB'

    publishDir path: "${params.output_dir}/module_1_outputs/rds/", mode: "copy", pattern: "*.rds", overwrite: true, saveAs : {filename -> "${sample_id}_${filename}"}

    input:
    tuple val(sample_id), val(genome_ix), val(sample_path)
    
    output:
    tuple val("${sample_id}"), path("${sample_id}_union_cellbarcodes.csv"), emit: cellbarcodes
    tuple val("${sample_id}"), path("output_rds.rds"), emit: rds

    script:

    """
    01_rds_and_molecule_h5_preparation.R \
      --intronic_h5_filtered ${sample_path}/${sample_id}_intronic-without-mono_${genome_ix}/filtered_feature_bc_matrix.h5 \
      --intronic_h5_raw ${sample_path}/${sample_id}_intronic-without-mono_${genome_ix}/raw_feature_bc_matrix.h5 \
      --exonic_h5_filtered ${sample_path}/${sample_id}_exonic_${genome_ix}/filtered_feature_bc_matrix.h5 \
      --exonic_h5_raw ${sample_path}/${sample_id}_exonic_${genome_ix}/raw_feature_bc_matrix.h5 \
      --genebody_h5_filtered ${sample_path}/${sample_id}_genebody_${genome_ix}/filtered_feature_bc_matrix.h5 \
      --genebody_h5_raw ${sample_path}/${sample_id}_genebody_${genome_ix}/raw_feature_bc_matrix.h5 \
      --intronic_withmono_h5_filtered ${sample_path}/${sample_id}_intronic-with-mono_${genome_ix}/filtered_feature_bc_matrix.h5 \
      --intronic_withmono_h5_raw ${sample_path}/${sample_id}_intronic-with-mono_${genome_ix}/raw_feature_bc_matrix.h5 \
      --output_rds "output_rds.rds" \
      --output_cellbarcodes "${sample_id}_union_cellbarcodes.csv"
    """

}
