#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

// TODO: check some parameters first
include { CHECK_DB } from './modules/cellranger_mkref.nf'

// samples_config_ch = channel.value(file(params.sample_options))
samples_ch = channel.value(file(params.module1.samples_general))

// channel with samples reads (sample_id, fastqdir)
samples_reads_ch = Channel
    .fromPath(params.module1.samples_general)
    .splitCsv(skip:1)
    .groupTuple(by:0)
    .map{it -> [it[0],              //sample_id
		it[1][0],           //chemistry
		it[-2].join(","),   //fastq_prefix
		it[-1].join(",")]}  //fastqdir

userdefined_params = Channel
    .fromPath(params.module1.sample_options)
    .splitCsv(skip: 1)

// MODULE 1
workflow module_1 {
    
    // GET index or create if it is not created before and put it to
    // the storeDir for permanent storing
    db_indexes = CHECK_DB(params.preprocessing.indexes)

    // run cellranger count
    // samples_reads_ch | cellranger_count
    samples_reads_ch.combine(db_indexes) | cellranger_count

    // TODO: make process which will collect all metrics after counting
    // in one file, sources here:
    // projectnb2/wax-dk/max/SCexp/calculate_metrics.R
    
    // Return the sample_dir paths for each sample
    samples_dirs_ch = cellranger_count.out.h5
	//.view()
	// .collect()
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

    // modifying molecule h5 files (sample_id, path to directory with
    // molecule info h5 file) only for genebody and intronic-with-mono GTF files
    molecule_h5_dirs = cellranger_count.out.h5
	.map{it->[it[0],it[-1]]}
	.filter{it -> it[1] =~ "genebody|with-mono"}

    // union of cellbarcodes for each sample id
    cellbarcodes = rds_and_h5_processing.out.cellbarcodes
    cellbarcodes.cross(molecule_h5_dirs)
	.map{it -> it.flatten()}
        // sample_id, cellbarcodes, output_name_prefix, path_to_molecule_info
	.map{it -> [it[0], it[1], it[3].getName(), it[3]]}  | modifying_molecule_h5
    
    // recalculate h5 matrix and cloupe files
    samples_dirs_ch.join(cellbarcodes) | cellranger_reanalyse

    rds = rds_and_h5_processing.out.rds
    userdefined_ch = rds.join(userdefined_params)

    // make plots for module I
    intronic_exonic_plot(userdefined_ch, params.module1.seurat_nclusters, params.module1.downstream_gtf)
}

// module1: calculate each sample separately
process cellranger_count {
    
    tag "${sample_index}"
    // echo true
    cpus 16
    memory '112 GB'
    time '24h'
    beforeScript 'source $HOME/.bashrc; module load bcl2fastq/2.20; module load cellranger/6.0.1'
    clusterOptions "-P ${params.scc_project} -l scratch_free=300G"

    
    storeDir "${params.output_dir}/raw_h5_and_cloupe_files/${sample_id}"

    input:
    
    tuple val(sample_id), val(chemistry), val(fastq_prefix), val(fastqdir), val(index_id)
        
    output:

    tuple val(sample_id), val(genome_ix), path("${sample_index}"), emit: h5
    // tuple val(sample_id), path("molecule_info.h5"), emit: h5_orig
    
    script:
    
    // ex: G183M1_exonic_mm10
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

    ## TODO --chemistry=threeprime --chemistry=SC3Pv3

    ## remove all files in outs directory which are not necessary now
    pushd \$PWD/${sample_index}/outs
    rm -rf *.bam *.bai raw_feature_bc_matrix filtered_feature_bc_matrix
    
    ## TODO: remove the following line after checking that modified
    ##  file works well in aggregation
    cp molecule_info.h5 molecule_info.h5.orig
    popd
    
    ## remove all cellranger temporary files excluding outs directory
    pushd ${sample_index}
    find . -mindepth 1 -name "*" | grep -v outs | xargs rm -rf
    popd    
    
    ## move all necessary files to a level up and remove outs directory  
    mv \$PWD/${sample_index}/outs/* \$PWD/${sample_index}/
    rm -rf \$PWD/${sample_index}/outs
    """
}


process rds_and_h5_processing {
    tag "${sample_id}"
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'
    
    cpus 4
    memory '32 GB'

    publishDir path: "${params.output_dir}/module_1_outputs/rds/", mode: "copy", pattern: "*.rds", overwrite: true, saveAs : {filename -> "${sample_id}_${filename}"}
    //publishDir path: "${params.output_dir}/module_1_outputs/${sample_id}/data/", mode: "copy", pattern: "${sample_id}_${gtfname_genome}_info.h5", overwrite: true

    input:
    tuple val(sample_id), val(genome_ix), val(sample_path)
    
    output:
    //tuple val("${sample_id}"), path("${sample_id}_${gtfname_genome}_info.h5"), emit: modified_h5
    tuple val("${sample_id}"), path("${sample_id}_union_cellbarcodes.csv"), emit: cellbarcodes
    tuple val("${sample_id}"), path("output_rds.rds"), emit: rds

    script:

    """
    ## pdfs with predefined and user defined parameters
    ## modifying molecule_info.h5 file for intronic+monoexonic genes

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

process modifying_molecule_h5 {

    executor 'local'
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'

    publishDir path: "${params.output_dir}/corrected_h5_and_cloupe_files/${sample_id}/${output_name}/", mode: "copy", pattern: "*.h5", overwrite: true
    
    input:
    tuple val(sample_id), val(cellbarcodes), val(output_name), val(molecule_info_dir_path)
    
    output:
    path("molecule_info.h5")

    script:
  
    """
    cp ${molecule_info_dir_path}/molecule_info.h5 ./

    03_modify_molecule_h5.R --molecule_info_h5 molecule_info.h5  --cellbarcodes ${cellbarcodes}
    """
}


process cellranger_reanalyse {
    // reanalize each raw h5 file with selected subset of cellbarcodes
    // The output contains corrected h5 matrix (not molecule h5 file) which can be used for
    // downstream analysis and  corrected cloupe for each individual sample
    // corrected - means that cloupe file will contain cellbarcodes from 4 GTFs
    
    tag "${sample_id}"
    cpus 4
    memory '32 GB'
    time '4h'
    
    beforeScript 'source $HOME/.bashrc; module load bcl2fastq/2.20; module load cellranger/6.1.2'


    publishDir path: "${params.output_dir}/corrected_h5_and_cloupe_files/${sample_id}", mode: "copy", overwrite: true
    // storeDir "${params.output_dir}/corrected_h5_and_cloupe_files/${sample_id}"

    input:
    
    tuple val(sample_id), val(genome_ix), val(sample_path), val(barcodes)
        
    output:

    path("${sample_id}_exonic_${genome_ix}")
    path("${sample_id}_intronic-with-mono_${genome_ix}")
    path("${sample_id}_intronic-without-mono_${genome_ix}")
    path("${sample_id}_genebody_${genome_ix}")
	 
    script:
    genebody = "${sample_id}_genebody_${genome_ix}"
    with_mono = "${sample_id}_intronic-with-mono_${genome_ix}"
    without_mono = "${sample_id}_intronic-without-mono_${genome_ix}"
    exonic = "${sample_id}_exonic_${genome_ix}"
    
    """
    cellranger reanalyze --id=exonic \
           --matrix=${sample_path}/${exonic}/raw_feature_bc_matrix.h5 \
	   --localcores=${task.cpus} \
	   --barcodes=${barcodes}

    cellranger reanalyze --id=with-mono \
           --matrix=${sample_path}/${with_mono}/raw_feature_bc_matrix.h5 \
	   --localcores=${task.cpus} \
	   --barcodes=${barcodes}

    cellranger reanalyze --id=without-mono\
           --matrix=${sample_path}/${without_mono}/raw_feature_bc_matrix.h5 \
	   --localcores=${task.cpus} \
	   --barcodes=${barcodes}

    cellranger reanalyze --id=genebody \
           --matrix=${sample_path}/${genebody}/raw_feature_bc_matrix.h5 \
	   --localcores=${task.cpus} \
	   --barcodes=${barcodes}

    mkdir -p ${exonic}
    cp exonic/outs/{cloupe.cloupe,filtered_feature_bc_matrix.h5} ${exonic}/

    mkdir -p ${without_mono}
    cp without-mono/outs/{cloupe.cloupe,filtered_feature_bc_matrix.h5} ${without_mono}/

    mkdir -p ${with_mono}
    cp with-mono/outs/{cloupe.cloupe,filtered_feature_bc_matrix.h5} ${with_mono}/

    mkdir -p ${genebody}
    cp genebody/outs/{cloupe.cloupe,filtered_feature_bc_matrix.h5} ${genebody}/

    ## find . -mindepth 1 -name "*" | grep -v ${sample_id} | xargs rm -rf
    """
}


process intronic_exonic_plot {
    tag "${sample_id}"
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'

    time '4h'
    cpus 4
    memory '32 GB'
    // executor 'local'

    publishDir path: "${params.output_dir}/module_1_outputs/${sample_id}/plots/", mode: "copy", pattern: "*.pdf", overwrite: true
    publishDir path: "${params.output_dir}/module_1_outputs/${sample_id}/data/", mode: "copy", pattern: "*.csv", overwrite: true
    publishDir path: "${params.output_dir}/module_1_outputs/${sample_id}/data/", mode: "copy", pattern: "${sample_id}_${gtfname_genome}_info.h5", overwrite: true

    input:
    tuple val(sample_id), val(rds), val(slope), val(intercept), val(mt_gtf), val(mt_percent), val(population)
    val(number_of_umap_clusters)
    val(downstream_gtf)

    output:
    path("*.pdf")
    path("*.csv")
    
    script:
    """
    02_intronic_exonic_plots.R \
      --sample_id ${sample_id}\
      --input_rds ${rds}\
      --number_of_umap_clusters ${number_of_umap_clusters}\
      --user_slope ${slope}\
      --user_intercept ${intercept}\
      --user_mt_gtf ${mt_gtf}\
      --user_mt_percent ${mt_percent}\
      --user_population ${population}\
      --downstream_gtf ${downstream_gtf}
    """
}

// END MODULE 1


// MODULE 2 PRECONFIGURATION

workflow module2_preconfiguration {
    // this subprocess is only required because users can modify the
    // aggregation.csv files and remove some of the samples
    
    // get all modified h5 files from all directories in output
    channel.fromPath("${params.output_dir}/corrected_h5_and_cloupe_files/**molecule_info.h5")
	.filter{it -> it =~ params.module1.downstream_gtf }
	.map{it -> [it.parent.parent.getName(),
		    it]}
	.toSortedList{a,b -> a[0] <=> b[0]}
	.map{it->it.collect{lst -> lst.join(",")}}
    | h5_modified_file

    // append information from samples_ch to aggregation file
    module_2_aggr_ch = create_aggregation_file(samples_ch, h5_modified_file.out.h5_list)
}

process get_user_defined_cellbarcodes {

    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'

    executor 'local'
    
    publishDir path: "${params.output_dir}/module_2_outputs", mode: "copy", pattern: "*.csv", overwrite: true
    
    input:
    path(aggr_file)
    path(userdefined_paths)
    
    output:
    path("aggr-user-defined-barcodes.csv")

    script:
    """
    04_get_user_defined_cellbarcodes.R --aggregation ${aggr_file} \
                                    --userdefined_meta ${userdefined_paths} \
                                    --output aggr-user-defined-barcodes.csv
    """
}

process create_aggregation_file {

    executor 'local'
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'
    // echo true
    // we use publish dir here because we will start module_2 as as
    // separate workflow
    publishDir path: "${projectDir}/configuration/", mode: "copy", overwrite: true
    
    input:
    val(samples)
    val(modified_h5)
    
    output:
    path("aggregation.csv"), emit: aggr_csv

    script:
    """
    05_combine_configs.R --samples ${samples} --modified_h5 ${modified_h5} --output aggregation.csv
    """
}

process h5_modified_file {

    executor 'local'
    input:
    val(path_lines)

    output:
    path("h5_modified_list.csv"), emit: h5_list

    script:
    def cmd="echo 'sample_id,molecule_h5' > h5_modified_list.csv\n"
    for( int i=0; i<path_lines.size(); i++ ) {
	cmd+="echo '${path_lines[i]}' >> h5_modified_list.csv\n"
    } 
    cmd
}

// END MODULE 2 PRECONFIGURATION


// MODULE 2

workflow module_2 {

    // here I am using aggregation.csv from the output directory, because
    // this file may have been edited, so the version from the working
    // directory is not suitable for us
    module_2_aggr_csv_ch = channel.value(projectDir+'/configuration/aggregation.csv')
    
    // create temporary file with all pathes to individual user-defined meta files
    module_2_meta_ch = channel.fromPath("${params.output_dir}/module_1_outputs/**user-defined_cloupe*.csv")
	.map{it -> it.toString()}
	.collectFile(name: 'samples_path.csv', newLine: true)

    // Problem: After aggregation all samples represented in aggregation.csv
    // will have specific numeration based on row position in this file.
    // ex. G183M1 --> cellbarcodes-1
    //     G183M2 --> cellbarcodes-2
    //     G183M4 --> cellbarcodes-3
       
    // it supposed that we can remove some lines from aggregation.csv
    // (removing low quality samples), that means we should reorder
    // 'dashes' (-1,-2,... for each sample) in user defined meta files
    // to associate cellbarcodes in right order

    // Combine all files with user-defined meta data to one file
    // and choose only cellbarcodes from samples which have only been
    // presented in aggregation.csv, changing 'dash' to correct one.
    // output: table with 1 column - CB (correct "dashed" cell barcodes)
    // copy it to /module_2_outputs/aggr-user-defined-barcodes.csv
    // because probably we would like to change it in feature
    get_user_defined_cellbarcodes(module_2_aggr_csv_ch, module_2_meta_ch)
    
    // userdefined barcodes which were collected from individual
    // samples barcodes and sliced by userdefined parameters (see configuration/sample_options.csv)
    user_barcodes_ch = channel.value(file(params.output_dir+'/module_2_outputs/aggr-user-defined-barcodes.csv'))

    // get h5 filtered matrix from aggregated samples
    module_2_aggr_csv_ch | cellranger_aggregate
    aggr_filtered_ch = cellranger_aggregate.out.filtered_matrix
    
    // only for debug
    // aggr_filtered_ch = channel.value(file("${params.output_dir}/module_2_outputs/aggregation/count/filtered_feature_bc_matrix.h5"))
    
    // postprocessing using matrix with the filtered barcodes and
    // table with user defined barcodes
    aggregation_postprocessing(aggr_filtered_ch, user_barcodes_ch)
    // aggregation_postprocessing(cellranger_aggregate.out.filtered_matrix, user_barcodes_ch)
}


process cellranger_aggregate {
    // cache false // for test only
    // echo true
    cpus 16
    // penv 'smp'
    memory '64 GB'

    time '24h'
    beforeScript 'source $HOME/.bashrc; module load bcl2fastq/2.20; module load cellranger/6.0.1'

    publishDir path: "${params.output_dir}/module_2_outputs/aggregation", mode: "copy", overwrite: true

    input:
    val(aggr_file)
    
    output:
    path("*")
    path("filtered_feature_bc_matrix.h5"), emit: filtered_matrix
    
    """
    cellranger aggr --id=aggregated --csv=${aggr_file}

    ## remove all files in outs directory which are not necessary now
    pushd aggregated/outs/count
    rm -rf raw_feature_bc_matrix filtered_feature_bc_matrix raw_feature_bc_matrix.h5
    popd

    ## remove all cellranger temporary files excluding outs directory
    find . -mindepth 1 -name "*" | grep -v aggregated | xargs rm -rf
    
    ## move all necessary files to a level up and remove outs directory  
    mv \$PWD/aggregated/outs/* \$PWD/
    mv \$PWD/count/filtered_feature_bc_matrix.h5 \$PWD/
    mv \$PWD/count/cloupe.cloupe \$PWD/
    rm aggregation.csv
    rm -rf \$PWD/aggregated/ 
    """
}

process aggregation_postprocessing {

    tag "aggr_postprocessing"
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'
    
    cpus 8
    memory '64 GB'

    publishDir path: "${params.output_dir}/module_2_outputs/postaggregation/plots", mode: "copy", pattern: "*.pdf", overwrite: true
    publishDir path: "${params.output_dir}/module_2_outputs/postaggregation/data", mode: "copy", pattern: "*.csv", overwrite: true
    publishDir path: "${params.output_dir}/module_2_outputs/postaggregation/rds", mode: "copy", pattern: "*.rds", overwrite: true
    
    input:
    val(filtered_mx)
    val(barcodes_mx)
    
    output:
    path("*.pdf")
    path("*.rds"), optional: true
    path("*.csv"), optional: true
    
    script:
    """
    06_postaggregation.R --filtered ${filtered_mx} \
		 --userdefined_cellbarcodes ${barcodes_mx}\
		 --nclusters ${params.module2.seurat_nclusters}\
		 --mt_percent ${params.module2.mt_percent}\
		 --min_genes ${params.module2.min_genes}\
		 --min_counts ${params.module2.min_counts}\
		 --dotplot_gene_list ${projectDir}/${params.module2.dotplot_gene_list}\
		 --featureplot_gene_list ${projectDir}/${params.module2.featureplot_gene_list}\
                 --seurat_cluster_labels ${projectDir}/${params.module2.seurat_cluster_labels}\
                 --labels_for_umap ${params.module2.labels_for_umap}\
                 --umap_min_dist ${params.module2.umap_min_dist}\
                 --umap_npcs ${params.module2.umap_npcs} \
                 --umap_kparam ${params.module2.umap_kparam} \
                 --umap_resolution ${params.module2.umap_resolution} \
                 --number_of_cores ${task.cpus}
    
    ## TODO: replace SCC installed imagemagick to conda
    ## for range of parameters (like --nclusters "6,7,8") represent all plots as 
    ## combined tiled plot

    if [[ -n \$(find . -name "*_aggr*pdf") ]]; then
       montage -geometry +0+0 -density 150 -tile 2x2 \$(find . -name "*_aggr*pdf" | sort -n | paste -s -d ' ') combined_plot.pdf

       rm -rf *_aggr*pdf
    fi
    """
}
