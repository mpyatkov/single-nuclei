#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

/**
 * Single-Nuclei RNA-seq Pipeline - Module 1
 * 
 * This pipeline maps FASTQ files against 4 different GTF reference files
 * to quantify gene expression using different counting strategies:
 * - exonic: Count only exonic reads
 * - genebody: Count pre-mRNA (exons + introns) for single nuclei
 * - intronic-with-mono: Count intronic reads + monoexonic gene reads
 * - intronic-without-mono: Count only intronic reads for multiexonic genes
 */

params {
    // SCC project identifier for job submission
    scc_project = 'wax-dk'
    
    // Output directory for all results
    output_dir = 'output'

    // Reference genome and annotation configuration
    preprocessing {
	// Four GTF counting strategies (comma-separated)
	// These define which genomic features to count for each reference
	indexes = 'exonic_mm10,intronic-without-mono_mm10,genebody_mm10,intronic-with-mono_mm10'
	
	// Base path for reference databases
	main_db_path = '/projectnb/wax-es/routines'
	
	// Directory containing reference genome FASTA files
	fasta_dir = "${params.preprocessing.main_db_path}/FASTA"
	
	// Directory containing GTF annotation files
	gtfs_dir = "${params.preprocessing.main_db_path}/SC_GTFS"
	
	// Directory for storing built Cell Ranger indexes
	indexes_output_dir = "${params.preprocessing.main_db_path}/SC_INDEXES"
    }

    // Sample configuration
    module1 {
	// samples.csv is auto-generated from samples.xlsx by pipeline.sh
	// Expected columns: sample_id, chemistry, condition, path_to_r1
	samples_general = 'samples.csv'
    }
    
    // Software module versions
    modules {
	cellranger = 'cellranger/6.0.1'
	bcl2fastq = 'bcl2fastq/2.20'
    }
}

/**
 * Extract sample prefix and directory from R1 FASTQ file paths
 * 
 * Example: "Sample1_S1_L001_R1_001.fastq.gz" -> prefix="Sample1", dir="/path/to/fastq"
 * Assumes filenames contain "_S" followed by sample number
 * 
 * @param path_to_r1 List of R1 FASTQ file paths
 * @return Tuple of [prefixes, fastq_directories] as lists
 */
def vget_prefix(path_to_r1) {
    
    // Extract sample prefix from filename (text before "_S")
    prefixes = path_to_r1.collect{it -> {
	fname = new File(it.toString()).getName()
	fname[0..fname.indexOf("_S")-1]
    }}

    // Extract parent directory containing FASTQ files
    fastqdirs = path_to_r1.collect{it -> {
	new File(it.toString()).parent
    }}
    return [prefixes, fastqdirs]
}

// Import workflow for checking/building Cell Ranger reference indexes
include { CHECK_DB } from './modules/cellranger_mkref.nf'

// Load sample configuration as a value channel
samples_ch = channel.value(file(params.module1.samples_general))

/**
 * Create channel with sample metadata from CSV file
 * 
 * Input CSV format (with header):
 *   sample_id, chemistry, condition, path_to_r1
 * 
 * Output channel format:
 *   [sample_id, chemistry, comma_separated_prefixes, comma_separated_fastq_dirs]
 * 
 * Multiple R1 files per sample are grouped together
 */
samples_reads_ch = Channel
    .fromPath(params.module1.samples_general)
    .ifEmpty{exit 1, "Cannot find ${params.module1.samples_general} configuration file"}
    .splitCsv(skip:1)  // Skip header row
    .groupTuple(by:0)   // Group by sample_id (column 0)
    .map{it -> [it[0],              // sample_id
		it[1][0],           // chemistry (take first occurrence)
		vget_prefix(it[-1])[0].join(","), // comma-separated prefixes
		vget_prefix(it[-1])[1].join(",") // comma-separated fastq directories
	]}

/**
 * Main workflow entry point
 * 
 * Pipeline flow:
 * 1. CHECK_DB: Verify/build Cell Ranger reference indexes for all 4 GTFs
 * 2. cellranger_count: Map FASTQ files against each reference index
 * 3. rds_and_h5_processing: Create RDS files from count matrices
 * 4. calc_metrics: Aggregate summary statistics
 */
workflow {
    
    // Step 1: Get or build reference indexes for all GTF files
    db_indexes = CHECK_DB(params.preprocessing.indexes)

    // Step 2: Run Cell Ranger count for each sample against each reference
    samples_reads_ch.combine(db_indexes) | cellranger_count

    // Step 3: Process output to extract sample metadata
    // Input: [sample_id, genome_index_name, h5_path, sample_directory]
    // Groups multiple GTF outputs per sample, selects first genome index as representative
    samples_dirs_ch = cellranger_count.out.h5
	.map{tuple -> [
	    tuple[0],                                    // sample_id
	    tuple[1].toString().split("_").last(),       // genome name (last part after underscore)
	    tuple[1].toString(),                         // full h5 output path
	    tuple[2].parent                              // parent directory
	]}
	.groupTuple(by:0)  // Group all outputs by sample_id
	.map{sample_id, genome_ixs, gtf_names, sample_paths -> [
	    sample_id,           // sample identifier
	    genome_ixs[0],       // representative genome index (first one)
	    sample_paths[0]      // path to sample output directory
	]}

    // Step 4: Create RDS files and extract cell barcodes
    rds_and_h5_processing(samples_dirs_ch)

    rds_and_h5_processing.out.rds | final_rds
    // Step 5: Calculate aggregate metrics from all samples
    summary_path = channel.value("${projectDir}/${params.output_dir}/raw_h5_and_cloupe_files/")
    calc_metrics(summary_path)
}

/**
 * Process: Run Cell Ranger count to quantify gene expression
 * 
 * Maps FASTQ files against a reference genome and generates:
 * - molecule_info.h5: UMI counts per gene/cell
 * - cloupe.cloupe: Loupe browser file
 * - web_summary.html: QC summary report
 */
process cellranger_count {
    
    tag "${sample_index}"
    cpus 16
    memory '112 GB'
    time '24h'
    
    // Load required software modules
    beforeScript "source \$HOME/.bashrc; module load ${params.modules.bcl2fastq}; module load ${params.modules.cellranger}"
    
    // Store results permanently (not in work directory)
    storeDir "${params.output_dir}/raw_h5_and_cloupe_files/${sample_id}"

    input:
    // Tuple: sample_id, chemistry, fastq_prefixes, fastq_dirs, reference_index
    tuple val(sample_id), val(chemistry), val(fastq_prefix), val(fastqdir), val(index_id)
        
    output:
    // Tuple: sample_id, genome_index_name, output_directory
    tuple val(sample_id), val(genome_ix), path("${sample_index}"), emit: h5
    
    script:
    
    // Extract genome name from index path
    genome_ix = index_id.getName()
    // Create unique sample+index identifier
    sample_index="${sample_id}_${genome_ix}"
    
    """
    # Run Cell Ranger count
    cellranger count \
           --id=${sample_index} \
           --sample=${fastq_prefix} \
           --fastqs=${fastqdir} \
           --transcriptome=${index_id} \
           --localcores=${task.cpus} \
           --chemistry=${chemistry} \
           --disable-ui

    # Clean up unnecessary large files to save space
    pushd \$PWD/${sample_index}/outs
    rm -rf *.bam *.bai raw_feature_bc_matrix filtered_feature_bc_matrix analysis
    popd
    
    # Remove temporary Cell Ranger files
    pushd ${sample_index}
    find . -mindepth 1 -name "*" | grep -v outs | xargs rm -rf
    popd    
      
    # Move outputs from outs/ to top-level directory
    mv \$PWD/${sample_index}/outs/* \$PWD/${sample_index}/
    rm -rf \$PWD/${sample_index}/outs
    """
}

/**
 * Process: Calculate aggregate metrics from all sample summaries
 * 
 * Collects web_summary.csv files from all samples and generates
 * a combined summary report with key QC statistics
 */
process calc_metrics {
    // Use local executor - lightweight task
    executor 'local'
    
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'
    
    publishDir path: "${params.output_dir}/module_1_outputs/", mode: "copy", pattern: "*.csv", overwrite: true
    

    input:
    val(summary_files_location)

    output:
    path("summary_by_samples.csv")
    

    script:

    """
    00_calculate_metrics.R --summary_path ${summary_files_location} --output summary_by_samples.csv
    """
}

/**
 * Process: Create RDS files from count matrices and extract cell barcodes
 * 
 * Reads filtered feature matrices from all 4 GTF strategies and:
 * - Creates Seurat objects
 * - Extracts union of cell barcodes across all strategies
 * - Outputs RDS file for downstream analysis
 */
process rds_and_h5_processing {
    tag "${sample_id}"
    
    beforeScript 'source $HOME/.bashrc; module load miniconda'
    conda '/projectnb2/wax-es/routines/condaenv/rlang4'
    
    cpus 8
    memory '32 GB'

    publishDir path: "${params.output_dir}/module_1_outputs/rds/", mode: "copy", pattern: "*.rds", overwrite: true, saveAs : {filename -> "${sample_id}_${filename}"}

    input:
    // Tuple: sample_id, representative_genome_index, sample_output_path
    tuple val(sample_id), val(genome_ix), val(sample_path)
    
    output:
    // Cell barcodes CSV and RDS file
    tuple val("${sample_id}"), path("${sample_id}_union_cellbarcodes.csv"), emit: cellbarcodes
    tuple val("${sample_id}"), path("output_rds.rds"), emit: rds

    script:

    """
    # Process all 4 GTF counting strategies
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

process final_rds {
  tag "${sample_id}"

    beforeScript 'source $HOME/.bashrc; module load miniconda'
    
    cpus 1
    memory '16 GB'

    publishDir path: "${params.output_dir}/module_1_outputs/final_rds/", mode: "copy", pattern: "*.rds", overwrite: true 

    input:
    tuple val(sample_id), path(rda)
    
    output:
    tuple val("${sample_id}"), path("*filtered.rds"), emit: cellbarcodes

    script:
    """
    02_rda_to_rds.R --rda ${rda} --sample_id ${sample_id}
    """
}
