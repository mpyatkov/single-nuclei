# Single-Nuclei RNA-seq Pipeline - Module 1

A Nextflow pipeline for processing single-nuclei RNA-seq data by mapping FASTQ files against multiple GTF reference files with different gene counting strategies.

## Overview

This pipeline performs the initial processing of single-nuclei RNA-seq data by:
1. Building or reusing Cell Ranger reference indexes for 4 different GTF annotation strategies
2. Quantifying gene expression using Cell Ranger count
3. Creating Seurat RDS objects for downstream analysis
4. Generating aggregate QC metrics across all samples

### Four GTF Counting Strategies

The pipeline maps reads against 4 different gene annotation files to capture different aspects of gene expression:

| Strategy | Description | Use Case |
|----------|-------------|----------|
| **exonic** | Counts only exonic reads | Traditional bulk RNA-seq, cytoplasmic RNA |
| **genebody** | Counts exons + introns (pre-mRNA) | Single nuclei RNA-seq (default for most analyses) |
| **intronic-with-mono** | Counts intronic reads + monoexonic genes | Captures nascent transcription + monoexonic genes |
| **intronic-without-mono** | Counts only intronic reads (multiexonic genes only) | Nascent transcription only |

## Installation

### Prerequisites

- Access to an HPC cluster with SGE (Sun Grid Engine) scheduler
- Nextflow (version 21.10.6 or compatible)
- Cell Ranger (version 6.0.1)
- R (version 4.4.3) with packages: readxl, argparser
- Conda environment at `/projectnb2/wax-es/routines/condaenv/rlang4`

### Reference Data Requirements

The pipeline expects reference files in the following structure:
```
/projectnb/wax-es/routines/
├── FASTA/
│   └── genome_mm10.fa          # Reference genome FASTA
├── SC_GTFS/
│   ├── exonic_mm10.gtf         # Exonic annotations
│   ├── genebody_mm10.gtf       # Gene body annotations
│   ├── intronic-with-mono_mm10.gtf
│   └── intronic-without-mono_mm10.gtf
└── SC_INDEXES/                 # Built indexes (auto-generated)
```

## Usage

### 1. Prepare Sample Configuration

Create a `samples.xlsx` file with the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| `sample_id` | Unique sample identifier | `Sample_01` |
| `chemistry` | 10x chemistry version | `ARC-v1`, `SC3Pv3` |
| `condition` | Experimental condition/group | `Control`, `Treatment` |
| `path_to_r1` | Full path to R1 FASTQ file | `/path/to/Sample1_S1_L001_R1_001.fastq.gz` |

**Notes:**
- Multiple R1 files per sample are supported (one row per FASTQ file)
- FASTQ filenames must contain "_S" followed by the sample number (e.g., `Sample1_S1_L001_R1_001.fastq.gz`)
- The corresponding R2 file must be in the same directory

### 2. Run the Pipeline

```bash
./pipeline.sh
```

This will:
1. Convert `samples.xlsx` to `samples.csv`
2. Build reference indexes (if not already present)
3. Run Cell Ranger count for each sample against all 4 GTF files
4. Create RDS files and calculate aggregate metrics

## Output Structure

```
output/
├── raw_h5_and_cloupe_files/
│   └── {sample_id}/
│       ├── {sample_id}_exonic_mm10/
│       ├── {sample_id}_genebody_mm10/
│       ├── {sample_id}_intronic-with-mono_mm10/
│       └── {sample_id}_intronic-without-mono_mm10/
│           ├── cloupe.cloupe
│           ├── molecule_info.h5
│           └── web_summary.html
└── module_1_outputs/
    ├── summary_by_samples.csv
    └── rds/
        └── {sample_id}_output_rds.rds
```

### Key Output Files

- **molecule_info.h5**: UMI counts per gene/cell (input for aggregation)
- **cloupe.cloupe**: Loupe Browser visualization file
- **web_summary.html**: QC summary with sequencing metrics
- **{sample_id}_output_rds.rds**: Seurat object for downstream analysis
- **summary_by_samples.csv**: Aggregate QC metrics across all samples

## Configuration

Pipeline parameters are defined in `main.nf`:

```groovy
params {
    scc_project = 'wax-dk'                    // SGE project name
    output_dir = 'output'                      // Output directory
    
    preprocessing {
        indexes = 'exonic_mm10,intronic-without-mono_mm10,genebody_mm10,intronic-with-mono_mm10'
        main_db_path = '/projectnb/wax-es/routines'
        fasta_dir = "${params.preprocessing.main_db_path}/FASTA"
        gtfs_dir = "${params.preprocessing.main_db_path}/SC_GTFS"
        indexes_output_dir = "${params.preprocessing.main_db_path}/SC_INDEXES"
    }
    
    modules {
        cellranger = 'cellranger/6.0.1'
        bcl2fastq = 'bcl2fastq/2.20'
    }
}
```

### Customizing Resource Requirements

Process resources can be adjusted in `main.nf`:

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| `cellranger_count` | 16 | 112 GB | 24h |
| `build_cellranger_index` | 16 | 64 GB | 24h |
| `rds_and_h5_processing` | 8 | 32 GB | - |
| `calc_metrics` | Local executor | - | - |

## Troubleshooting

### Common Issues

**1. "Cannot find samples.csv configuration file"**
- Ensure `samples.xlsx` exists in the repository root
- Check that the XLSX file has the correct column names

**2. "File not found" errors for FASTA or GTF**
- Verify reference files exist at the expected paths
- Check that `params.preprocessing.main_db_path` is set correctly

**3. Cell Ranger out of memory**
- Increase `memory` directive in the `cellranger_count` process
- Check SGE cluster available resources

**4. FASTQ filename parsing errors**
- Ensure filenames contain "_S" pattern (e.g., `Sample_S1_L001_R1_001.fastq.gz`)
- Verify R1 and R2 files are in the same directory

## Pipeline Flow

```
samples.xlsx
    ↓
[Convert to CSV]
    ↓
samples.csv
    ↓
[Parse samples]
    ↓
Channel: [sample_id, chemistry, prefixes, fastq_dirs]
    ↓
[CHECK_DB] → Build/check reference indexes (4 GTF files)
    ↓
[cellranger_count] → Map against each index
    ↓
[rds_and_h5_processing] → Create RDS files
    ↓
[calc_metrics] → Aggregate QC metrics
    ↓
Results in output/
```

## Technical Details

### Workflow Steps

1. **CHECK_DB**: Verifies reference indexes exist, builds if needed
2. **cellranger_count**: Runs Cell Ranger for each sample × GTF combination
3. **rds_and_h5_processing**: Processes all 4 counting strategies, creates Seurat objects
4. **calc_metrics**: Aggregates web_summary.csv files from all samples

### Dependencies

- **Nextflow**: Workflow orchestration
- **Cell Ranger**: Read mapping and UMI counting
- **R packages**: Seurat, tidyverse, readxl (for XLSX conversion)
- **Conda**: R environment with required packages

## License

This pipeline is for internal research use at Waxman Lab, Boston University.

## Contact

For questions or issues, please contact the Waxman Lab or open an issue in the repository.
