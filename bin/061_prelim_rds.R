#!/usr/bin/env Rscript

DEBUG <- FALSE

library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Extract the user-defined cellbarcodes for all samples from aggregated.csv ')
  
  p <- add_argument(p, '--filtered', help="path to filtered matrix intronic-with-mono/genebody",
                    default='filtered_feature_bc_matrix.h5')
  
  p <- add_argument(p, '--userdefined_cellbarcodes', help='user defined cellbarcodes',
                    default="aggr-user-defined-barcodes.csv")

  p <- add_argument(p, '--mt_percent', 
                    help='mt percent (should be less than user-defined percent otherwise will no effect)',
                    default=10.0)
  
  p <- add_argument(p, '--min_genes', 
                    help='min number of genes detected in each cell',
                    default=400)
  
  p <- add_argument(p, '--min_counts', 
                    help='min number of molecules detected within each cell',
                    default=400)
  
  p <- add_argument(p,'--output_rds', default="module2_aggregated_seurat.rds", help="output RDS filename")
  
  return(parse_args(p))
}

argv <- ParseArguments()

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)

## read cellranger results from filtered dataset for intronic counting (without monoexonic genes)
read_h5_matrix <- function(matrix_path, matrix_name) {
  print(paste0("Loading: ", matrix_name," - " ,matrix_path))
  mx <- Read10X_h5(matrix_path)
  mx <- CreateSeuratObject(mx)
  mx@meta.data$samples <- matrix_name
  mx@meta.data$CB <- rownames(mx@meta.data)
  mx
}

recalculate_seurat_object_aggr <- function(
    seurat_object,
    meta_data,
    mt_percent,
    min_genes,
    min_counts,
    vars_to_regress = "percent.mt") 
{
  
  stats <- list()
  
  meta_data <- as.data.frame(meta_data)
  rownames(meta_data) <- meta_data$CB
  
  print("Number of rows in MetaData")
  print(nrow(meta_data)) 
  stats$userdefined_cb <- nrow(meta_data)
  
  print("==>RUN seurat object<==")
  new_seurat <- seurat_object[,meta_data$CB]
  
  print("==>RUN add meta data<==")
  new_seurat <- AddMetaData(new_seurat, meta_data) 
  
  print("==>RUN mt_percent, min_count, min_genes filtering<==")
  new_seurat[["percent.mt"]] <- PercentageFeatureSet(new_seurat, pattern = "^mt-")
  
  new_seurat <- subset(new_seurat, subset = percent.mt < mt_percent)
  stats$after1_mt_filter <- as.vector(table(Idents(new_seurat)))
  
  #new_seurat <- subset(new_seurat, subset = nFeature_RNA > min_genes & nCount_RNA > min_counts )
  new_seurat <- subset(new_seurat, subset = nFeature_RNA > min_genes )
  stats$after2_min_genes_filter <- as.vector(table(Idents(new_seurat)))
  
  new_seurat <- subset(new_seurat, subset = nCount_RNA > min_counts )
  stats$after3_min_counts_filter <- as.vector(table(Idents(new_seurat)))
  
  print("==>RUN SCTransform<==")
  new_seurat <- SCTransform(new_seurat, vars.to.regres = vars_to_regress, conserve.memory = T)
  gc()
  
  new_seurat <- RunPCA(new_seurat, npcs = 8) %>%
    RunHarmony(
      assay.use = "SCT",
      reduction = "pca",
      dims.use = 1:8,
      group.by.vars = "sample_id"
    ) %>%
    RunUMAP(reduction = "harmony", dims = 1:8, min.dist = 0.001) %>%
    FindNeighbors(reduction = "harmony", dims = 1:8)
  
  new_seurat <- FindClusters(new_seurat, resolution = 0.1) %>% 
    NormalizeData(., assay = "RNA")
  
  list(sobj = new_seurat,
       stats = stats)
}

################
#### PARAMETERS ####
################
arg_filtered_mx <- argv$filtered
arg_userdefined <- argv$userdefined_cellbarcodes

## base parameters
mt_percent <- argv$mt_percent
min_genes <- argv$min_genes
min_counts <- argv$min_counts

## maybe vector parameters

if (DEBUG) {
  arg_filtered_mx <- "/projectnb/wax-dk/max/globus/single_cell/single_nuclei/work/a6/18c969c9ded78ea63457265fe5aec2/filtered_feature_bc_matrix.h5"
  arg_userdefined <- "/projectnb/wax-dk/max/globus/single_cell/single_nuclei/work/e7/965eba18d228b2bfb1cc6bb73ebd09/with-mono_aggr-user-defined-barcodes.csv"
  mt_percent <- 10
  min_genes <- 400
  min_counts <- 400
}

## MAIN
userdefined_meta <- read_csv(arg_userdefined, col_names = T, show_col_types = FALSE) %>% 
  mutate(batch = sample_id)

userdefined_mx <- read_h5_matrix(arg_filtered_mx, "userdefined")


stats_sobj <- recalculate_seurat_object_aggr(seurat_object = userdefined_mx, 
                                              meta_data = userdefined_meta,
                                              mt_percent = mt_percent,
                                              min_genes = min_genes,
                                              min_counts = min_counts,
                                              vars_to_regress = c("percent.mt","nCount_RNA"))

print("STATS:")
print(stats_sobj$stats)
userdefined_mx_raw <- stats_sobj$sobj
gc()

if (!DEBUG){
  saveRDS(userdefined_mx_raw, argv$output_rds)  ## "module2_aggregated_seurat.rds"
}




