#!/usr/bin/env Rscript
 
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(patchwork)
library(argparser)


library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Intronic/Exonic plot')
    p <- add_argument(p,'--rda', default = 'FALSE', help='add to RDS all raw files')
    p <- add_argument(p,'--sample_id', default = 'SAMPLE_ID', help='sample id')
    
    return(parse_args(p))
}

argv <- ParseArguments()

print(argv)

seurat_mark_doublets <- function(seurat_obj){
  ## input seurat object without information about singlet/doublet
  ## output seurat object with dblfinder column in meta.data with singlet/doublet
  
  ## do nothing if dblfinder column exists, return object back
  if ("dblfinder" %in% colnames(seurat_obj@meta.data)){
    return(seurat_obj)
  }
  if(nrow(seurat_obj@meta.data) < 20) {
    print(paste0("WARNING: only ",nrow(seurat_obj@meta.data)," cellbarcodes available"))
    print("DO NOT CALCULATING DOUBLETS. SET ALL CELLBARCODES AS SINGLETS")
    
    CB <- seurat_obj@meta.data$CB
    dblfinder <- data.frame(CB = CB, dblfinder = "singlet")
    rownames(dblfinder) <- CB
  } else {
    print("==>MARK DOUBLETS<==")
    ## find doublets
    sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA") %>%
      scDblFinder(.)
    
    dblfinder <- data.frame(
      CB = colnames(sce),
      dblfinder = colData(sce)$scDblFinder.class
    )
    
    rownames(dblfinder) <- dblfinder$CB
    dblfinder$CB <- NULL
  }
  
  new_seurat <- AddMetaData(seurat_obj, dblfinder)
  
  new_seurat
}

add_population <- function(sobj) {
  df <- sobj@meta.data %>%
    mutate(Intronic_Exonic_Population = case_when(summ2 == 4 | intronic_umi > exonic_umi ~ "I", 
                                                  .default = "II"),
           Intronic_Exonic_Population_v2 = case_when(intronic_umi > exonic_umi ~ "IA",
                                                     intronic_umi <= exonic_umi & summ2 == 4 ~ "IB",
                                                     .default = "II"))
  
  rownames(df) <- df$CB
  df <- select(df, -CB)
  
  AddMetaData(sobj, df)
}

# plot_qc <- function(sobj, sample_id) {
  
  
#   doublets <- table(sobj@meta.data$dblfinder)
#   sing <- pluck(doublets, "singlet", .default = 0)
#   dbl <- pluck(doublets, "doublet", .default = 0)
#   title <- str_glue("Singlets/Doublets for {sample_id},\n singlets = {sing}, doublets = {dbl}")
  
#   p1 <- ggplot(sobj@meta.data, aes(x = log10(exonic_umi), y = log10(intronic_umi))) +
#     geom_point(data = subset(sobj@meta.data, dblfinder == "singlet"), color = "gray", alpha = 0.7, size = 1) +
#     geom_point(data = subset(sobj@meta.data, dblfinder == "doublet"), color = "red3", alpha = 0.7, size = 1) +
#     geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", linewidth = 1)+
#     ggtitle(title)
  
#   pop <- table(sobj@meta.data$Intronic_Exonic_Population)
#   above <- pluck(pop, "I", .default = 0)  
#   below <- pluck(pop, "II", .default = 0) 
  
#   title <- str_glue("Selected population for {sample_id},\n I = {above}, II = {below}")
#   p2 <- ggplot(sobj@meta.data, aes(x = log10(exonic_umi), y = log10(intronic_umi))) +
#     geom_point(data = subset(sobj@meta.data, Intronic_Exonic_Population == "II"), color = "gray", alpha = 0.7, size = 1) +
#     geom_point(data = subset(sobj@meta.data, Intronic_Exonic_Population == "I"), color = "red3", alpha = 0.7, size = 1) +
#     geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", linewidth = 1)+
#     ggtitle(title)
  
#   df <- sobj@meta.data %>% filter(dblfinder == "singlet")
#   pop <- table(df$Intronic_Exonic_Population)
#   above <- pluck(pop, "I", .default = 0)  
#   below <- pluck(pop, "II", .default = 0) 
#   title <- str_glue("After doublet filtration for {sample_id},\n I = {above}, II = {below}")
#   p3 <- ggplot(df , aes(x = log10(exonic_umi), y = log10(intronic_umi))) +
#     geom_point(data = subset(df, Intronic_Exonic_Population == "II"), color = "gray", alpha = 0.7, size = 1) +
#     geom_point(data = subset(df, Intronic_Exonic_Population == "I"), color = "red3", alpha = 0.7, size = 1) +
#     geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", linewidth = 1)+
#     ggtitle(title)
  
  
#   return(p1+p2+p3+plot_layout(guides = "collect"))
# }

read_rds_v2 <- function(fpath){
  load(fpath)
  rm(genebody_raw)
  gc()
  
  only_filtered <- only_filtered %>% tibble::column_to_rownames(var = "CB") 
  only_filtered$CB <- NULL
  only_filtered$summ <- NULL
  withmono_raw <- AddMetaData(withmono_raw, only_filtered)
  withmono_raw
}

filter_low_quality <- function(obj){
  
  # First filter genes (min.cells = 3)
  counts <- GetAssayData(obj, layer = "counts")
  genes_to_keep <- rownames(counts)[rowSums(counts > 0) >= 3]
  
  # Then filter cells (min.features = 400) and genes together
  obj <- subset(obj, 
                features = genes_to_keep,
                subset = nFeature_RNA >= 400)
  obj
}

save_rds_file_v2 <- function(rda_path, sample_id) {
  new_sample_id <- ifelse(str_detect(sample_id,"_"), sample_id, str_replace(sample_id,"M","_M"))
  sobj <- read_rds_v2(rda_path)
  sobj <- filter_low_quality(sobj)
  sobj <- seurat_mark_doublets(sobj)
  sobj <- add_population(sobj)
  sobj$orig.ident <- new_sample_id
  Idents(sobj) <- "orig.ident"
  #write_csv(sobj@meta.data %>% mutate(sample_id = new_sample_id), str_c(new_sample_id, "_metadata_20260203.csv"))
  #plots <- plot_qc(sobj,new_sample_id)
  #ggsave(str_c(new_sample_id,"_qc.pdf"),plots, width = 18, height = 9)

  ## rename cellbarcodes from AAA..AAA-1 to AAA..AAA-1_<sample_id>
  new_barcodes <- paste0(colnames(sobj), "_", sobj@meta.data$orig.ident)
  sobj <- RenameCells(sobj, new.names = new_barcodes)

  saveRDS(sobj, str_c(new_sample_id,"_filtered.rds"))
  rm(sobj)
  gc()
}


save_rds_file_v2(argv$rda, argv$sample_id)
