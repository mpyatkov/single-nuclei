#!/usr/bin/env Rscript

## detects important markers for all clusters for each individual sample

DEBUG <- FALSE
DOTPLOT_NTOP_MARKERS <- 60     ## max number of markers for sample specific dotplots
SAMPLE_CLUSTER_NTOP <- 200     ## top N sumrank markers for each pair of sample_id - cluster_id
AVG_EXP_THR <- 0.2             ## threshold for DotPlot function
UNION_INTERSECTION_NTOP <- 30  ## top N markers for union/intersection dotplots from each sample


library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Detecting all markers')
  p <- add_argument(p,'--input_rds', default="input_rds.rds", help="load precalculated R object")
  return(parse_args(p))
}

argv <- ParseArguments()
print(argv)


## install lab libs
library(devtools)
devtools::install_github("mpyatkov/FindMarkersLoupe", upgrade = "never")
library(FindMarkersLoupe)


library(tidyverse)
library(stringi)
library(gridExtra)
library(grid)
library(Seurat)
library(harmony)
library(patchwork)
library(cowplot)
library(foreach)

if (DEBUG){
  source("/projectnb/wax-dk/max/SCexp/G190_G183/bin/UmapPlot.R")
  source("/projectnb/wax-dk/max/SCexp/35K_HEP_201022/findMarkersPairwise.R")
} else {
  source("UmapPlot.R")
  source("findMarkersPairwise.R")
}

#### SUPPLEMENTARY FUNCTIONS
## Wrapper for mp_dotplot with discrete gradient legend
plot_dotplot_for_any_cluster <- function(sobj, genes, title,...){
  #cols <- c("#fee08b","#ffffbf","#e6f598","#fdae61","#f46d43","#d53e4f","#abdda4","#66c2a5","#3288bd")
  #cols <- c("#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd")
  
  #cols <- c("#de2d26","#fc9272","#fee0d2","#deebf7","#9ecae1","#3182bd")
  cols <- c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061")
  
  mp_dotplot(sobj, genes = genes, title = title, scale = F) +
    #geom_point(aes(size=pct.exp), shape = 21, colour="gray", stroke=0.2)+ ## add edge stroke
    #scale_size(range = c(1,7))+                                           ## change default (1,6) sizes of the DotPlot circles
    scale_color_stepsn(colours = cols, nice.breaks = T,breaks = scales::breaks_pretty(n = 8))+
    theme(legend.key.width = unit(1.5, "cm"))
}

#### END FUNCTIONS

sobject_name <- "findmarkersloop_output"

if (DEBUG) {
  argv$input_rds <- "/projectnb/wax-dk/max/SCexp/35K_HEP_201022/hep_3samples_26K_recalculated.rds"
  sobject_name <- "TEST_TMP_REMOVE"
}

# load rds file
sobject <- readRDS(argv$input_rds)


## add sample_id - cluster_id column to sobject
sobject <- add_id_cluster(sobject)

## extract all possible combination of clusters
all_clusters <- pairwise_cluster_list(sobject)
## run FindMarkersLoupe
if (DEBUG) {
  res <- read_csv(str_glue("/projectnb/wax-dk/max/SCexp/35K_HEP_201022/ALL10_clusters_TEST_pairwise_clusters.csv"), col_names = T)
} else {
  res <- FindAllMarkersPairwise(seurat_obj = sobject, seurat_clusters = all_clusters)  
  write_csv(res, str_glue("{sobject_name}_pairwise_clusters.csv"))
}

## split sobject by samples
Idents(sobject) <- sobject@meta.data$seurat_clusters
sobject_samples <- SplitObject(sobject, split.by = "sample_id")
for (s in names(sobject_samples)) {
  sobject_samples[[s]] <- NormalizeData(sobject_samples[[s]], assay = "RNA")
}

######
###### Significant markers for each individual sample-cluster
######

## find preliminary list of markers for each individual sample
## to detect significance used 'sumrank'
## topgenes - top 20 genes for each pair sample_id - cluster_id
## it is possible that genes will be less than 20
prelim_list_final_markers <- map(names(sobject_samples), function(sample_id){
  temporary <- extract_significant(res, sample_id, topgenes = 200, adj_pval_thr = 0.001, log2fc_thr = 1)
  temporary$final
}) %>%
  set_names(names(sobject_samples))


## filter markers by avg.expression level
## last parameter (thresh_avgexp) threshold of scaled avg.exp calculated by Seurat::DotPlot function
filtered_markers <- generic_filtering_by_avgexp(sobject_samples, prelim_list_final_markers, thresh_avgexp = AVG_EXP_THR)

## drop markers which is not markers at all and has low avg.expr for specific cluster
final_res <- left_join(res, filtered_markers, by=c("gname","c1")) %>%
  drop_na()

## recalculate sumranks and extract significant markers again
filtered_adjusted_markers <- map(names(sobject_samples), function(sample_id){
  temporary <- extract_significant(final_res, sample_id, topgenes = SAMPLE_CLUSTER_NTOP)
  temporary$final
}) %>%
  set_names(names(sobject_samples))


## create markers plots for each individual sample
markers_plots <- map(names(sobject_samples), function(sample){
  
  current_clusters <- levels(Idents(sobject_samples[[sample]]))
  
  map(current_clusters, function(clst){
    
    title <- str_glue("Cluster {clst} ({sample})")
    
    # sort genes here
    filtered_genes <- filtered_adjusted_markers[[sample]] %>% filter(c1 == clst) %>% pull(gname) %>% head(DOTPLOT_NTOP_MARKERS) ## <--extract top N markers
    
    # arrange by avg.expr * expression percent
    sorted_genes <- DotPlot(sobject_samples[[sample]], features = filtered_genes)$data %>%
      filter(id == clst) %>%
      arrange(desc(avg.exp.scaled*pct.exp)) %>%
      pull(features.plot)
    
    plot_dotplot_for_any_cluster(sobject_samples[[sample]], genes = sorted_genes, title = title)})
}) %>% 
  set_names(., names(sobject_samples))

## save pdf for each sample (all clusters for a one sample)
map(names(markers_plots), function(sample){
  output_filename <- str_glue("{sample}_{DOTPLOT_NTOP_MARKERS}top_markers_individual_samples.pdf")
  markers_plots[[sample]] %>%
    map(function(p){plot_grid(wrap_plots(p))}) %>%
    marrangeGrob(nrow = 4, ncol = 1) %>%
    ggsave(filename = output_filename, width = 12, height = 18)
})


#####
##### Known biomarkers for each individual sample
#####

## create markers plots for each individual sample
knowm_markers_plots <- map(names(sobject_samples), function(sample){
  
  current_clusters <- levels(Idents(sobject_samples[[sample]]))
  
  title <- str_glue("{sample}")
  
  biomarkers <- list(
    PC = c("Glul","Gulo","Oat","Cyp2e1"),
    PP = c("Pck1","Cyp2f2","Hal"),
    Kupffer = c("Clec4f","Csf1r"),
    Immune = c("Ptprc"),
    HSC=c("Colec11","Dcn","Ecm1"),
    Endo=c("Stab2","Gpr182","Kdr","Fcgr2b","Aqp1"),
    Div=c("Top2a"),
    Cholang=c("Epcam","Krt19","Krt7","Sox9")
  )
  
  plot_dotplot_for_any_cluster(sobject_samples[[sample]], genes = biomarkers, title = title)
}) 


knowm_markers_plots %>%
  map(function(p){plot_grid(wrap_plots(p))}) %>%
  marrangeGrob(nrow = 4, ncol = 1) %>%
  ggsave(filename = "known_biomarkers_for_each_sample.pdf", width = 12, height = 18)


##### 
##### UNION/INTERSECT MARKERS FOR EACH CLUSTER BY INDIVIDUAL SAMPLES
##### 

Idents(sobject) <- sobject@meta.data$seurat_clusters

## create dotplots for union/intersection of top markers for each sample

union_intersect_dotplots <- map(levels(Idents(sobject)), function(cluster) {
  ## how many top genes to extract
  ntop <- UNION_INTERSECTION_NTOP
  
  ## extract data for specific cluster only
  tst_cells <- WhichCells(sobject, idents = cluster)
  tst_seu <- subset(sobject, cells = tst_cells) %>%
    NormalizeData(., assay = "RNA")
  
  ## TOPN markers for specific cluster
  mrks <- multiple_sample_markers(filtered_adjusted_markers, cluster, ntop = ntop)
  
  Idents(tst_seu) <- tst_seu@meta.data$id_cluster
  
  cols <- c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061")
  if (length(mrks$sig_intersect) == 0) {
    (mp_dotplot(tst_seu, genes = mrks$sig_union, scale = F)+
       scale_color_stepsn(colours = cols, nice.breaks = T))&
      plot_annotation(title = str_glue("No intersection of genes. Union of {ntop} most significant genes from each sample"))
  } else {
    (mp_dotplot(tst_seu, genes = mrks$sig_intersect, scale = F)+
       geom_point(aes(size=pct.exp), shape = 21, colour="gray", stroke=0.3)+
       theme(legend.key.width = unit(1.5, "cm"))+
       scale_color_stepsn(colours = cols, nice.breaks = T, breaks = scales::breaks_pretty(n = 8))) /
      (mp_dotplot(tst_seu, genes = mrks$sig_union, scale = F)+
         geom_point(aes(size=pct.exp), shape = 21, colour="gray", stroke=0.3)+
         theme(legend.key.width = unit(1.5, "cm"))+
         scale_color_stepsn(colours = cols, nice.breaks = T,breaks = scales::breaks_pretty(n = 8)))&
      plot_annotation(title = str_glue("Top intersection / bottom union of top {ntop} most significant genes from each sample"))
  }
  
  
  # if (length(mrks$sig_intersect) == 0) {
  #   (mp_dotplot(tst_seu, genes = mrks$sig_union)+
  #      geom_point(aes(size=pct.exp), shape = 21, colour="gray", stroke=0.3)+
  #      scale_size(range = c(1,7))+
  #      #scale_color_stepsn(colours = cols, nice.breaks = T,breaks = scales::breaks_pretty(n = 8))+
  #      theme(legend.key.width = unit(1.5, "cm"))+
  #      scale_color_stepsn(colours = cols, nice.breaks = T,breaks = scales::breaks_pretty(n = 8)))&
  #     plot_annotation(title = str_glue("Union of {ntop} most significant genes from each sample"))
  # } else {
  #   (mp_dotplot(tst_seu, genes = mrks$sig_intersect)+
  #      geom_point(aes(size=pct.exp), shape = 21, colour="gray", stroke=0.3)+
  #      scale_size(range = c(1,7))+
  #      theme(legend.key.width = unit(1.5, "cm"))+
  #      scale_color_stepsn(colours = cols, nice.breaks = T,breaks = scales::breaks_pretty(n = 8))) /
  #     (mp_dotplot(tst_seu, genes = mrks$sig_union)+
  #        geom_point(aes(size=pct.exp), shape = 21, colour="gray", stroke=0.3)+
  #        scale_size(range = c(1,7))+
  #        theme(legend.key.width = unit(1.5, "cm"))+
  #        scale_color_stepsn(colours = cols, nice.breaks = T,breaks = scales::breaks_pretty(n = 8)))&
  #     plot_annotation(title = str_glue("Top intersection / bottom union of top {ntop} most significant genes from each sample"))
  # }
  
})

## TODO: calculate nrow based on number of samples
union_intersect_dotplots %>%
  map(function(p){plot_grid(p)}) %>%
  marrangeGrob(nrow = 3, ncol = 1) %>%
  ggsave(filename = "significant_markers_all_samples_union_intersection.pdf", width = 12, height = 18)
