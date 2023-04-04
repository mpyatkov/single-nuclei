#!/usr/bin/env Rscript

## TMP: original location FindMarkersLoupe, mm10 converter and Utils.R: /projectnb/wax-dk/max/RSRC/G190

DEBUG <- FALSE

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Extraction/Combining')
    p <- add_argument(p,'--input_rds', default="input_rds.rds", help="load precalculated R object")
    p <- add_argument(p, '--recluster', help='if TRUE makes reclustering extracted cluster using "ncluster" and "umap_resolution" parameters', default="TRUE")
    p <- add_argument(p, '--nclusters', help='number of cluster on UMAP', default="12")
    p <- add_argument(p, '--umap_min_dist',
                      help='This controls how tightly the embedding is allowed compress points together. Large values - global structure. Smaller values - local structure (0.3 to 0.001)',
                      default="0.001")
    p <- add_argument(p, '--umap_npcs', help='Total Number of PCs to compute and store', default="8")
    
    p <- add_argument(p, '--umap_resolution', help='This parameters will overwrite "nclusters" parameter if is not "auto"', default="0.1,0.3")
    p <- add_argument(p, '--extract_combine_config', help='configuration file with a header. Contains info how to exctact/combine clusters in RDS', 
                      default="extract_combine_config.csv")
    
    p  <- add_argument(p, '--number_of_cores', 
                       help='some part of the script can be runned in parallel',
                       default=2)

    p <- add_argument(p,'--output_rds', default="output_rds.rds", help="output RDS filename")
    

    return(parse_args(p))
}

argv <- ParseArguments()

### LOAD LIBRARIES ###
library(tidyverse)
library(gridExtra)
library(grid)
library(Seurat)
library(harmony)
library(patchwork)
library(cowplot)
#source("Utils.R")
source("UmapPlot.R")
################
### SUPP FUNC ##

## parse arguments
## "1,2,3,4,5" -c(1,2,3,4,5)
str_to_vec <- function(s){
  str_split(s,",", simplify = T) %>% 
    str_trim() %>% 
    discard(~.=="") %>% 
    as.numeric()
}


################
#### PARAMETERS ####
################
## input RDS with Seurat object
arg_input_rds <- argv$input_rds

arg_output_rds <- argv$output_rds

recluster_bool <- ifelse(argv$recluster %in% c('TRUE', 'true'), TRUE, FALSE)

## maybe vector parameters
## does not work with ifelse
umap_resolution <- if(argv$umap_resolution == "auto"){"auto"} else {str_to_vec(argv$umap_resolution)}
number_of_clusters <- str_to_vec(argv$nclusters)
umap_min_dist <- str_to_vec(argv$umap_min_dist)
umap_npcs <- str_to_vec(argv$umap_npcs)

## configuration file
extract_combine_config_path <- argv$extract_combine_config

if (DEBUG) {
  # ONLY FOR DEBUG
  arg_input_rds <- "/projectnb2/wax-dk/max/SCexp/G190_G183/output/module_2_outputs/postaggregation/with-mono_matrix_with-mono_cellbarcodes/rds/module2_aggregated_seurat.rds"
  ## configuration file
  extract_combine_config_path <- "/projectnb2/wax-dk/max/SCexp/G190_G183/configuration/combine_extract_config.csv"
  umap_npcs <- c(8,10)
}

### MAIN ###

## only for debug (default_genes_list)
genes_list <- c('Cyp2e1', 'Glul', 'Oat', 'Gulo', 'Ass1', 'Hamp', 'Gstp1', 'Ubb', 'Selenbp2', 'Sds', 'Cyp2f2', 'Pck1',
                'Hal', 'Il2rb', 'Cxcr6', 'Gzma', 'Csf3r', 'S100a6', 'Pkhd1', 'Sox9', 'Epcam', 'Krt7', 'Krt19', 'Irf8',
                'Itgax', 'Clec4f', 'Csf1r', 'Jchain', 'Cd79a', 'Cd79b', 'Top2a', 'Stab2', 'Kdr', 'Aqp1', 'Fcgr2b',
                'Gpr182', 'Ebf1', 'Skap1', 'Ptprc', 'Ank3', 'Dcn', 'Colec11', 'Ecm1', 'Alb', 'Ttr', 'Apoa1', 'Serpina1c',
                'mt-Atp6', 'mt-Atp8', 'mt-Co3', 'mt-Co1', 'mt-Co2', 'mt-Nd2', 'mt-Nd4', 'mt-Cytb', 'mt-Nd1', 'mt-Nd3',
                'mt-Nd4l', 'mt-Nd5', 'mt-Nd6')

## load config file
## headers: c(ORIGINAL_CLUSTER_ID, COMBINING, EXTRACTING)
## TODO: check input format, otherwise produce error
extract_combine_config <- read_csv(extract_combine_config_path, col_names = T,
                                   show_col_types = FALSE,
                                   col_types = list(ORIGINAL_CLUSTER_ID = "c",
                                                    COMBINING = "c",
                                                    EXTRACTING = "c"))
# load rds file
# TODO: check if seurat file exists
input_seurat <- readRDS(arg_input_rds)

print(argv)

## function with side effect, produces error and terminates script if config is not correct
## if OK, return list
checking_configuration_file <- function(seurat_obj, config_file){

  ## make a table with extracting/combining
  tmp_clusters <- tibble(ORIGINAL_CLUSTER_ID = as.character(levels(Idents(seurat_obj)))) %>% 
    full_join(., config_file)  
  
  ## checking for the absence of unknown clusters in combine/extract config
  unknown_clusters <- setdiff(tmp_clusters$ORIGINAL_CLUSTER_ID, as.numeric(levels(seurat_obj))) 
  
  ## exit if unknown cluster is presented
  if (!is_empty(unknown_clusters)) {
    print(glue::glue("Cannot find the following clusters: '{unknown_clusters}' in Seurat object"))
    if (!DEBUG){ q("no", 1, FALSE) }
  }
  
  ## coalesce for combining clusters with ORIGINAL_CLUSTER_ID if NA. Set Exporting = 0 if NA
  tmp_clusters <- tmp_clusters %>% 
    rowwise() %>% 
    mutate(COMBINING = coalesce(COMBINING,ORIGINAL_CLUSTER_ID)) %>% 
    ungroup() %>% 
    replace_na(list(EXTRACTING = "0"))
  
  tmp_clusters
}

checked_config <- checking_configuration_file(input_seurat, extract_combine_config)


#### COMBINING


# if original_cluster_id == COMBINING skip step
# else return modified Seurat obj
if (any(checked_config$ORIGINAL_CLUSTER_ID != checked_config$COMBINING)) {
  print("COMBINING CLUSTERS...")
  
  ## vector values - new cluster names
  ## vector names - old cluster names
  new_cluster_names <- checked_config$COMBINING
  names(new_cluster_names) <- checked_config$ORIGINAL_CLUSTER_ID

  print("Old names/New names")
  new_cluster_names %>% print
  
  ## rename clusters
  input_seurat <- RenameIdents(input_seurat, new_cluster_names)
  
  ## need to rename seurat_clusters column in meta.data too
  input_seurat@meta.data$seurat_clusters <- Idents(input_seurat)
}


#### EXTRACTING


## checking for non-zero clusters for extraction
if (all(checked_config$EXTRACTING == 0)) {
  "There is no clusters for extraction in configuration file" %>% print
  if (!DEBUG){ q("no", 1, FALSE) }
}

clusters_to_extract <- checked_config %>% 
  filter(EXTRACTING == 1) %>% 
  pull(COMBINING) %>% 
  unique

input_seurat <- subset(input_seurat, subset = seurat_clusters %in% clusters_to_extract)


#### DRAW UMAP/UMAPS and recluster if recluster_bool set to TRUE

## if multiple_parameters is TRUE do not export meta.data and projections
multiple_parameters <- list(umap_resolution, number_of_clusters, umap_min_dist, umap_npcs) %>%
  map(length) %>%
  flatten_int() %>%
  discard(~.==1) %>%
  length() %>%
  map(`!=`,0) %>% unlist()

#### COMBINATIONS OF ALL PARAMETERS 
#### IF THEY ARE PRESENTED AS VECTORS
params <- purrr::cross(list(
  npcs = umap_npcs, 
  min.dist=umap_min_dist,
  number_of_clusters = number_of_clusters,
  resolution = umap_resolution)) %>%
  map(lift(list))


### add ix for each sublist starting from 001,002 .. length(params)
params <- map2(seq(length(params)), params,
               function(ix,l) list_modify(l,ix = str_pad(ix, 3, pad="0")))

print("Current params:")
print(params)

paste("Multiple params? ", multiple_parameters) %>% print


if (recluster_bool){
  
  if (multiple_parameters) {
    
    ## The output contains multiple UMAP plots, NO RDS
    print("Make multiple plot")
    
    library(foreach)
    library(doParallel)
    registerDoParallel(as.numeric(argv$number_of_cores))
    
    foreach(i = params) %dopar% {
    #foreach(i = params) %do% {
      
      ## RECLUSTERING
      print("reclustering")
      new_seurat <- recalc_umap_2(input_seurat,
                                  npcs = i$npcs,
                                  min.dist = i$min.dist,
                                  number_of_clusters = i$number_of_clusters,
                                  opt_resolution =  i$resolution,
                                  runsct = T)
      
      gc()

      #### PREPARING TITLE
      ## extract resolution parameter
      print("preparing title")
      
      print("Colnames")
      print(colnames(new_seurat@meta.data))
      
      resolution_text <- if (i$resolution == "auto") {
        grep("SCT_snn_res.", colnames(new_seurat@meta.data), value = T) %>%
          gsub("SCT_snn_res.","", .) %>% tail(., n=1)
      } else {
        i$resolution
      }
      
      TITLE <- str_interp("${i$ix} Aggregated. Resolution: ${resolution_text}, min.dist: ${i$min.dist}, #PCs: ${i$npcs}, #clusters: ${i$number_of_clusters}" )

      #### MAKING PLOT
      print("creating UmapDotplotTableHeatmap plot")
      final.umap.plot <- UmapDotplotTableHeatmapPlot(new_seurat, genes_list = genes_list, sample_id = "all", TITLE = TITLE)
      
      #### SAVE PLOT
      print("save to file")
      plot_name <- str_interp("${i$ix}_aggr_res${resolution_text}_pc${i$npcs}_dst${i$min.dist}_cl${i$number_of_clusters}.pdf")
      final.umap.plot %>% ggsave(filename = plot_name, width = 20, height = 15.50)
    } ## end of loop
    
  } else {
    
    print("Make individual plot")
    i <- flatten(params)
    new_seurat <- recalc_umap_2(input_seurat, 
                                npcs = i$npcs, 
                                min.dist = i$min.dist, 
                                number_of_clusters = i$number_of_clusters, 
                                opt_resolution =  i$resolution,
                                runsct = F)
    
    gc()
    
    #### PREPARING TITLE
    ## extract resolution parameter
    resolution_text <- grep("SCT_snn_res.", colnames(new_seurat@meta.data), value = T) %>%
      gsub("SCT_snn_res.","", .) %>% tail(., n=1)
    
    TITLE <- str_interp("Aggregated. Resolution: ${resolution_text}, min.dist: ${i$min.dist}, #PCs: ${i$npcs}, #clusters: ${i$number_of_clusters}" )
    
    #### MAKING PLOT
    all_samples_combined <- UmapDotplotTableHeatmapPlot(new_seurat, genes_list = genes_list, sample_id = "all", TITLE = TITLE)
    
    #### MAKING INDIVIDIAL PLOTS
    ob.list <- SplitObject(new_seurat, split.by = "sample_id")
    
    all_samples_by_sample_id <- lapply(X = names(ob.list), FUN = function(x) {
      tmp <- UmapDotplotTableHeatmapPlot(ob.list[[x]], 
                                         genes_list = genes_list, 
                                         sample_id = x, 
                                         TITLE = TITLE)
      cowplot::plot_grid(tmp)
    })
    
    #### COMBINE AND SAVE ALL PLOTS 
    append(list(cowplot::plot_grid(all_samples_combined)), all_samples_by_sample_id) %>%
    # all_samples_combined %>% 
      marrangeGrob(nrow =1, ncol=1) %>% 
      ggsave(filename = "aggr_umaps.pdf", width = 20.0, height = 15.50)
    gc()
    
    if (!DEBUG){
      saveRDS(new_seurat, arg_output_rds)  
    }
    
  }
  
} else {
  print("RECLUSTERING FALSE: Exporting Seurat object 'as is', ignoring UMAP specific parameters")
  
  #### PREPARING TITLE
  ## extract resolution parameter
  resolution_text <- grep("SCT_snn_res.", colnames(input_seurat@meta.data), value = T) %>%
    gsub("SCT_snn_res.","", .) %>% tail(., n=1)
  
  num_of_clusters <- Idents(input_seurat) %>% levels %>% length
  num_of_pc <- input_seurat@reductions$pca@stdev %>% length
  ## min.dist - UNKNOWN because Seurat object looks like does not store any distances
  TITLE <- str_interp("RECLUSTERING=FALSE. Resolution: ${resolution_text}, min.dist: Unknown, #PCs: ${num_of_pc}, #clusters: ${num_of_clusters}" )
  
  #### MAKING PLOT
  all_samples_combined <- UmapDotplotTableHeatmapPlot(input_seurat, genes_list = genes_list, sample_id = "all", TITLE = TITLE, rename_umap_labels = F)
  
  #### MAKING INDIVIDIAL PLOTS
  ob.list <- SplitObject(input_seurat, split.by = "sample_id")
  
  all_samples_by_sample_id <- lapply(X = names(ob.list), FUN = function(x) {
    tmp <- UmapDotplotTableHeatmapPlot(ob.list[[x]], 
                                       genes_list = genes_list, 
                                       sample_id = x, 
                                       TITLE = TITLE)
    cowplot::plot_grid(tmp)
  })
  
  #### COMBINE AND SAVE ALL PLOTS 
  append(list(cowplot::plot_grid(all_samples_combined)), all_samples_by_sample_id) %>% 
    marrangeGrob(nrow =1, ncol=1) %>% 
    ggsave(filename = "aggr_umaps.pdf", width = 20.0, height = 15.50)
  gc()
  
  if (!DEBUG){
    saveRDS(input_seurat, arg_output_rds)  
  }
}





