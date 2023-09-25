#!/usr/bin/env Rscript

DEBUG <- FALSE

library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Extract the user-defined cellbarcodes for all samples from aggregated.csv ')
  
  p <- add_argument(p,'--input_rds', default="input_rds.rds", help="load precalculated R object")
  
  p <- add_argument(p, '--nclusters', 
                    help='number of cluster on UMAP',
                    default="12")
  
  p <- add_argument(p, '--umap_min_dist', 
                    help='This controls how tightly the embedding is allowed compress points together. Large values - global structure. Smaller values - local structure (0.3 to 0.001)', default="0.01")
  
  p <- add_argument(p, '--umap_npcs', 
                    help='Total Number of PCs to compute and store',
                    default="10")
  
  p <- add_argument(p, '--umap_kparam', 
                    help='Total Number of PCs to compute and store',
                    default="20")
  
  p <- add_argument(p, '--umap_resolution', 
                    help='This parameters will overwrite "nclusters" parameter if is not "auto"',
                    default="auto")
  
  p  <- add_argument(p, '--number_of_cores', 
                     help='some part of the script can be runned in parallel',
                     default=2)
  
  p <- add_argument(p, '--dotplot_gene_list', 
                    help='multiple column csv with a header, each column is list of genes for dotplot, header - plot title',
                    default="dotplot_gene_list.csv")
  
  p <- add_argument(p, '--featureplot_gene_list', 
                    help='one column csv with a header, each column list of genes for dotplot',
                    default="featureplot_gene_list.csv")
  
  return(parse_args(p))
}

argv <- ParseArguments()

library(gridExtra)
library(grid)
library(Seurat)
library(harmony)
library(patchwork)
library(dplyr)
library(readr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(purrr)
library(gridExtra)
library(stringr)
source("UmapPlot.R")

## parse arguments
## "1,2,3,4,5" -c(1,2,3,4,5)
str_to_vec <- function(s){
  str_split(s,",", simplify = T) %>% 
    str_trim() %>% 
    discard(~.=="") %>% 
    as.numeric()
}


# ### extract clusters statistics from seurat_object
# extractClusterStats <- function(seurat_obj, labels = NULL) {
#   ## labels is data.frame with 2 columns (seurat_clusters, cluster_labels)
#   
#   cluster_stats <- seurat_obj@meta.data %>% 
#     select(sample_id, seurat_clusters) %>% 
#     mutate(all = n()) %>% 
#     
#     group_by(seurat_clusters) %>% 
#     mutate(n_cl = n()) %>% 
#     ungroup() %>% 
#     
#     group_by(sample_id) %>% 
#     mutate(n_sample = n()) %>% 
#     ungroup() %>% 
#     
#     group_by(seurat_clusters, sample_id) %>% 
#     summarise(ncells = n(), pt_per_sample=round(100*ncells/n_sample,2), pt_per_cluster = round(100*ncells/n_cl,2)) %>% 
#     distinct() %>% 
#     ungroup() 
#   
#   ## if labels exist
#   if (!labels %>% is.null()) {
#     cluster_stats <- left_join(cluster_stats, labels, by = "seurat_clusters") %>% 
#       relocate(cluster_labels, .after = seurat_clusters) %>%
#       rowwise() %>%
#       mutate(cluster_labels = paste0(cluster_labels,"(",seurat_clusters,")"))
#   }
#   
#   cluster_stats
#   
# }
# 
# ## rearrange columns for cluster_stats table
# prepareMetaData <- function(cluster_stats) {
#   spec <- build_wider_spec(cluster_stats, 
#                            names_from = sample_id, 
#                            values_from = c("ncells","pt_per_sample", "pt_per_cluster")) %>% 
#     arrange(sample_id, .name)
#   
#   pivot_wider_spec(cluster_stats, spec)
# }
# 

  # ## calculate cluster stats information from seurat object
  # cluster_stats <- extractClusterStats(seurat_obj, labels = df_cluster_names)
  #   ## export cluster stats (seurat_obj, cluster labels) for all clusters
  # prepareMetaData(cluster_stats) %>% write_csv("aggr_cluster_stats.csv")
  # 
  #  ## export meta-data (clusters, cellbarcodes, batch and other columns)
  # tmp_meta <- seurat_obj@meta.data %>% 
  #   select(CB,sample_id:percent.mt, -percent.mt, seurat_clusters) 
  # 
  # ## add information about user_defined clusters
  # if(!is.null(user_defined_labels)){
  #   tmp_meta <- tmp_meta %>% 
  #     left_join(., user_defined_labels$new_labels %>% tibble::enframe(name="seurat_clusters", value = "cluster_labels", by="seurat_clusters"))
  # }
  # tmp_meta %>% write_csv("aggr_meta_data.csv")
  # 
  # #### projections export ####
  # seurat_obj@reductions$umap@cell.embeddings %>% 
  #   as.data.frame(.) %>% 
  #   tibble::rownames_to_column(var = "Barcode") %>% 
  #   write_csv("aggr_projections.csv")
  
  

ExportDotplotTables <- function(dotplot_df){
  
  ## raw expression without any scaling
  dt_raw_expression <- dotplot_df %>% 
    select(gene = features.plot, avg_expression = avg.exp, cluster_id = id) %>% 
    pivot_wider(names_from = cluster_id,
                values_from = avg_expression)
  
  ## scaled expression by maximum expression for specific gene 
  dt_scaled_expression <- dotplot_df %>% 
    select(gene = features.plot, avg_expression = avg.exp, cluster_id = id) %>% 
    group_by(gene) %>% 
    mutate(max_in_col = max(avg_expression),
           scaled_expression = avg_expression/max_in_col) %>% 
    ungroup() %>% 
    select(gene, cluster_id,scaled_expression) %>% 
    pivot_wider(names_from = cluster_id,
                values_from = scaled_expression)
  
  ## raw percent of how many cell in cluster expressed specific gene
  dt_raw_percent_of_cells <- dotplot_df %>% 
    select(gene = features.plot, pct_cells_expressed = pct.exp, cluster_id = id) %>% 
    pivot_wider(names_from = cluster_id,
                values_from = pct_cells_expressed)
  
  return(list(
    list(
      name = "raw_expression",
      data = dt_raw_expression
    ), 
    list(
      name = "scaled_expression",
      data = dt_scaled_expression
    ),
    list(
      name = "percent_of_cells",
      data = dt_raw_percent_of_cells
    )
  ))
}


################
#### PARAMETERS ####
################

## files with additional info
arg_input_rds <- argv$input_rds
dotplot_gene_list_arg <- argv$dotplot_gene_list
feature_gene_list_arg <- argv$featureplot_gene_list
# seurat_cluster_labels <- argv$seurat_cluster_labels

#labels_for_umap <- ifelse(argv$labels_for_umap %in% c('Y', 'y'), TRUE, FALSE)

## maybe vector parameters
umap_resolution <- ifelse(argv$umap_resolution == "auto", "auto", str_to_vec(argv$umap_resolution))
number_of_clusters <- str_to_vec(argv$nclusters)
umap_kparam <- str_to_vec(argv$umap_kparam)
umap_min_dist <- str_to_vec(argv$umap_min_dist)
umap_npcs <- str_to_vec(argv$umap_npcs)
parallel <- TRUE

if (DEBUG) {
  arg_input_rds <- "/projectnb/wax-dk/max/ATAC_TCPO/G183_G193_SC/output/module_2_outputs/original_results_before_hybrid_cellbarcodes/postaggregation/with-mono_matrix_with-mono_cellbarcodes/rds/module2_aggregated_seurat.rds"

  ## files with additional info
  dotplot_gene_list_arg <- "/projectnb/wax-dk/max/globus/single_cell/single_nuclei/configuration/dotplot_gene_list.csv"
  feature_gene_list_arg <- "/projectnb/wax-dk/max/globus/single_cell/single_nuclei/configuration/featureplot_gene_list.csv"

  ## maybe vector parameters
  umap_resolution <- "auto"
  number_of_clusters <- 8
  umap_kparam <- 20
  umap_min_dist <- 0.001
  umap_npcs <- 8
  parallel <- FALSE
}

########################
#### MAIN ####
########################

userdefined_mx_raw <- readRDS(arg_input_rds)

## ------start reading configs-------
## read gene list and remove NA because in most case the gene list
## will not be the same by size (default_genes)
full_genes_list <- read_csv(dotplot_gene_list_arg, col_names = T, show_col_types = FALSE) %>% 
  as.list() %>% 
  map(discard, is.na)

## read genes for featureplots
feature_genes_list <- read_csv(feature_gene_list_arg, col_names = T) %>% pull(genes)

## ------end reading configs-------

## if multiple_parameters is TRUE do not export meta.data and projections
multiple_parameters <- list(umap_resolution, number_of_clusters, umap_kparam, umap_min_dist, umap_npcs) %>%
  map(length) %>%
  flatten_int() %>%
  discard(~.==1) %>%
  length() %>%
  map(`!=`,0) %>% unlist()

paste("Multiple params? ", multiple_parameters) %>% print

## read gene list and remove NA because in most case the gene list
## will not be the same by size (default_genes)
full_genes_list <- read_csv(dotplot_gene_list_arg, col_names = T, show_col_types = FALSE) %>% 
  as.list() %>% 
  map(discard, is.na)

## genes list for individual dotplots
genes_by_groups_wo_default <- names(full_genes_list) %>% 
  discard(.=="default_genes") %>% 
  full_genes_list[.]


#### COMBINATIONS OF ALL PARAMETERS 
#### IF THEY ARE PRESENTED AS VECTORS
params <- purrr::cross(list(
  npcs = umap_npcs, 
  min.dist=umap_min_dist,
  number_of_clusters = number_of_clusters,
  k.param = umap_kparam,
  resolution = umap_resolution)) %>%
  map(lift(list))


### add ix for each sublist starting from 001,002 .. length(params)
params <- map2(seq(length(params)), params,
               function(ix,l) list_modify(l,ix = str_pad(ix, 3, pad="0")))

if (multiple_parameters) {
  print("Make multiple plot")
  
  library(foreach)
  library(doParallel)
  
  `%uniform_do%` <- `%do%`
  if (parallel) {
    `%uniform_do%` <- `%dopar%`
  }
  
  registerDoParallel(as.numeric(argv$number_of_cores))
  
  foreach(i = params) %uniform_do% {

    print("reclustering")
    new_seurat <- recalc_umap_2(userdefined_mx_raw,
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
    
    ## "Export projections for temporary plots")
    as_tibble(new_seurat@reductions$umap@cell.embeddings, rownames = "CB") %>%
      left_join(., new_seurat@meta.data %>% select(CB, seurat_clusters, sample_id), by = "CB") %>%
      write_csv(., str_glue("{i$ix}_aggr_res{resolution_text}_pc{i$npcs}_dst{i$min.dist}_cl{i$number_of_clusters}.csv"), col_names = T)
    
    TITLE <- str_interp("${i$ix} Aggregated. Resolution: ${resolution_text}, min.dist: ${i$min.dist}, #PCs: ${i$npcs}, #clusters: ${i$number_of_clusters}" )
    
    #### MAKING PLOT
    print("creating UmapDotplotTableHeatmap plot")
    final.umap.plot <- UmapDotplotTableHeatmapPlot_v2(seurat_obj = new_seurat, 
                                                      genes_list = DEFAULT_MARKER_LIST,
                                                      sample_id = "all",
                                                      additional_genes_list = CELL_PAPER_MARKERS, 
                                                      rename_umap_labels = T,
                                                      TITLE = TITLE)
    
    #### SAVE PLOT
    print("save to file")
    plot_name <- str_interp("${i$ix}_aggr_res${resolution_text}_pc${i$npcs}_dst${i$min.dist}_cl${i$number_of_clusters}.pdf")
    final.umap.plot %>% ggsave(filename = plot_name, width = 23.0, height = 17.05)
  }
}else {
  print("Make individual plot")
  
  i <- flatten(params)
  
  ## I rewrite original userdefined_mx_raw because the raw object 
  ## is not required for downstream analysis
  new_seurat <- recalc_umap_2(userdefined_mx_raw, 
                              npcs = i$npcs, 
                              min.dist = i$min.dist, 
                              number_of_clusters = i$number_of_clusters, 
                              opt_resolution =  i$resolution,
                              runsct = T)
  gc()
  
  #### PREPARING TITLE
  ## extract resolution parameter
  resolution_text <- grep("SCT_snn_res.", colnames(new_seurat@meta.data), value = T) %>%
    gsub("SCT_snn_res.","", .) %>% tail(., n=1)
  
  TITLE <- str_interp("Aggregated. Resolution: ${resolution_text}, min.dist: ${i$min.dist}, #PCs: ${i$npcs}, #clusters: ${i$number_of_clusters}" )
  
  all_samples_combined <- UmapDotplotTableHeatmapPlot_v2(seurat_obj = new_seurat, 
                                                         genes_list = DEFAULT_MARKER_LIST,
                                                         sample_id = "all",
                                                         additional_genes_list = CELL_PAPER_MARKERS, 
                                                         rename_umap_labels = T,
                                                         TITLE = TITLE)
  ## export projections for aggregated plot
  as_tibble(new_seurat@reductions$umap@cell.embeddings, rownames = "CB") %>%
    left_join(., new_seurat@meta.data %>% select(CB, seurat_clusters, sample_id), by = "CB") %>%
    write_csv(., str_glue("aggregated_projections.csv"), col_names = T)
  
  
  #### MAKING INDIVIDIAL PLOTS
  ob.list <- SplitObject(new_seurat, split.by = "sample_id")
  
  all_samples_by_sample_id <- lapply(X = names(ob.list), FUN = function(x) {
    
    tmp <- UmapDotplotTableHeatmapPlot_v2(seurat_obj = ob.list[[x]], 
                                          genes_list = DEFAULT_MARKER_LIST,
                                          sample_id = x,
                                          additional_genes_list = CELL_PAPER_MARKERS, 
                                          rename_umap_labels = T,
                                          TITLE = TITLE)
    
    ## export projections for separate samples
    ## as_tibble(ob.list[[x]]@reductions$umap@cell.embeddings, rownames = "CB") %>%
    ##     left_join(., ob.list[[x]]@meta.data %>% select(CB, seurat_clusters), by = "CB") %>%
    ##     write_csv(., str_glue("{x}_projections.csv"), col_names = T)
    
    cowplot::plot_grid(tmp)
  })
  
  #### COMBINE AND SAVE ALL PLOTS 
  append(list(cowplot::plot_grid(all_samples_combined)), all_samples_by_sample_id) %>%
    # all_samples_combined %>% 
    marrangeGrob(nrow =1, ncol=1) %>% 
    ggsave(filename = "aggr_umaps.pdf", width = 23.0, height = 17.05)

  ####
  ## individual dotplots
  ####
  print("Export individual dotplots")
  names(genes_by_groups_wo_default) %>% 
    map(function(name) {
      genes <- genes_by_groups_wo_default[[name]]
      ngenes <- length(genes)
      mp_dotplot(new_seurat, genes, name)
    }) %>% 
    marrangeGrob(nrow = 4, ncol = 1) %>%
    ggsave(filename = "aggr_dotplots_by_groups.pdf", width = 13.9, height = 18)
  
  ####
  ## individual featureplots
  ####
  print("Export Featureplots")
  DefaultAssay(new_seurat) <- "RNA"
  feature_genes_list %>% 
    map(function(gname) {
      FeaturePlot(new_seurat, features = gname, pt.size = 1, order = T, raster = T)
    }) %>% 
    marrangeGrob(nrow = 2, ncol = 2, layout_matrix = matrix(1:4, 2, 2, TRUE)) %>%
    #grid.arrange(grobs = ., ncol = 2, as.table = FALSE) %>% 
    ggsave(filename = "aggr_featureplots.pdf", width = 11, height = 8.5)
  
  ####
  ## Export dotplot as tables for all possible dotplot lists
  ## For all aggregated samples together
  print("Export DotPlots as tables")
  names(full_genes_list) %>%
    map(function(name){
      genes <- full_genes_list[[name]]
      ngenes <- length(genes)
      tmpdtp <-DotPlot(new_seurat, features = genes, assay = "RNA")
      dt_tables <- ExportDotplotTables(tmpdtp$data)
      walk(dt_tables, function(l) {
        write_csv(l$data, file = paste0("all_samples_",name,"_",l$name,".csv"), col_names = T)
      })
    })
  
  ## TODO: export dotplot for each separate sample
  ## for each sample in all.split.by.sample
  ##      export 3 tables for sample
  ####
  
  if (!DEBUG){
    saveRDS(new_seurat, "module2_postprocessed_seurat.rds")  
  }
  
}

fn <- "Rplots.pdf"
if (file.exists(fn)) {file.remove(fn)}


