### creates combined plot with layout (umap/table/dotplot)
combinedUmapTableDotplot <- function(seurat_obj, TITLE, labels_for_umap = NULL, genes_list, cluster_stats, sample_id = "all") {
  
  #### table for plots
  table_for_plots <- cluster_stats %>% 
    select(seurat_clusters, sample_id, ncells)
  
  ## dotplot for all samples together
  DefaultAssay(seurat_obj) <- "RNA"
  gc()
  
  
  if (!is.null(labels_for_umap)) {
    ## rename clusters (set new_labels)
    seurat_obj <- RenameIdents(seurat_obj, labels_for_umap$new_labels)
  }
  
  ## Plot all samples together on one UMAP without labels
  aggr_umap_all_clusters <- mp_plot_umap(seurat_obj, title = TITLE, draw_labels = T)
  
  if (!is.null(labels_for_umap)) {
    ## rename clusters back (set old_labels)
    seurat_obj <- RenameIdents(seurat_obj, labels_for_umap$old_labels)
  }
  
  ## TODO: at the moment all genes will be represented on the dotplots
  aggr_umap_all_clusters_dotplot <- mp_dotplot(seurat_obj, genes_list)
  
  
  clusters_stats_table_plot <- plot_table(get_table_by_id(sample_id, table_for_plots), title = "")
  
  ## set layout for 11 plots
  layout <- "
    AAAAC
    AAAAC
    AAAAC
    BBBBB
    "
  combined_all_clusters <- cowplot::plot_grid((aggr_umap_all_clusters+theme(plot.margin=margin(0,1,0,0,"cm")))+
                                                aggr_umap_all_clusters_dotplot+
                                                clusters_stats_table_plot+
                                                plot_layout(design = layout)+plot_annotation(title = sample_id))
  combined_all_clusters
}

#' Supplementary function to create UMAP plots
#'
#' @param seurat_obj - Seurat object with predefined clusters
#' @param draw_labels - bool, if TRUE draw labels
#' @param title - Add title to plot
#' @param ... - Additional parameters
#' @return A patchworked ggplot object if combine = TRUE; otherwise, a list of ggplot objects
#' @export
#'
#' @examples
mp_plot_umap <- function(seurat_obj, title = "", draw_labels = T,...){
  tmp <- DimPlot(seurat_obj, reduction = "umap", ...) + 
    ggtitle(title)
  
  ## If labels is required add them to umap plot
  if (draw_labels) {
    
    ## to make labels the same color as legend colors use the following code
    # aggr_umap_all_clusters <- LabelClusters(..., color = unique(ggplot_build(aggr_umap_all_clusters)$data[[1]]$colour), ...)
    
    tmp <- LabelClusters(tmp, 
                         id = "ident", 
                         size = 5, 
                         repel = T,  
                         box.padding = 1, 
                         segment.size = 0.2,
                         max.overlaps = 10)
  }
  tmp
}

## supplementary function to create dotplot
mp_dotplot <- function(sobj, genes, title = "", ...) {
  DotPlot(sobj, features = genes, assay = "RNA")+
    RotatedAxis()+ 
    {if(title !="") ggtitle(title)} +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal", 
          legend.justification = "center")
}

# ## supplementary function to create dotplot
# mp_dotplot <- function(sobj, genes, title = "", ...) {
#   DotPlot(sobj, features = genes, assay = "RNA", ...)+
#     RotatedAxis()+ 
#     {if(title !="") ggtitle(title)} +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank())
# }

#' Create table with title as grob object
#'
#' @param df - any data.frame
#' @param title - title
#'
#' @return grob object
#' @export
#'
#' @examples
mp_plot_table <- function(df, title){
  table <- gridExtra::tableGrob(df, rows = NULL, theme = ttheme_minimal(base_size = 12))
  h <- grobHeight(table)
  w <- grobWidth(table)
  title <- textGrob(title, y=unit(0.5,"npc") + 1.5*h, vjust=0, gp=gpar(fontsize=20))
  table <- gTree(children = gList(title, table))
  table    
}


### simplified version of previous function. 
### Using it only in case we have multiple options for some parameters
### example 'pc = "6,10,20"'
### generates multiple pdfs with umap/table/dotplot for each combination of parameters

# print(npcs)
# #seurat_obj <- userdefined_mx_raw
# system.time(
#   # function(new_seurat, npcs, var.genes.no, harmony_vars = "batch", min.dist, number_of_clusters) {
#   seurat_obj <- recalc_umap(seurat_obj, 
#                             npcs = npcs, 
#                             min.dist = min.dist, 
#                             number_of_clusters = number_of_clusters, 
#                             opt_resolution = resolution,
#                             k.param = k.param)
#   
# )


# #### scaling and normalization for RNA assay ####
# seurat_obj <- seurat_obj %>% 
#   NormalizeData(., assay = "RNA") 

# 
# gc()
# 
# 





UmapDotplotTableHeatmapPlot <- function(seurat_obj, 
                                 genes_list,
                                 counter = "000", 
                                 sample_id = "all", 
                                 TITLE = "") {
  
  #### PREDICT CELL TYPES
  ### Assign predicted cell types labels
  ## get new_labels,old_labels,heatmap_table
  celltypes_list <- FindCellTypesByMarkers(sobj = seurat_obj)
  celltypes_heatmap <- CreateCellTypesHeatmap(celltypes_list$heatmap_table)
  
  # ## dotplot for all samples together
  # DefaultAssay(seurat_obj) <- "RNA"

  #### ADD LABELS FOR UMAP IF EXISTS
  ## rename clusters (set new_labels)
  seurat_obj <- RenameIdents(seurat_obj, celltypes_list$new_labels)

    ## Plot all samples together on one UMAP without labels
  aggr_umap_all_clusters <- mp_plot_umap(seurat_obj, title = TITLE, draw_labels = T)
  
  ## rename clusters back (set old_labels)
  seurat_obj <- RenameIdents(seurat_obj, celltypes_list$old_labels)
  
  
  #### CREATE DOTPLOT
  ## TODO: at the moment all genes will be represented on the dotplots
  aggr_umap_all_clusters_dotplot <- mp_dotplot(seurat_obj, genes_list)
  
  #### CREATE STATS TABLE WITH CLUSTERS
  clusters_stats_table_plot <- mp_plot_table(list_sample_cluster_stats(seurat_obj)[[sample_id]], title = "")
  
  layout <- "
    AAAAC
    AAAAC
    AAAAC
    BBBBB
    "
  all_samples_combined <- cowplot::plot_grid((aggr_umap_all_clusters+theme(plot.margin=margin(0,1,0,0,"cm")))+
                                                aggr_umap_all_clusters_dotplot+
                                                clusters_stats_table_plot+
                                                plot_layout(design = layout)+plot_annotation(title = sample_id))
  
  #### ATTACH HEATMAP TO UmapStatsDotplot plot
  umap_heatmap_layout <- "
AAAAA####
AAAAA#BBB
AAAAA#BBB
AAAAA####
"
  umap_heatmap_plot <- all_samples_combined+celltypes_heatmap+plot_layout(design = umap_heatmap_layout)
  
  
  umap_heatmap_plot
  
}

## Predict celltypes by set of markers
FindCellTypesByMarkers <- function(sobj) {
  
  ## If function can not find any significant genes for several cluster
  ## then heatmap should show all zeros for such clusters for each celltype
  ## and such cluster will not be associated with any cell type, preserving the initial index
  ## this option works because RenameIdents can utilize partitial labels like c(`1` = "AAA", `4` ="BBB")
  ## !! If it was detected only one celltype and even if the expression is negative for this celltye
  ## the procedure anyway assign this celltype for cluster because other is zero
  
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
  
  DefaultAssay(sobj) <- "RNA"
  
  tmp <- map_dfr(names(biomarkers), function(name) {
    mm <- FindAllMarkers(sobj,features = biomarkers[[name]] ,logfc.threshold = -10, min.pct = 0.0) 
    
    if (is_empty(mm)){
      mm <- tibble(pval = 1, 
                   avg_log2FC = 0, 
                   pct.1 = 0, 
                   pct.2 = 0, 
                   p_val_adj = 1,
                   cluster = as.factor(0), 
                   gene = biomarkers[[name]][1])
    }
    mm %>% 
      group_by(cluster) %>% 
      summarise(score = mean(avg_log2FC)) %>%
      mutate(celltype = name)
  }) 
  
  DefaultAssay(sobj) <- "SCT"
  
  ## FindAllMarkers skips low expressed genes and produces as output empty rows for some celltypes
  ## To make full table without empty rows we have to append celltypes with zero scores for such of celltypes
  
  tmp <- left_join(tmp %>% expand(cluster,celltype), tmp) %>% 
    replace_na(list(score = 0))
  
  labels <- tmp %>% 
    group_by(cluster) %>%
    summarise(prediction = ifelse(score[which.max(score)]> 0, celltype[which.max(score)], "Unknown" )) %>% 
    ungroup() %>% 
    mutate(prediction = paste0(prediction,"(",cluster,")"))
  
  new_labels <- labels$prediction
  names(new_labels) <- labels$cluster
  
  ## swap names and values
  old_labels <- setNames(names(new_labels), new_labels)
  
  list(new_labels = new_labels, old_labels = old_labels, heatmap_table = tmp)   
}

#' Extract number and percent of cells info for each cluster about 
#'
#' @param seurat_obj - Seurat object with defined clusters
#'
#' @return data.frame with statistics
#' @export
#'
#' @examples
list_sample_cluster_stats <- function(seurat_obj){
  
  
  get_table_by_id <- function(sid, dff) {
    tmp <- dff
    
    ## for all samples together
    if (sid != "all") {
      tmp <- dff %>%  
        filter(sample_id == sid)
    } else {
      tmp <- dff %>% 
        group_by(seurat_clusters) %>% 
        summarise(ncells = sum(ncells)) %>% 
        ungroup()
    }
    
    tmp <- tmp %>% 
      mutate(all = sum(ncells), 
             cluster = as.character(seurat_clusters),
             pct = round(100*ncells/all,2)) %>% 
      select(cluster, ncells, pct)
    
    bind_rows(tmp, tmp %>% 
                summarise(cluster = "total", 
                          ncells = sum(ncells), 
                          pct = sum(pct)))
  }
  
  all_cells <- seurat_obj@meta.data %>% 
    select(seurat_clusters, sample_id) %>% 
    group_by_all() %>% 
    summarise(ncells = n()) %>% 
    ungroup()
  
  samples_names <- c("all", seurat_obj@meta.data$sample_id %>% unique() %>% sort())
  
  res <- map(samples_names, function(name){
    get_table_by_id(sid = name, dff = all_cells)
  }) %>% purrr::set_names(samples_names) 
  
  res
}

## create heatmap by table (cluster, celltype, score)
CreateCellTypesHeatmap <- function(df){
  ggplot(df, aes(x = celltype, y = as.factor(cluster))) +
    geom_tile(aes(fill = score),color= "gray50",size = 0.1)+
    scale_fill_gradient2(low = "blue", mid="white", high = "tomato")+
    geom_text(aes(label=round(score,2)), size = 6)+
    scale_x_discrete(position = "top") +
    xlab("")+
    ylab("")+
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(color = "black",size = 15),
          axis.text.y = element_text(color = "black", size = 15))
}

#' find optimal number of cluster using binary search for resolution parameter
#'
#' @param expected_clusters - Int number of clusters we would like to get
#' @param seurat_object Seurat object with default assay SCT (normalization using SCTransform)
#'
#' @return (float) new resolution if method did not find exact number of clusters in 15 iteration return closest resolution
#' @export
#'
#' @examples
binary_search_resolution <- function(expected_clusters, seurat_object) {
  x0 <- 0; x1 <- 1; new_resolution <- 0
  number_of_clusters <- 0
  counter <- 0 
  max_counter <- 15 ## search until max counter otherwise it could be infinite loop for some samples
  
  while (number_of_clusters != expected_clusters && counter < max_counter)
  {
    new_resolution <- x0+(x1-x0)/2
    
    number_of_clusters <- FindClusters(seurat_object, resolution = new_resolution) %>% 
      .@meta.data %>% 
      select(paste0("SCT_snn_res.",new_resolution)) %>% 
      group_by(paste0("SCT_snn_res.",new_resolution)) %>% 
      distinct() %>%
      summarise(n=n()) %>%
      pull(n)
    
    # print(new_resolution)
    # print(number_of_clusters)
    # print(str_interp("resolution: ${new_resolution} -> number_of_clusters: ${number_of_clusters} x0:${x0} x1:${x1}"))
    
    if (number_of_clusters < expected_clusters) {
      x0 <- new_resolution
    } 
    else {
      x1 <- new_resolution
    }
    
    counter <- counter+1
  }
  
  cat(paste0("==>BINARY SEARCH COUNTER: ",counter,"<=="))
  
  return(new_resolution)
}

## recalculate UMAP (PCA,HARMONY,UMAP,FIND.NEIGHBOURS,BIN.SEARCH,FIND.CLUSTERS)
recalc_umap_2 <- function(seurat_obj, npcs = 8, min.dist = 0.001, number_of_clusters = 8, opt_resolution = "auto", k.param = 20, runsct = F) {
  
  if (runsct) {
    seurat_obj <- SCTransform(seurat_obj, conserve.memory = T)
  }
  
  print("Sct - Find neighbours")
  seurat_obj <- RunPCA(seurat_obj, npcs = npcs) %>%
    RunHarmony(
      assay.use = "SCT",
      reduction = "pca",
      dims.use = 1:npcs,
      group.by.vars = "sample_id"
    ) %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs, min.dist = min.dist) %>%
    FindNeighbors(reduction = "harmony", dims = 1:npcs)
    
  print("calculate resolution")
  if (opt_resolution == "auto") {
    opt_resolution <- binary_search_resolution(expected_clusters = number_of_clusters, 
                                               seurat_object = seurat_obj)
  }
  
  print("compute clusters")
  seurat_obj <- FindClusters(seurat_obj, resolution = opt_resolution) %>% 
    NormalizeData(., assay = "RNA")
  
  seurat_obj
}

