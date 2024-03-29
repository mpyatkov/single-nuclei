## "Single-nucleus RNA sequencing of pre-malignant liver reveals disease-associated hepatocyte state with HCC prognostic potential"
## Carlessi et al., 2023, Cell Genomics 3

CELL_PAPER_MARKERS <- list(
  `Hep`= c("Egfr","Pck1","Slc7a2","Cps1","Cyp7b1"),
  `Mesenchy/HSC` = c("Dcn","Lsamp","Ank3","Col14a1", "Bmp5"),
  `Endo` = c("Stab2","F8","Fgd5","Ptprb","Bmp6"),
  `BECs/Cholang` = c("Pkhd1","Ctnnd2","Erbb4","Hnf1b"),
  `Myeloid` = c("Cd5l","Adgre1","Myo1f","Clec4f","Sirpa"),
  `B cells` = c("Ebf1","Pax5","Aff3","Prkcb"),
  `T/NK Cells` = c("Skap1","Bcl11b","Lef1","Cd247"),
  `pDCs` = c("Siglech","Runx2","Mctp2","Card11"),
  `Mesothelial` = c("Bnc2","Sulf1","Grip1","Mast4"),
  `daHep` = c("Abcc4","Cachd1","Ehhadh","lnc12608","Nt5e","Gsta1","Samd4","Gclc","Plin2","Nrg1",
              "Spag5", "Taco1", "Klf6", "Mrnip","Gm28153","Mapt","Glis3","Robo1",
              "Peak1","Igsf11","Dst","Dhx40","Slc7a11","Krt8") # lnc12608 = Pvt1
)

DEFAULT_MARKER_LIST = list(
  `Standard biomarkers` = c('Cyp2e1', 'Glul', 'Oat', 'Gulo', 'Ass1', 'Hamp', 'Gstp1', 'Ubb', 'Selenbp2', 'Sds', 'Cyp2f2', 'Pck1',
                         'Hal', 'Il2rb', 'Cxcr6', 'Gzma', 'Csf3r', 'S100a6', 'Pkhd1', 'Sox9', 'Epcam', 'Krt7', 'Krt19', 'Irf8',
                         'Itgax', 'Clec4f', 'Csf1r', 'Jchain', 'Cd79a', 'Cd79b', 'Top2a', 'Stab2', 'Kdr', 'Aqp1', 'Fcgr2b',
                         'Gpr182', 'Ebf1', 'Skap1', 'Ptprc', 'Ank3', 'Dcn', 'Colec11', 'Ecm1', 'Alb', 'Ttr', 'Apoa1', 'Serpina1c',
                         'mt-Atp6', 'mt-Atp8', 'mt-Co3', 'mt-Co1', 'mt-Co2', 'mt-Nd2', 'mt-Nd4', 'mt-Cytb', 'mt-Nd1', 'mt-Nd3',
                         'mt-Nd4l', 'mt-Nd5', 'mt-Nd6'),
  
  `PC(extreme)` = c("Slc1a2", "Slc16a10", "Tox", "Lhpp", "Tenm3", "Rcan2", "Slc22a1", "lnc13654"),
  `HSC` = c("Ntm", "Nrxn1", "Ccbe1", "Gpc6", "Egfem1", "Hgf", "Zfp804b", "Nkain2"),
  `Endo` = c("Meis2", "Fbxl7", "Plekhg1"), 
  `Cholang` = c("Glis3", "Pdgfd"), 
  `Dividing` = c("Cenpp","Ect2", "Tpx2", "Hmmr","Anln", "Diaph3"),
  `Immune` = c("Inpp4b", "Bank1", "Fyb", "Hdac9", "Slc8a1") 
)

DEFAULT_HEATMAP_BIOMARKERS <- list(
  PC = c("Glul","Gulo","Oat","Cyp2e1"),
  PP = c("Pck1","Cyp2f2","Hal"),
  Endo=c("Stab2","Meis2","Fbxl7","Ptprb","F8","Plekhg1", "Fgd5", "Bmp6"),
  HSC=c("Dcn","Lsamp","Ank3","Col14a1", "Bmp5", "Ntm", "Nrxn1", "Ccbe1", "Gpc6", "Egfem1", "Hgf", "Zfp804b", "Nkain2"),
  Kupffer = c("Clec4f","Csf1r","Bank1","Fyb", "Hdac9", "Slc8a1"),
  Immune = c("Ebf1","Skap1","Inpp4b"), ##was only Ptprc
  Cholang=c("Pkhd1","Glis3","Pdgfd","Erbb4", "Ctnnd2", "Hnf1b"), ## added Pkhd1  
  Div=c("Top2a","Cenpp","Ect2","Tpx2","Hmmr","Anln","Diaph3"),
  `PC(Ext)`=c("Slc1a2","Slc16a10","Tox","Lhpp","Tenm3","Rcan2","Slc22a1","lnc13654"),
  `B cells` = c("Ebf1","Pax5","Aff3","Prkcb"),
  `T/NK Cells` = c("Skap1","Bcl11b","Lef1","Cd247"),
  `pDCs` = c("Siglech","Runx2","Mctp2","Card11"),
  `Mesothelial` = c("Bnc2","Sulf1","Grip1","Mast4")
)



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

combinedUmapTable_v2 <- function(seurat_obj, TITLE, labels_for_umap = NULL, cluster_stats, sample_id = "all") {
  
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
  

  clusters_stats_table_plot <- plot_table(get_table_by_id(sample_id, table_for_plots), title = "")
  
  ## set layout for 11 plots
  layout <- "
    AAAAB
    AAAAB
    AAAAB
    "
  combined_all_clusters <- cowplot::plot_grid((aggr_umap_all_clusters+theme(plot.margin=margin(0,1,0,0,"cm")))+
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
mp_plot_umap <- function(seurat_obj, title = "", draw_labels = T, ...){
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
  DotPlot(sobj, features = genes, assay = "RNA", ...)+
    RotatedAxis()+ 
    {if(title !="") ggtitle(title)} +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal", 
          legend.justification = "center",
          legend.text=element_text(size=8),
          legend.title = element_text(size = 10))
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
                                 TITLE = "",
                                 heatmap_markers = NULL,
                                 rename_umap_labels = T) {
  
  #### PREDICT CELL TYPES
  ### Assign predicted cell types labels
  ## get new_labels,old_labels,heatmap_table
  
  if(is.null(heatmap_markers)){
    heatmap_markers <-  DEFAULT_HEATMAP_BIOMARKERS
  }
  
  celltypes_list <- FindCellTypesByMarkers(sobj = seurat_obj, biomarkers = heatmap_markers)
  celltypes_heatmap <- CreateCellTypesHeatmap(celltypes_list$heatmap_table)
  
  # ## dotplot for all samples together
  # DefaultAssay(seurat_obj) <- "RNA"

  if (rename_umap_labels) {
      #### ADD LABELS FOR UMAP IF EXISTS
      ## rename clusters (set new_labels)
      seurat_obj <- RenameIdents(seurat_obj, celltypes_list$new_labels)

      ## Plot all samples together on one UMAP without labels
      aggr_umap_all_clusters <- mp_plot_umap(seurat_obj, title = TITLE, draw_labels = T)
      
      ## rename clusters back (set old_labels)
      seurat_obj <- RenameIdents(seurat_obj, celltypes_list$old_labels)
  } else {
      ## just plot UMAP
      aggr_umap_all_clusters <- mp_plot_umap(seurat_obj, title = TITLE, draw_labels = T)
  }
    
  #### CREATE DOTPLOT
  ## TODO: at the moment all genes will be represented on the dotplots
  aggr_umap_all_clusters_dotplot <- mp_dotplot(seurat_obj, genes_list)
  
  #### CREATE STATS TABLE WITH CLUSTERS
  clusters_stats_table_plot <- mp_plot_table(list_sample_cluster_stats(seurat_obj)[[sample_id]], title = "")
  
  layout <- "
    AAAA#C
    AAAA#C
    AAAA#C
    BBBBBB
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

UmapDotplotTableHeatmapPlot_v2 <- function(seurat_obj, 
                                 genes_list = NULL,
                                 additional_genes_list = NULL,
                                 sample_id = "all", 
                                 TITLE = "",
                                 heatmap_markers = NULL,
                                 rename_umap_labels = T) {
  
  #### PREDICT CELL TYPES
  ### Assign predicted cell types labels
  ## get new_labels,old_labels,heatmap_table

  if(is.null(heatmap_markers)){
    heatmap_markers <-  DEFAULT_HEATMAP_BIOMARKERS
  }
  
  if (is.null(genes_list)) {
    genes_list <- DEFAULT_MARKER_LIST
  }
  
  celltypes_list <- FindCellTypesByMarkers(sobj = seurat_obj, biomarkers = heatmap_markers)
  celltypes_heatmap <- CreateCellTypesHeatmap(celltypes_list$heatmap_table, celltype_order = names(heatmap_markers), rotate_x = T)
  
  ## CREATE UMAP
  DefaultAssay(seurat_obj) <- "RNA"

  if (rename_umap_labels) {
      #### ADD LABELS FOR UMAP IF EXISTS
      ## rename clusters (set new_labels)
      seurat_obj <- RenameIdents(seurat_obj, celltypes_list$new_labels)

      ## Plot all samples together on one UMAP without labels
      aggr_umap_all_clusters <- mp_plot_umap(seurat_obj, title = TITLE, draw_labels = T)
      
      ## rename clusters back (set old_labels)
      seurat_obj <- RenameIdents(seurat_obj, celltypes_list$old_labels)
  } else {
      ## just plot UMAP
      aggr_umap_all_clusters <- mp_plot_umap(seurat_obj, title = TITLE, draw_labels = T)
  }
    
  #### CREATE DOTPLOT
  ## TODO: at the moment all genes will be represented on the dotplots
  default_dotplot <- mp_dotplot(seurat_obj, genes_list)
  
  ##
  additional_dotplot <- NULL
  ## Additional dotplot
  if (!is.null(additional_genes_list)){
    additional_dotplot <- mp_dotplot(seurat_obj, additional_genes_list)
  }
  
  #### CREATE STATS TABLE WITH CLUSTERS
  clusters_stats_table_plot <- mp_plot_table(list_sample_cluster_stats(seurat_obj)[[sample_id]], title = "")
  
  ## LAYOUTS
  umap_stats_layout <- "
    AAAA#B
    AAAA#B
    AAAA#B
    "
  
  ## UMAP-STATS
  umap_stats_plot <- cowplot::plot_grid((aggr_umap_all_clusters+theme(plot.margin=margin(0,1,0,0,"cm")))+
                                                clusters_stats_table_plot+
                                                plot_layout(design = umap_stats_layout)+plot_annotation(title = sample_id))
  
  #### ATTACH HEATMAP TO UmapStatsDotplot plot
  umap_stats_heatmap_layout1 <- "
AAAAA#####
AAAAA#BBBB
AAAAA#BBBB
AAAAA#####
CCCCCCCCCC
"

  ### for additional dotplot
  umap_stats_heatmap_layout2 <- "
AAAAA#####
AAAAA#BBBB
AAAAA#BBBB
AAAAA#BBBB
AAAAA#####
CCCCCCCCCC
DDDDDDDDDD
"

  umap_stats_heatmap_dotplot_plot <- if(is.null(additional_genes_list)){
    umap_stats_plot+
      celltypes_heatmap+
      default_dotplot+
      plot_layout(design = umap_stats_heatmap_layout1)

  } else {
    umap_stats_plot+
      celltypes_heatmap+
      default_dotplot+
      additional_dotplot+
      plot_layout(design = umap_stats_heatmap_layout2)
  }

  umap_stats_heatmap_dotplot_plot
}

## Predict celltypes by set of markers
FindCellTypesByMarkers <- function(sobj, biomarkers = NULL) {
  
  if(is.null(biomarkers)) {
    stop("parameter biomarkers should be specied (now NULL)")
  }
  
  ## If function can not find any significant genes for several cluster
  ## then heatmap should show all zeros for such clusters for each celltype
  ## and such cluster will not be associated with any cell type, preserving the initial index
  ## this option works because RenameIdents can utilize partitial labels like c(`1` = "AAA", `4` ="BBB")
  ## !! If it was detected only one celltype and even if the expression is negative for this celltye
  ## the procedure anyway assign this celltype for cluster because other is zero
  
  # biomarkers <- heatmap_markers
  # 
  # if (is.null(heatmap_markers)) {
  #   biomarkers <- list(
  #     PC = c("Glul","Gulo","Oat","Cyp2e1"),
  #     PP = c("Pck1","Cyp2f2","Hal"),
  #     Kupffer = c("Clec4f","Csf1r","Bank1","Fyb", "Hdac9", "Slc8a1"),
  #     Immune = c("Ebf1","Skap1","Inpp4b"), ##was only Ptprc
  #     HSC=c("Colec11","Dcn","Ecm1","Ntm","Nrxn1","Ccbe1","Zfpm2","Gpc6","Egfem1","Lsamp","Hgf","Zfp804b","Bmp5","Nkain2"), 
  #     Endo=c("Stab2","Gpr182","Kdr","Fcgr2b","Aqp1","Meis2","Fbxl7","Ptprb","F8","Plekhg1"),
  #     Div=c("Top2a","Cenpp","Ect2","Tpx2","Hmmr","Anln","Diaph3"),
  #     Cholang=c("Epcam","Krt19","Krt7","Sox9","Pkhd1","Bicc1","Glis3","Pdgfd","Erbb4","Thsd4"), ## added Pkhd1
  #     Ext_PC=c("Slc1a2","Slc16a10","Tox","Lhpp","Tenm3","Rcan2","Slc22a1","lnc13654")
  #   )
  # }
  # 

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
  
  tmp <- left_join(tmp %>% tidyr::expand(cluster,celltype), tmp) %>% 
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
        summarise(ncells = sum(ncells), 
                  avg.counts = round(mean(avg.counts), 0),
                  avg.genes = round(mean(avg.genes),0 )) %>% 
        ungroup()
    }
    
    tmp <- tmp %>% 
      mutate(all = sum(ncells), 
             cluster = as.character(seurat_clusters),
             pct = round(100*ncells/all,2)) %>% 
      select(cluster, ncells, pct, avg.counts, avg.genes)
    
    bind_rows(tmp, tmp %>% 
                summarise(cluster = "total", 
                          ncells = sum(ncells), 
                          pct = sum(pct),
                          avg.counts = round(mean(avg.counts),0), 
                          avg.genes = round(mean(avg.genes),0)))
  }
  
  all_cells <- seurat_obj@meta.data %>% 
    select(seurat_clusters, sample_id, counts = nCount_RNA, genes = nFeature_RNA) %>% 
    group_by(seurat_clusters, sample_id) %>% 
    summarise(ncells = n(),
              avg.counts = round(mean(counts),0),
              avg.genes = round(mean(genes),0)) %>% 
    ungroup()
  
  samples_names <- c("all", seurat_obj@meta.data$sample_id %>% unique() %>% sort())
  
  res <- map(samples_names, function(name){
    get_table_by_id(sid = name, dff = all_cells)
  }) %>% purrr::set_names(samples_names) 
  
  res
}

## create heatmap by table (cluster, celltype, score)
CreateCellTypesHeatmap <- function(df, labels_text_size = 6, xaxis_text_size = 15, yaxis_text_size = 15, rotate_x = FALSE, celltype_order = NULL){
  
  if(!is.null(celltype_order)){
    df$celltype = factor(df$celltype, levels = celltype_order, ordered = TRUE)
  }
  
  tmp <- ggplot(df, aes(x = celltype, y = as.factor(cluster))) +
    geom_tile(aes(fill = score),color= "gray50",size = 0.1)+
    scale_fill_gradient2(low = "blue", mid="white", high = "tomato")+
    geom_text(aes(label=round(score,2)), size = labels_text_size)+
    scale_x_discrete(position = "top") +
    xlab("Markers")+
    ylab("Clusters")+
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(color = "black", size = xaxis_text_size),
          axis.text.y = element_text(color = "black", size = yaxis_text_size))
  
  if (rotate_x){
    # tmp <- tmp+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = "black", size = xaxis_text_size))
    tmp <- tmp+theme(axis.text.x = element_text(angle = 90, hjust=0, color = "black", size = xaxis_text_size),
                     axis.text.x.top = element_text(vjust = 0.5))
  }
  tmp
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

mp_plot_samples <- function(hep_seurat, sorted_genes, samples.list1, title) {
  #cols <- c("white", "#F7EEF7","#E5CCFF", "red3","red4") ## red palette
  cols <- rev(c("#225ea8","#2171b5","#6baed6","#eff3ff","#fbb4b9","#fbb4b9"))# blue,red colors
  
  p1 <- wrap_elements(mp_dotplot(hep_seurat, sorted_genes, title = "ALL") + scale_colour_gradientn(colors = cols))
  
  p2 <- map(names(samples.list1), function (name) {
    wrap_elements(mp_dotplot(samples.list1[[name]], sorted_genes, title = name) + scale_colour_gradientn(colors = cols))
  }) %>% reduce(`/`)
  
  p3 <- wrap_elements(p1/p2+plot_layout(heights = c(1,4)))&
    plot_annotation(title = title, 
                    theme = theme(plot.title = element_text(size = 18)))
  p3
}
