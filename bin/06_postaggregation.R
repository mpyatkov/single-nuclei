#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Extract the user-defined cellbarcodes for all samples from aggregated.csv ')
    
    p <- add_argument(p, '--filtered', help="path to filtered matrix intronic-with-mono/genebody",
                      default='filtered_feature_bc_matrix.h5')
    
    p <- add_argument(p, '--userdefined_cellbarcodes', help='user defined cellbarcodes',
                      default="aggr-user-defined-barcodes.csv")
    
    p <- add_argument(p, '--nclusters', 
                      help='number of cluster on UMAP',
                      default="12")
    
    p <- add_argument(p, '--mt_percent', 
                      help='mt percent (should be less than user-defined percent otherwise will no effect)',
                      default=10.0)
    
    p <- add_argument(p, '--min_genes', 
                      help='min number of genes detected in each cell',
                      default=400)
    
    p <- add_argument(p, '--min_counts', 
                      help='min number of molecules detected within each cell',
                      default=400)

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

    p  <- add_argument(p, '--labels_for_umap', 
                       help='activate labels for UMAP plots {Yy, Nn}',
                       default="Y")

    p  <- add_argument(p, '--number_of_cores', 
                       help='some part of the script can be runned in parallel',
                       default=2)
        
    p <- add_argument(p, '--dotplot_gene_list', 
                      help='multiple column csv with a header, each column is list of genes for dotplot, header - plot title',
                      default="dotplot_gene_list.csv")
    
    p <- add_argument(p, '--featureplot_gene_list', 
                      help='one column csv with a header, each column list of genes for dotplot',
                      default="featureplot_gene_list.csv")
    
    p <- add_argument(p, '--seurat_cluster_labels', 
                      help='two column csv with a header (seurat_clusters, cluster_labels)',
                      default="seurat_cluster_labels.csv")
    
    
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

## parse arguments
## "1,2,3,4,5" -c(1,2,3,4,5)
str_to_vec <- function(s){
  str_split(s,",", simplify = T) %>% 
    str_trim() %>% 
    discard(~.=="") %>% 
    as.numeric()
}

## read cellranger results from filtered dataset for intronic counting (without monoexonic genes)
read_h5_matrix <- function(matrix_path, matrix_name) {
    print(paste0("Loading: ", matrix_name," - " ,matrix_path))
    mx <- Read10X_h5(matrix_path)
    mx <- CreateSeuratObject(mx)
    mx@meta.data$samples <- matrix_name
    mx@meta.data$CB <- rownames(mx@meta.data)
    mx
}

## find optimal number of cluster using binary search for resolution parameter
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

recalculate_seurat_object_aggr <- function(
    seurat_object,
    meta_data,
    mt_percent,
    min_genes,
    min_counts,
    vars_to_regress = "percent.mt") 
{
    meta_data <- as.data.frame(meta_data)
    rownames(meta_data) <- meta_data$CB
    
    print("==>RUN seurat object<==")
    new_seurat <- seurat_object[,meta_data$CB]
    
    print("==>RUN add meta data<==")
    new_seurat <- AddMetaData(new_seurat, meta_data) 
    
    print("==>RUN mt_percent, min_count, min_genes filtering<==")
    new_seurat[["percent.mt"]] <- PercentageFeatureSet(new_seurat, pattern = "^mt-")
    new_seurat <- subset(new_seurat, subset = percent.mt < mt_percent)
    new_seurat <- subset(new_seurat, subset = nFeature_RNA > min_genes & nCount_RNA > min_counts )
    
    print("==>RUN SCTransform<==")
    new_seurat <- SCTransform(new_seurat, vars.to.regres = vars_to_regress, conserve.memory = T)
    gc()
    
    new_seurat
}

## recalculate UMAP (PCA,HARMONY,UMAP,FIND.NEIGHBOURS,BIN.SEARCH,FIND.CLUSTERS)
recalc_umap <- function(new_seurat, npcs, harmony_vars = "batch", min.dist, number_of_clusters, opt_resolution = "auto", k.param = 20) {
    
    ## number_of_clusters will be ignored if opt_resolution != "auto"
    
    mito.genes <- grep(pattern = "^(MT|Mt|mt)-", x = rownames(x = new_seurat@assays$SCT@data), value = TRUE)
    var.genes.no <- dplyr::setdiff(new_seurat@assays$SCT@var.features, mito.genes)
    
    ## we do not need to use ScaleData on SCTranform data
    ## https://github.com/satijalab/seurat/discussions/4259
    # print("==>RUN ScaleData (probably is not required)<==")
    #new_seurat <- ScaleData(new_seurat, vars.to.regress = vars_to_regress)
    
    print("==>RUN PCA<==")
    
    print(npcs)
    new_seurat <- RunPCA(new_seurat, npcs = npcs, features = var.genes.no)
    
    print("==>RUN HARMONY<==")
    new_seurat <- RunHarmony(new_seurat, group.by.vars = harmony_vars, assay="SCT")
    
    print("==>RUN UMAP<==")
    new_seurat <- RunUMAP(new_seurat, reduction = "harmony", dims = 1:npcs, min.dist = min.dist) #, min.dist = 0.01
    #new_seurat <- RunUMAP(new_seurat, reduction = "harmony", dims = 1:10)
    new_seurat <- FindNeighbors(new_seurat, reduction = "harmony", dims = 1:npcs, k.param = k.param)
    
    print("==>RUN OPT RESOLUTION<==")
    if (opt_resolution == "auto") {
        opt_resolution <- binary_search_resolution(expected_clusters = number_of_clusters, 
                                                   seurat_object = new_seurat)
    }
    str_interp("Resolution: ${opt_resolution}") %>% print()
    new_seurat <- FindClusters(new_seurat, resolution =opt_resolution) 
    
    new_seurat
}

### extract clusters statistics from seurat_object
extractClusterStats <- function(seurat_obj, labels = NULL) {
  ## labels is data.frame with 2 columns (seurat_clusters, cluster_labels)
  
  cluster_stats <- seurat_obj@meta.data %>% 
    select(sample_id, seurat_clusters) %>% 
    mutate(all = n()) %>% 
    
    group_by(seurat_clusters) %>% 
    mutate(n_cl = n()) %>% 
    ungroup() %>% 
    
    group_by(sample_id) %>% 
    mutate(n_sample = n()) %>% 
    ungroup() %>% 
    
    group_by(seurat_clusters, sample_id) %>% 
    summarise(ncells = n(), pt_per_sample=round(100*ncells/n_sample,2), pt_per_cluster = round(100*ncells/n_cl,2)) %>% 
    distinct() %>% 
      ungroup() 
  
  ## if labels exist
  if (!labels %>% is.null()) {
    cluster_stats <- left_join(cluster_stats, labels, by = "seurat_clusters") %>% 
        relocate(cluster_labels, .after = seurat_clusters) %>%
        rowwise() %>%
        mutate(cluster_labels = paste0(cluster_labels,"(",seurat_clusters,")"))
  }
  
  cluster_stats
  
}

# get_table_by_id("G183M1", table_for_plots)
## extract table by sample id from table with clusters statistics
get_table_by_id <- function(sid, dff) {
  tmp <- dff
  
  if (sid != "all") {
    tmp <- dff %>%  
      filter(sample_id == sid)
  } else {
    tmp <- dff %>% 
      group_by(seurat_clusters) %>% 
      summarise(ncells = sum(ncells)) %>% 
      ungroup()
    
  }
  
  tmp <- tmp %>% ungroup() %>% mutate(all = sum(ncells), 
                        cluster = as.character(seurat_clusters),
                        pct = round(100*ncells/all,2)) %>% 
    select(cluster, ncells, pct)
  
  bind_rows(tmp, tmp %>% summarise(cluster = "total", ncells = sum(ncells), pct = sum(pct)))
}

## supplementary function to create umap
plot_umap <- function(obj, title, labels = NULL){
  tmp <- DimPlot(obj, reduction = "umap") + 
    ggtitle(title)
  
  ## If labels is required add them to umap plot
  if (!is.null(labels_for_umap)) {
    
    ## to make labels the same color as legend colors use the following code
    # aggr_umap_all_clusters <- LabelClusters(..., color = unique(ggplot_build(aggr_umap_all_clusters)$data[[1]]$colour), ...)
    
    tmp <- LabelClusters(tmp, 
                         id = "ident", 
                         size = 7, 
                         repel = T,  
                         box.padding = 1, 
                         segment.size = 0.3,
                         max.overlaps = 10)
  }
  tmp
}

## supplementary function to create dotplot
dotplot <- function(sobj, genes, title = "") {
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

## supplementary function to create table on the ggplot
plot_table <- function(df, title){
  table <- gridExtra::tableGrob(df, rows = NULL, theme = ttheme_minimal(base_size = 18))
  h <- grobHeight(table)
  w <- grobWidth(table)
  title <- textGrob(title, y=unit(0.5,"npc") + 1.5*h, vjust=0, gp=gpar(fontsize=20))
  table <- gTree(children = gList(title, table))
  table    
}

## rearrange columns for cluster_stats table
prepareMetaData <- function(cluster_stats) {
  spec <- build_wider_spec(cluster_stats, 
                           names_from = sample_id, 
                           values_from = c("ncells","pt_per_sample", "pt_per_cluster")) %>% 
    arrange(sample_id, .name)
  
  pivot_wider_spec(cluster_stats, spec)
}

### creates combined plot with layout (umap/table/dotplot)
combinedUmapTableDotplot <- function(seurat_obj, TITLE, labels_for_umap, genes_list, cluster_stats, sample_id = "all") {
  
  #### table for plots
  table_for_plots <- cluster_stats %>% 
    select(seurat_clusters, sample_id, ncells)
  
  ## dotplot for all samples together
  DefaultAssay(seurat_obj) <- "RNA"
  gc()
  
  
  ## rename clusters (set new_labels)
  seurat_obj <- RenameIdents(seurat_obj, labels_for_umap$new_labels)

  ## Plot all samples together on one UMAP without labels
  aggr_umap_all_clusters <- plot_umap(seurat_obj, title = TITLE, labels = labels_for_umap)
  
  ## rename clusters back (set old_labels)
  seurat_obj <- RenameIdents(seurat_obj, labels_for_umap$old_labels)
  
    
  ## TODO: at the moment all genes will be represented on the dotplots
  aggr_umap_all_clusters_dotplot <- dotplot(seurat_obj, genes_list)
  
  
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

### If the set of parameter for each individual parameter contains only one value
### we are going to use this function
### save meta.data, projections and umap plots for each individual sample (batch)
saveUmapDotplotTableIndividual <- function(seurat_obj, 
                                           genes_list,
                                           number_of_clusters, 
                                           npcs, 
                                           min.dist, 
                                           resolution = "auto", 
                                           k.param = 20,
                                           df_cluster_names = NULL) {
  
  #seurat_obj <- userdefined_mx_raw
  system.time(
    # function(new_seurat, npcs, var.genes.no, harmony_vars = "batch", min.dist, number_of_clusters) {
    seurat_obj <- recalc_umap(seurat_obj, 
                              npcs = npcs, 
                              min.dist = min.dist, 
                              number_of_clusters = number_of_clusters, 
                              opt_resolution = resolution,
                              k.param = k.param)
    
  )
  
  gc()
    
  #### normalization for RNA assay ####
  seurat_obj <- seurat_obj %>% 
    NormalizeData(., assay = "RNA") 
  
  ## calculate cluster stats information from seurat object
  cluster_stats <- extractClusterStats(seurat_obj, labels = df_cluster_names)
  
  #### attaching cluster labels if exist ####
  ## user_defined_labels = NULL OR list(old_labels,new_labels)
  user_defined_labels <- if (!is.null(df_cluster_names)) {  
      ## preparing cluster labels
      newlabels <- unique(cluster_stats$cluster_labels)
      names(newlabels) <- levels(Idents(seurat_obj))

      ## remove NA names
      newlabels <- newlabels[!is.na(names(newlabels))]
      
      ## make a copy of old labels, will be recovered latter
      oldlabels <- setNames(names(newlabels), newlabels)
      
      ## assign new cluster labels
      list(new_labels = newlabels, old_labels = old_labels)
  }
    
  ## export cluster stats (seurat_obj, cluster labels) for all clusters
  prepareMetaData(cluster_stats) %>% write_csv("aggr_cluster_stats.csv")
  
  
  ## export meta-data (clusters, cellbarcodes, batch and other columns)
  tmp_meta <- seurat_obj@meta.data %>% 
    select(CB,sample_id:percent.mt, -percent.mt, seurat_clusters) 
  
  ## add information about user_defined clusters
  if(!is.null(user_defined_labels)){
    tmp_meta <- tmp_meta %>% 
      left_join(., user_defined_labels$new_labels %>% tibble::enframe(name="seurat_clusters", value = "cluster_labels", by="seurat_clusters"))
  }
  
  tmp_meta %>% write_csv("aggr_meta_data.csv")

  
  #### projections export ####
  seurat_obj@reductions$umap@cell.embeddings %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(var = "Barcode") %>% 
    write_csv("aggr_projections.csv")
  
  
  
  ## loop over all samples
  ## create umap/dotplot/table plot for each sample separately
  
  ### Assign predicted cell types labels
  ## get new_labels,old_labels,heatmap_table
  celltypes_list <- FindCellTypesByMarkers(sobj = seurat_obj)

  ## if the user has provided custom labels, use them
  if (!is.null(user_defined_labels)) {  
    celltypes_list$new_labels <- user_defined_labels$new_labels
    celltypes_list$old_label <- user_defined_labels$old_labels 
  }
  
  celltypes_heatmap <- CreateCellTypesHeatmap(celltypes_list$heatmap_table)
  umap_heatmap_layout <- "
AAAAA####
AAAAA#BBB
AAAAA#BBB
AAAAA####
"
  #######################################
  #### umap/dotplot/table plot for all CB
  #######################################
  
  #### umap plot with all clusters ####
  #### make title
  resolution_text <- grep("SCT_snn_res.", colnames(seurat_obj@meta.data), value = T) %>% 
    gsub("SCT_snn_res.","", .) %>% tail(., n=1)
  
  TITLE <- str_interp("Aggregated. Resolution: ${resolution_text}, min.dist: ${min.dist}, #PCs: ${npcs}, #clusters: ${number_of_clusters} k.param: ${k.param}" )
  
  all_samples_combined <- combinedUmapTableDotplot(seurat_obj = seurat_obj, 
                                                   TITLE = TITLE, 
                                                   labels_for_umap = celltypes_list,
                                                   cluster_stats = cluster_stats,
                                                   sample_id = "all",
                                                   genes_list = genes_list)
  all_samples_combined <- cowplot::plot_grid(all_samples_combined+celltypes_heatmap+plot_layout(design = umap_heatmap_layout))
  
  ###########################################
  #### Plot each batch separately with dotplot
  ###########################################
  ## https://github.com/satijalab/seurat/issues/1825
  ob.list <- SplitObject(seurat_obj, split.by = "batch")
  all_samples_by_sample_id <- lapply(X = names(ob.list), FUN = function(x) {
    celltypes_list <- FindCellTypesByMarkers(sobj = ob.list[[x]])
    celltypes_heatmap <- CreateCellTypesHeatmap(celltypes_list$heatmap_table)
    
    tmp <- combinedUmapTableDotplot(seurat_obj = ob.list[[x]], 
                             TITLE = TITLE, 
                             labels_for_umap = celltypes_list,
                             cluster_stats = cluster_stats,
                             sample_id = x,
                             genes_list = genes_list)
    umap_heatmap_plot <- cowplot::plot_grid(tmp+celltypes_heatmap+plot_layout(design = umap_heatmap_layout))
    # umap <- plot_umap(ob.list[[x]],"empty title", labels = labels_for_umap)
    # dtplt <- dotplot(ob.list[[x]], unlist(flatten(genes_list)))
    # table <- plot_table(df = get_table_by_id(x, table_for_plots), title = x)
    # tmp <- umap+dtplt+table+plot_layout(design = layout)+plot_annotation(title = x)
    # cowplot::plot_grid(tmp)
  })
  
  #### TODO (from postaggregation.R)
  ## ? create umap for batch 
  ## create dotplot separately for each group of genes (list of genes, seurat_obj, file)
  ## create featureplot for each group of genes (list of genes, seurat_obj, file)
  
  #list(aggr_umap_all_clusters, aggr_group_by_batch, aggr_split_by_batch) %>%
  append(list(all_samples_combined), all_samples_by_sample_id) %>% 
    marrangeGrob(nrow =1, ncol=1) %>% 
    ggsave(filename = "aggr_umaps.pdf", width = 20.0, height = 15.50)
  gc()
  
  return(seurat_obj)
}

### simplified version of previous function. 
### Using it only in case we have multiple options for some parameters
### example 'pc = "6,10,20"'
### generates multiple pdfs with umap/table/dotplot for each combination of parameters
saveUmapDotplotTable <- function(seurat_obj, 
                                 genes_list,
                                 number_of_clusters, 
                                 npcs, 
                                 min.dist, 
                                 resolution = "auto", 
                                 k.param = 20,
                                 counter = "000") {
  
  print(npcs)
  #seurat_obj <- userdefined_mx_raw
  system.time(
    # function(new_seurat, npcs, var.genes.no, harmony_vars = "batch", min.dist, number_of_clusters) {
    seurat_obj <- recalc_umap(seurat_obj, 
                              npcs = npcs, 
                              min.dist = min.dist, 
                              number_of_clusters = number_of_clusters, 
                              opt_resolution = resolution,
                              k.param = k.param)
    
  )
  
  gc()
  
  #DefaultAssay(userdefined_mx) <- "SCT"
  #userdefined_mx <- recalculate_clusters(13, userdefined_mx)
  print("Extract clusters")
  cluster_stats <- extractClusterStats(seurat_obj)

  ## TODO: export by enabled flag
  
  ###########################################
  #### umap plot with all clusters ####
  ###########################################
  
  #### scaling and normalization for RNA assay ####
  seurat_obj <- seurat_obj %>% 
    NormalizeData(., assay = "RNA") 
  
  
  #### make title
  resolution_text <- grep("SCT_snn_res.", colnames(seurat_obj@meta.data), value = T) %>% 
    gsub("SCT_snn_res.","", .) %>% tail(., n=1)
  
  TITLE <- str_interp("${counter} Aggregated. Resolution: ${resolution_text}, min.dist: ${min.dist}, #PCs: ${npcs}, #clusters: ${number_of_clusters} k.param: ${k.param}" )
  
  ### Assign predicted cell types labels
  ## get new_labels,old_labels,heatmap_table
  celltypes_list <- FindCellTypesByMarkers(sobj = seurat_obj)
  celltypes_heatmap <- CreateCellTypesHeatmap(celltypes_list$heatmap_table)
  
  umap_heatmap_layout <- "
AAAAA####
AAAAA#BBB
AAAAA#BBB
AAAAA####
"
  
  all_samples_combined <- combinedUmapTableDotplot(seurat_obj = seurat_obj, 
                                                   TITLE = TITLE, 
                                                   labels_for_umap = celltypes_list,
                                                   cluster_stats = cluster_stats,
                                                   sample_id = "all",
                                                   genes_list = genes_list)
  
  umap_heatmap_plot <- all_samples_combined+celltypes_heatmap+plot_layout(design = umap_heatmap_layout)

  name <- str_interp("${counter}_aggr_res${resolution_text}_pc${npcs}_dst${min.dist}_cl${number_of_clusters}_kp${k.param}.pdf")
  
  umap_heatmap_plot %>% ggsave(filename = name, width = 20, height = 15.50)

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
  
  ## FindAllMarkers skips low expressed genes and produces as output empty rows for some celltypes
  ## To make full table without empty rows we have to append celltypes with zero scores for such of celltypes
  ## TODO: make it to tmp also, replace celltype[which.max(score)] to function which will return
  ## celltype only for positive log2FC, and "unknown" for the rest 
    
  heatmap_table <- left_join(tmp %>% expand(cluster,celltype), tmp) %>% 
    replace_na(list(score = 0))
  
  labels <- tmp %>% 
    group_by(cluster) %>% 
    summarise(prediction = celltype[which.max(score)]) %>% 
    ungroup() %>% 
    mutate(prediction = paste0(prediction,"(",cluster,")"))
  
  new_labels <- labels$prediction
  names(new_labels) <- labels$cluster
  
  ## swap names and values
  old_labels <- setNames(names(new_labels), new_labels)
  
  list(new_labels = new_labels, old_labels = old_labels, heatmap_table = heatmap_table)   
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


################
#### PARAMETERS ####
################
arg_filtered_mx <- argv$filtered
arg_userdefined <- argv$userdefined_cellbarcodes

## files with additional info
dotplot_gene_list_arg <- argv$dotplot_gene_list
feature_gene_list_arg <- argv$featureplot_gene_list
seurat_cluster_labels <- argv$seurat_cluster_labels

## base parameters
mt_percent <- argv$mt_percent
min_genes <- argv$min_genes
min_counts <- argv$min_counts

labels_for_umap <- ifelse(argv$labels_for_umap %in% c('Y', 'y'), TRUE, FALSE)

## maybe vector parameters
umap_resolution <- ifelse(argv$umap_resolution == "auto", "auto", str_to_vec(argv$umap_resolution))
number_of_clusters <- str_to_vec(argv$nclusters)
umap_kparam <- str_to_vec(argv$umap_kparam)
umap_min_dist <- str_to_vec(argv$umap_min_dist)
umap_npcs <- str_to_vec(argv$umap_npcs)

## if multiple_parameters is TRUE do not export meta.data and projections
multiple_parameters <- list(umap_resolution, number_of_clusters, umap_kparam, umap_min_dist, umap_npcs) %>%
    map(length) %>%
    flatten_int() %>%
    discard(~.==1) %>%
    length() %>%
    map(`!=`,0) %>% unlist()

paste("Multiple params? ", multiple_parameters) %>% print
## print(list(umap_resolution, number_of_clusters, umap_kparam, umap_min_dist, umap_npcs))
## print(list(umap_resolution, number_of_clusters, umap_kparam, umap_min_dist, umap_npcs) %>%  map(length))
## print(list(umap_resolution, number_of_clusters, umap_kparam, umap_min_dist, umap_npcs) %>%   map(class))
## print(umap_resolution)
## print(number_of_clusters)
## print(umap_kparam)
## print(umap_min_dist)
## print(umap_npcs)

# arg_filtered_mx <- "./filtered_feature_bc_matrix.h5"
# arg_userdefined <- "./aggr-user-defined-barcodes.csv"
# number_of_clusters <- 9
# mt_percent <- 10.0
# min_genes <- 400
# min_counts <- 400
# dotplot_gene_list_arg <- "./dotplot_gene_list.csv"
# feature_gene_list_arg <- "./featureplot_gene_list.csv"
# seurat_cluster_labels <- "seurat_cluster_labels.csv"
# labels_for_umap <- "Y"

########################
#### MAIN ####
########################

###########################################
#### TODO: read gene list
###########################################

## only for debug (default_genes_list)
# genes_list <- c('Cyp2e1', 'Glul', 'Oat', 'Gulo', 'Ass1', 'Hamp', 'Gstp1', 'Ubb', 'Selenbp2', 'Sds', 'Cyp2f2', 'Pck1', 
#                 'Hal', 'Il2rb', 'Cxcr6', 'Gzma', 'Csf3r', 'S100a6', 'Pkhd1', 'Sox9', 'Epcam', 'Krt7', 'Krt19', 'Irf8', 
#                 'Itgax', 'Clec4f', 'Csf1r', 'Jchain', 'Cd79a', 'Cd79b', 'Top2a', 'Stab2', 'Kdr', 'Aqp1', 'Fcgr2b', 
#                 'Gpr182', 'Ebf1', 'Skap1', 'Ptprc', 'Ank3', 'Dcn', 'Colec11', 'Ecm1', 'Alb', 'Ttr', 'Apoa1', 'Serpina1c', 
#                 'mt-Atp6', 'mt-Atp8', 'mt-Co3', 'mt-Co1', 'mt-Co2', 'mt-Nd2', 'mt-Nd4', 'mt-Cytb', 'mt-Nd1', 'mt-Nd3', 
#                 'mt-Nd4l', 'mt-Nd5', 'mt-Nd6')

## read gene list and remove NA because in most case the gene list
## will not be the same by size (default_genes)
full_genes_list <- read_csv(dotplot_gene_list_arg, col_names = T, show_col_types = FALSE) %>% 
  as.list() %>% 
  map(discard, is.na)

## default genes_list for dotplot with umaps
genes_list <- full_genes_list$default_genes

## genes list for individual dotplots
genes_by_groups_wo_default <- names(full_genes_list) %>% 
  discard(.=="default_genes") %>% 
  full_genes_list[.]

## read genes for featureplots
feature_genes_list <- read_csv(feature_gene_list_arg, col_names = T) %>% pull(genes)

userdefined_meta <- read_csv(arg_userdefined, col_names = T, show_col_types = FALSE) %>% 
    mutate(batch = sample_id)

userdefined_mx <- read_h5_matrix(arg_filtered_mx, "userdefined")

system.time(
    userdefined_mx_raw <- recalculate_seurat_object_aggr(seurat_object = userdefined_mx, 
                                                     meta_data = userdefined_meta,
                                                     mt_percent = mt_percent,
                                                     min_genes = min_genes,
                                                     min_counts = min_counts,
                                                     vars_to_regress = c("percent.mt","nCount_RNA"))
)

gc()

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

# tictoc::tic("Whole time:")

if (multiple_parameters) {
  print("Make multiple plot")
  library(foreach)
  library(doParallel)
  registerDoParallel(as.numeric(argv$number_of_cores))
  foreach(i = params) %dopar% {
  
  ## foreach(i = params) %do% {

      ## print(length(genes_list))
      ## i$number_of_clusters %>% print
      ## i$npcs %>% print
      ## i$resolution %>% print
      ## i$k.param %>% print
      ## i$min.dist %>% print
      
    saveUmapDotplotTable(seurat_obj = userdefined_mx_raw,
                         genes_list = genes_list,
                         number_of_clusters = i$number_of_clusters,
                         npcs = i$npcs,
                         min.dist = i$min.dist,
                         resolution = i$resolution,
                         k.param = i$k.param,
                         counter = i$ix)
    
    gc()
  }
}else {
  print("Make individual plot")

  # ## read file with cluster names (2 columns seurat_clusters, seurat_names)
  df_cluster_names <- read_csv(seurat_cluster_labels, col_names = T, show_col_types = FALSE) %>% 
      mutate(seurat_clusters = as.factor(seurat_clusters))


  ## if seurat_clusters == cluster_labels, do nothing
  if (all(df_cluster_names[,1] == df_cluster_names[,2])) {
      df_cluster_names <- NULL
  }

  local_params <- flatten(params)
  
  ## I rewrite original userdefined_mx_raw because the raw object 
  ## is not required for downstream analysis
  userdefined_mx_raw <- saveUmapDotplotTableIndividual(seurat_obj = userdefined_mx_raw, 
                                 genes_list  = genes_list,
                                 number_of_clusters = local_params$number_of_clusters, 
                                 npcs =local_params$npcs, 
                                 min.dist = local_params$min.dist, 
                                 resolution = local_params$resolution, 
                                 k.param = local_params$k.param,
                                 df_cluster_names = df_cluster_names) 
  gc()
  
  ####
  ## individual dotplots
  ####
  names(genes_by_groups_wo_default) %>% 
    map(function(name) {
      genes <- genes_by_groups_wo_default[[name]]
      ngenes <- length(genes)
      dotplot(userdefined_mx_raw, genes, name)
    }) %>% 
    marrangeGrob(nrow = 4, ncol = 1) %>%
    ggsave(filename = "aggr_dotplots_by_groups.pdf", width = 13.9, height = 18)
  
  ####
  ## individual featureplots
  ####
  DefaultAssay(userdefined_mx_raw) <- "RNA"
  feature_genes_list %>% 
    map(function(gname) {
      FeaturePlot(userdefined_mx_raw, features = gname, pt.size = 1, order = T, raster = T)
    }) %>% 
    marrangeGrob(nrow = 2, ncol = 2, layout_matrix = matrix(1:4, 2, 2, TRUE)) %>%
    #grid.arrange(grobs = ., ncol = 2, as.table = FALSE) %>% 
    ggsave(filename = "aggr_featureplots.pdf", width = 11, height = 8.5)
  DefaultAssay(userdefined_mx_raw) <- "SCT"

  ## save rds
  saveRDS(userdefined_mx_raw, "module2_aggregated_seurat.rds")
  
}

  # tictoc::toc()

fn <- "Rplots.pdf"
if (file.exists(fn)) {file.remove(fn)}


