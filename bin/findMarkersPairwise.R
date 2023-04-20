#' Add id_cluster column to sobject meta.data
#'
#' @param seurat_obj (Seurat object)
#'
#' @return Seurat object
#' @export
#'
#' @examples
add_id_cluster <- function(seurat_obj) {
  new_identities <- seurat_obj@meta.data %>% select(CB, sample_id, seurat_clusters) %>% 
    mutate(id_cluster = str_glue("{sample_id}_{seurat_clusters}")) %>% 
    select(CB, id_cluster)
  
  seurat_obj <- AddMetaData(seurat_obj, new_identities)
  Idents(seurat_obj) <- seurat_obj@meta.data$id_cluster
  
  seurat_obj
}

#' Generate all pairs of clusters where first pair component always lower than 
#' second component
#'
#' for 3 clusters (0,1), (0,2), (1,2) -- n(n-1)/2 - pairs
#'
#' @param seurat_obj 
#'
#' @return list of all pairs
#' @export
#'
#' @examples
pairwise_cluster_list <- function(seurat_obj){
  all_clusters <- cross(list(c1 = levels(seurat_obj), c2 = levels(seurat_obj))) %>% 
    keep(., function(x) {
      n1 <- str_extract(x$c1, "(?<=_)([[:alnum:]]+)")
      n2 <- str_extract(x$c2, "(?<=_)([[:alnum:]]+)")
      s1 <- str_extract(x$c1, "([[:alnum:]]+)")
      s2 <- str_extract(x$c2, "([[:alnum:]]+)")
      x$c1 != x$c2 && s1==s2 && n1 < n2
    })
  
  all_clusters
}

#' Find all markers by pairwise comparisons all vs all clusters
#'
#' @param seurat_obj Seurat object
#' @param features Vector of genes
#' @param logfc.threshold Log2FC threshold
#' @param min.pct Minimal percent of expression in cell threshold
#' @param topn How many top genes to return for each comparison (upper limit)
#'
#' @return Dataframe with \code{topn} expressed genes for each cluster
#' @export
#'
#' @examples
#' markers <- FindAllMarkersPairwiseCustomMultiple(seurat_object, debug = T)
FindAllMarkersPairwise <- function(seurat_obj, 
                                   features = NULL, 
                                   seurat_clusters=NULL,
                                   #logfc.threshold = 0.25, 
                                   #min.pct = 0.1,
                                   debug = F,
                                   #topn = 50,
                                   supp_id = NULL,
                                   parallel = F) 
{
  # extract_n_top_bottom <- function(t,n){
  #   rbind(slice_head(t,n = n), slice_tail(t, n = n)) %>% 
  #     distinct()
  # }
  
  `%uniform_do%` <- `%do%`
  if (parallel) {
    ## TODO: register parallel number of cores
    `%uniform_do%` <- `%dopar%`
  }
  
  ## if cluster does not specified just make cross of all existed clusters
  if (is.null(seurat_clusters)) {
    seurat_clusters <- cross(list(c1 = levels(seurat_obj), c2 = levels(seurat_obj)), .filter=function(x,y){x>=y})
  }

  #print(length(seurat_clusters))
  output <- foreach::foreach(i = seurat_clusters, .combine = rbind) %uniform_do% #%do% #%dopar% 
    {
      name <- paste0(i$c1,'vs', i$c2)
      if (debug) { name %>% print }
      FindMarkersLoupe(seurat_obj, id.1 = i$c1, id.2 = i$c2) %>% 
        #arrange(avg_log2FC) %>% 
        filter(abs(log2_fold_change) > 1 & ((id.1.intensity+id.2.intensity)/2) > 0.0000001 & adjusted_p_value < 1) %>% 
        # extract_n_top_bottom(., n = topn) %>% 
        mutate(id = name, c1 = i$c1, c2=i$c2)
        # %>% tibble::rownames_to_column(var = "gname")
    }
  
  output_mirror <- output %>% 
    mutate(tmp = c1, c1 = c2, c2 = tmp,
           #avg_log2FC = -avg_log2FC,
           log2_fold_change = -log2_fold_change,
           tmp_pct = id.1.intensity,
           id.1.intensity = id.2.intensity,
           id.2.intensity = tmp_pct) %>% 
    select(-tmp, -tmp_pct)
  
  tmp <- bind_rows(output, output_mirror)
  
  if (!is.null(supp_id)) {
    tmp$supp_id <- supp_id
  }
  
  tmp
}

#' Extract significant markers using sumrank
#'
#' @param df - (data.frame) - output from FindMarkersLoupe function for each 
#' pair of clusters
#' @param sample_id - (chr) 
#' @param topgenes - (int) number of top genes extracted for each cluster-cluster comparison
#' @param adj_pval_thr - (float) threshold by adjuster pvalue
#' @param log2fc_thr - (float) threshold by log2fc
#'
#' @return list with two data.frames
#' initial - all_filtered_genes - data.frame filtered by thresholds
#' final - all_filtered_genes.final - data.frame with sumranks
#' @export
#'
#' @examples
extract_significant <- function(df, sample_id, topgenes = 200, adj_pval_thr = 0.001, log2fc_thr = 0){ 
  
  all_filtered_genes <- df %>% filter(str_detect(c1, sample_id)) %>% 
    filter(adjusted_p_value < adj_pval_thr, log2_fold_change > log2fc_thr) %>% 
    mutate(id = str_extract(c1, "([[:alnum:]]+)"),
           c1 = str_extract(c1, "(?<=_)([[:alnum:]]+)"),
           c2 = str_extract(c2, "(?<=_)([[:alnum:]]+)")) %>% 
           c1 = str_extract(c1, "(?<=_)([[:print:]]+)"),
           c2 = str_extract(c2, "(?<=_)([[:print:]]+)")) %>% 
    group_by(c1,c2) %>% 
    slice_min(order_by = adjusted_p_value, n = topgenes) %>% ## extracting top N for each comparison
    arrange(adjusted_p_value, .by_group = T) %>% 
    #mutate(rank = 21-20*row_number()/n()) ## range from 20 to 1, if lines less than 20 gradually decrese rank
    mutate(rank = topgenes+1-row_number()) ## range from 20 to 1, if lines less than 20 gradually decrese rank
  
  all_filtered_genes.sumranks <- all_filtered_genes %>% 
    group_by(gname,c1) %>% 
    summarise(sumrank = sum(rank))
  
  all_filtered_genes.final <- all_filtered_genes %>% 
    group_by(gname,c1) %>% 
    summarise(comparison_power = n()) %>% ## number of c2 clusters which told us that this gene was significant (by pvalue and log2fc)
    left_join(., all_filtered_genes.sumranks, by = c("gname","c1")) %>%  ## attaching ranks
    ungroup() %>% group_by(c1) %>% 
    arrange(desc(sumrank), desc(comparison_power), .by_group = T) %>% 
    #slice_head(n = 20) %>% ## if it is required to trim final genes
    ungroup() %>% 
    group_by(gname) %>% 
    mutate(significant_for_clusters = n()) ## in how many clusters this specific gene significant
  
  
  all_filtered_genes.final # %>% rename(cluster1 = c1)
  list(final = all_filtered_genes.final, initial = all_filtered_genes)
}


## process all samples presented separately in splitted list
#' Using Seurat's Dotplot function for filtration genes by average expression
#' If average expression for some specific gene less than thresh_avgexp 
#' then drop this gene from markers
#'
#' @param list_seurat_split - list of Seurat object splitted by sample_id
#' @param list_final_markers - list final markers for specific sample with
#' @param thresh_avgexp 
#'
#' @return 
#' @export
#'
#' @examples
generic_filtering_by_avgexp <- function(list_seurat_split, list_final_markers, thresh_avgexp = 0.1){
  
  ## filter by average expression specific cluster in specific sample
  filter_by_avgexp <- function(sobj, df_markers, cl1, sample_id, thresh_avgexp = 0.1){
    tmp <- DotPlot(sobj, features = df_markers %>% filter(c1 == cl1) %>% pull(gname))$data %>%
      select(gname = features.plot, pct.exp, c1 = id, avg.exp.scaled) %>%
      mutate_if(is.factor,as.character) %>% 
      filter(c1 == cl1 & avg.exp.scaled > thresh_avgexp) %>%
      mutate(c1 = str_glue("{sample_id}_{cl1}"))
    tmp
  }
  
  ## process all clusters in specific sample
  filter_by_avgexp_onesample <- function(sobj, df_markers, sample_id, thresh_avgexp = 0.1) {
    map_dfr(levels(Idents(sobj)), function(cl_id){
      filter_by_avgexp(sobj, df_markers = df_markers, cl_id, sample_id, thresh_avgexp)
    })
  }

  sample_names <- names(list_seurat_split)
  
  map_dfr(sample_names, function(sample_id) {
    filter_by_avgexp_onesample(list_seurat_split[[sample_id]], list_final_markers[[sample_id]], sample_id, thresh_avgexp = thresh_avgexp)
  })
}

## for union/intersect of significant markers
## input - list [G193M1, G193M3,..], cluster_id
multiple_sample_markers <- function(list_of_markers, cluster_id, ntop = 30) {
  #sample_cluster_id <- str_glue("{sample_id}_{cluster_id}")
  markers_for_all_samples_specific_clusters <- map(names(list_of_markers), function(sample_id) {
    list_of_markers[[sample_id]] %>% filter(c1 == cluster_id) %>% head(ntop) %>%  pull(gname)
  })
  
  sig_union <- reduce(markers_for_all_samples_specific_clusters, union)
  sig_intersect <- reduce(markers_for_all_samples_specific_clusters, intersect)
  
  list(sig_union = sig_union,
       sig_intersect = sig_intersect)
}

# ## --- working example
# 
# sobject_name <- ""
# sobject <- add_id_cluster(sobject)
# all_clusters <- pairwise_cluster_list(sobject)
# system.time(res <- FindAllMarkersPairwise(seurat_obj = sobject, seurat_clusters = all_clusters))
# write_csv(res, str_glue("{sobject_name}_pairwise_clusters.csv"))
# res <- read_csv("{sobject_name}_pairwise_clusters.csv", col_names = T)
# 
# Idents(sobject) <- sobject@meta.data$seurat_clusters
# sobject_samples <- SplitObject(sobject, split.by = "sample_id")
# 
# for (s in names(sobject_samples)) {
#   sobject_samples[[s]] <- NormalizeData(sobject_samples[[s]], assay = "RNA")
# }
# 
# prelim_list_final_markers <- map(names(sobject_samples), function(sample_id){ 
#   temporary <- extract_significant(res, sample_id)
#   temporary$final
# }) %>% 
#   set_names(names(sobject_samples))
# 
# ## last parameter (0.2) threshold by avg.exp 
# filtered_markers <- generic_filtering_by_avgexp(sobject_samples, prelim_list_final_markers, 0.2)
# 
# ## drop markers which is not markers at all and has low avg.expr for specific cluster
# final_res <- left_join(res, filtered_markers, by=c("gname","c1")) %>% 
#   drop_na()
# 
# filtered_adjusted_markers <- map(names(sobject_samples), function(sample_id){ 
#   temporary <- extract_significant(final_res, sample_id)
#   temporary$final
# }) %>% 
#   set_names(names(samples))
# 
# 
# ntop_markers <- 60
# newplots <- map(names(sobject_samples), function(sample){
#   map(0:9,{function(clst){
#     title <- str_glue("Cluster {clst} ({sample})")
#     
#     # sort genes here
#     filtered_genes <- adjusted_markers[[sample]] %>% filter(c1 == clst) %>% pull(gname) %>% head(ntop_markers) ## <--extract top N markers
#     
#     # arrange by avg.expr * expression percent
#     sorted_genes <- DotPlot(sobject_samples[[sample]], features = filtered_genes)$data %>% 
#       filter(id == clst) %>% 
#       arrange(desc(avg.exp.scaled*pct.exp)) %>%
#       pull(features.plot)
#     
#     plot_dotplot_for_any_cluster(sobject_samples[[sample]], genes = sorted_genes, title = title)
#   }})
# }) %>% set_names(., names(sobject_samples))
# 
# 
# ## save pdf for each sample (all cluster for sample)
# ## 1 pdf = 1 sample all clusters
# map(names(newplots), function(sample){
#   output_filename <- str_glue("{sample}_{ntop_markers}top_markers_individual samples.pdf")
#   newplots[[sample]] %>% 
#     map(function(p){plot_grid(wrap_plots(p))}) %>% 
#     marrangeGrob(nrow = 4, ncol = 1) %>%
#     ggsave(filename = output_filename, width = 12, height = 18)  
# })
# 
# 
# 
# Idents(sobject) <- sobject@meta.data$seurat_clusters
# 
# #multiple_sample_markers(list_final_markers1, "9")
# 
# 
# ## create dotplots for union/intersection of top markers for each sample
# ## only for 
# union_intersect_dotplots <- map(levels(Idents(sobject)), function(cluster) {
#   ## how many top genes to extract
#   ntop <- 30
#   
#   ## extract data for specific cluster only
#   tst_cells <- WhichCells(sobject, idents = cluster)
#   tst_seu <- subset(sobject, cells = tst_cells) %>% 
#     NormalizeData(., assay = "RNA")
#   
#   mrks <- multiple_sample_markers(filtered_adjusted_markers, cluster, ntop = ntop)
#   
#   Idents(tst_seu) <- tst_seu@meta.data$id_cluster
# 
#   if (length(mrks$sig_intersect) == 0) {
#     (mp_dotplot(tst_seu, genes = mrks$sig_union)+ 
#        scale_colour_gradient2(low = "red", mid = "white", high = "blue"))&
#       plot_annotation(title = str_glue("Union of {ntop} most significant genes from each sample"))
#   } else {
#     (mp_dotplot(tst_seu, genes = mrks$sig_intersect)+ 
#        scale_colour_gradient2(low = "red", mid = "white", high = "blue")) / 
#       (mp_dotplot(tst_seu, genes = mrks$sig_union)+
#          scale_colour_gradient2(low = "red", mid = "white", high = "blue"))&
#       plot_annotation(title = str_glue("Top intersection / bottom union of top {ntop} most significant genes from each sample"))
#   }
# })
# 
# # union_intersect_dotplots[1]

