#!/usr/bin/env Rscript

## TMP: original location FindMarkersLoupe, mm10 converter and Utils.R: /projectnb/wax-dk/max/RSRC/G190

DEBUG <- FALSE

library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Extraction/Combining')
  p <- add_argument(p,'--input_rds', default="input_rds.rds", help="load precalculated R object")
  p <- add_argument(p, '--de_config', help='configuration file with a header. Contains info which samples/clusters should be compared to each other in RDS', 
                    default="de_config.csv")
  return(parse_args(p))
}

argv <- ParseArguments()
print(argv)

## install lab libs
library(devtools)
devtools::install_github("mpyatkov/FindMarkersLoupe", upgrade = "never")
devtools::install_github("mpyatkov/NotationConverter", upgrade = "never")

# detach("package:FindMarkersLoupe", unload=TRUE)

library(FindMarkersLoupe)
library(NotationConverter)

library(tidyverse)
library(gridExtra)
library(grid)
library(Seurat)
library(harmony)
library(patchwork)
library(cowplot)
source("UmapPlot.R")


if (DEBUG) {
 argv$de_config <- "/projectnb2/wax-dk/max/SCexp/G190_G183/configuration/de_config.csv" 
 argv$input_rds <- "/projectnb/wax-dk/max/SCexp/G190_G183/output/module_3_outputs/extract_combine/rds/output.rds"
}

# load rds file
input_seurat <- readRDS(argv$input_rds)


## load de_config
de_config <- read_csv(argv$de_config, col_names = T, show_col_types = FALSE) %>% 
  distinct()

##### validate_config
## AGGR,SAMPLE1,1,1 - same cluster if we compare AGGR vs SAMPLE does not have any sense
## ANY,ANY,1,1 - same sample/cluster pair does not have any sense

validation <- function(s1,s2,c1,c2) {
  
  init_row <- tibble(SAMPLE_ID_1 = s1, SAMPLE_ID_2 = s2, CLUSTER_ID_1 = c1, CLUSTER_ID_2 = c2)
  
  if (("AGGR" %in% c(s1,s2)) && c1 == c2) {
    res <- bind_cols(init_row, tibble(status = "FAIL", message = str_glue("(config error here: {s1},{s2},{c1},{c2}) - If one sample is AGGR then clusters CLUSTER_ID_1 ({c1}) and CLUSTER_ID_2 ({c2}) should not be the same")))
    return(res)
  } 
  if (s1 == s2 && c1 == c2) {
    res <- bind_cols(init_row, tibble(status = "FAIL", message = str_glue("(config error here: {s1},{s2},{c1},{c2}) - It does not make any sense to compare the same sample/cluster pair")))
    return(res)
  }
  
  bind_cols(tibble(SAMPLE_ID_1 = s1,
                   SAMPLE_ID_2 = s2, 
                   CLUSTER_ID_1 = c1, 
                   CLUSTER_ID_2 = c2), 
            tibble(status = "OK", message = "OK"))
}

## add status and message to original configuration file
de_config_validated <- pmap_dfr(de_config, ~ validation(..1, ..2, ..3, ..4)) %>% 
  left_join(de_config, .)

## quit if we have any "FAIL" status 
if (any(de_config_validated$status == "FAIL")) {
  ixs <- which(de_config_validated$status == "FAIL")
  map(ixs, function(ix){
    print(de_config_validated$message[ix])
  })
  if (!DEBUG){ q("no", 1, FALSE) }
}

## validation that config file and Seurat object have the same samples
all_sample_ids <- input_seurat@meta.data$sample_id %>% unique
user_defined_samples <- c(de_config_validated$SAMPLE_ID_1, de_config_validated$SAMPLE_ID_2) %>% 
  unique %>% 
  discard(`==`,"AGGR")
any_garbage <- setdiff(user_defined_samples,all_sample_ids)
if (!is_empty(any_garbage)) {
  print(str_glue("config error: Cannot find the following samples inside Seurat object: {any_garbage}"))
  if (!DEBUG){ q("no", 1, FALSE) }
}

## validation that config file and Seurat object have the same clusters id
all_cluster_ids <- input_seurat@meta.data$seurat_clusters %>% unique
user_defined_clusters <- c(de_config_validated$CLUSTER_ID_1, de_config_validated$CLUSTER_ID_2) %>% unique
any_garbage <- setdiff(user_defined_clusters,all_cluster_ids)
if (!is_empty(any_garbage)) {
  print(str_glue("config error: Cannot find the following clusters inside Seurat object: {any_garbage}"))
  if (!DEBUG){ q("no", 1, FALSE) }
}

if (DEBUG) {
  de_config_validated <- de_config_validated %>% filter(status != "FAIL")
}

## add index to each row
de_config_validated <- de_config_validated %>% 
  mutate(ix = row_number()) %>% 
  select(-status,-message)


#######################################################
## preparing sampleid_clusterid Idents for Seurat object
#######################################################

## Extract top/bottom n rows function
mp_extract_n_top_bottom <- function(t,n){
  rbind(slice_head(t,n = n), slice_tail(t, n = n)) %>% 
    distinct()
}

#' Calculate pairwise clusters comparisons
#'
#' @param seurat_obj 
#' @param s1 sample id 1 
#' @param s2 sample id 2
#' @param c1 cluster id 1
#' @param c2 cluster id 2
#'
#' @return list of the following objects:
#' - segex_output table with differential expression of clusters
#' - segex_filename output file for segex filename
#' - pdf with dotplot
#' 
#' @export
#'
#' @examples
compute_de <- function(seurat_obj, s1,c1,s2,c2, ix) {
  
  id.1 <- NULL
  id.2 <- NULL
  new_idents_df <- NULL
  id.1n <- str_glue("{s1}_{c1}")
  id.2n <- str_glue("{s2}_{c2}")

  ## create new Idents and use them 

  if (s1 == "AGGR" && s2 == "AGGR") {
    ## processing AGGR/AGGR
    ## using "seurat_clusters" as Idents
    
    new_idents_df <- input_seurat@meta.data %>% select(CB,new_idents = seurat_clusters) %>% 
      mutate(new_idents = case_when(
        new_idents == c1 ~ id.1n,
        new_idents == c2 ~ id.2n,
        TRUE ~ "Other"
      ))
    
    id.1 <- id.1n
    id.2 <- id.2n
    
  } else if (s1=="AGGR" || s2 == "AGGR") {
    
    ## processing of AGGR/sample or sample/AGGR
    ## using modified "sampleid_clusterid" as Idents
    ## preliminary prepare new Idents: AGGR_cluster_id and sampleid_clusterid
         
    ## need to specify which s1 or s2 are AGGR
    aggr_sid <- ifelse (s1 == "AGGR", s1, s2)
    aggr_cluster <- ifelse (s1 == "AGGR", c1, c2)

    not_aggr_sid <- ifelse (s1 == "AGGR", s2, s1)
    not_aggr_cluster <- ifelse (s1 == "AGGR", c2, c1)

    new_idents_df <- input_seurat@meta.data %>% 
      select(CB, sample_id, seurat_clusters) %>% 
      #rowwise() %>% 
      mutate(new_idents = case_when(
        sample_id == not_aggr_sid & seurat_clusters == not_aggr_cluster ~ glue::glue("{not_aggr_sid}_{not_aggr_cluster}"),
        seurat_clusters == aggr_cluster ~ glue::glue("AGGR_{aggr_cluster}"),
      TRUE ~ "Other"
      )) %>% 
      select(CB,new_idents) 
 
    id.1 <- id.1n
    id.2 <- id.2n
  }
  else {
    ## processing of sample/sample
    ## using "sampleid_clusterid"
    new_idents_df <- input_seurat@meta.data %>% select(CB, sample_id,seurat_clusters) %>% 
      mutate(new_idents = glue::glue("{sample_id}_{seurat_clusters}")) %>% 
      mutate(new_idents = case_when(
        sample_id == s1 & seurat_clusters == c1 ~ id.1n,
        sample_id == s2 & seurat_clusters == c2 ~ id.2n,
        TRUE ~ "Other"
      )) %>% 
      select(CB,new_idents) 

    id.1 <- id.1n
    id.2 <- id.2n
  }
  
  print(str_glue("ID1: {id.1}, ID2: {id.2}"))
  print(table(new_idents_df$new_idents))
  
  ## add new Idents to meta.data
  input_seurat <- AddMetaData(input_seurat, new_idents_df)
  
  ## activate new Idents
  Idents(object = input_seurat) <- input_seurat@meta.data$new_idents
  
  ## PROCESSING
  
  ## findMarkersLoupe (we cannot use the "short" Seurat object, because the entire 
  ## matrix must be used to calculate the intensities)
  markers_short <- FindMarkersLoupe(input_seurat, id.1 = id.1, id.2 = id.2, formatted = "short")
  
  if (DEBUG) {
    markers_short_seurat <- FindMarkers(input_seurat, ident.1 = id.1, ident.2 = id.2, logfc.threshold = 0) %>% 
      tibble::rownames_to_column(var="gname") 
  }

  ## prepare Segex output data.frame
  segex_output <- exportToSegex(input_df = markers_short)
  segex_fn <- str_glue("{ix}_scLoupe_{id.1}_vs_{id.2}_DiffExp_IntronicMonoExonic.tsv")

  
  ## DOTPLOT
  ## subset of only required clusters (need to create 2 lines DotPlot)
  only_clusters_CB <- WhichCells(input_seurat, idents = c(id.1, id.2))
  ## short 'Seurat' object only CB related to clusters
  only_clusters_seurat <- subset(input_seurat, cells = only_clusters_CB)
  
  ## get top 30 genes from both sides
  ## TODO: double plot top +/-30 genes and +/-30 lncRNA
  ## TODO: need to create histograms which show shift of pvalue if we compare two
  ## different is size clusters (small clusters will not have any significant 
  ## genes by pvalue, but some genes will have good Log2FC)

  top60genes <- markers_short %>% 
    #filter( 0.5*((id.2.intensity+1)+ (id.2.intensity+1)) > 1, !grepl("lnc", gname))
    filter(!grepl("lnc", gname)) %>% 
    arrange(desc(log2_fold_change)) %>% 
    mp_extract_n_top_bottom(., n = 30) %>% 
    pull(gname)
  
  top60lncrna <- markers_short %>% 
    filter(grepl("lnc",gname)) %>% 
    arrange(desc(log2_fold_change)) %>% 
    mp_extract_n_top_bottom(.,  n = 30) %>% 
    pull(gname)
  
  ncells <- table(new_idents_df$new_idents)
  
  id.1.ncells <- ncells[[id.1n]]
  id.2.ncells <- ncells[[id.2n]]
  
  cols <- rev(c("#225ea8","#6baed6","#eff3ff","#fbb4b9","#fbb4b9"))
  top60genes_dotplot <- wrap_elements(mp_dotplot(only_clusters_seurat, 
                                            top60genes, 
                                            title = str_glue("{id.1}({id.1.ncells} cells) vs {id.2}({id.2.ncells} cells). Top 30 genes with average intensity > 1")) + 
                                   scale_colour_gradientn(colors = cols))
  
  top60lncrna_dotplot <- wrap_elements(mp_dotplot(only_clusters_seurat, 
                                                  top60lncrna, 
                                                 title = str_glue("{id.1}({id.1.ncells} cells) vs {id.2}({id.2.ncells} cells). Top 30 lncRNA")) + 
                                        scale_colour_gradientn(colors = cols))
  
  top60_dotplot <- top60genes_dotplot/top60lncrna_dotplot+plot_layout(heights = c(5,4))
  
  ## back to usual Idents and clean meta.data
  Idents(object = input_seurat) <- input_seurat@meta.data$seurat_clusters
  input_seurat@meta.data$new_idents <- NULL
  
  if (DEBUG) {
    list(segex_output = segex_output$segex,
         lp = markers_short,
         sr = markers_short_seurat,
         segex_filename = segex_fn,
         pdf = top60_dotplot)
  } else {
    list(segex_output = segex_output$segex,
         segex_filename = segex_fn,
         pdf = top60_dotplot)
  }
}
 
# TODO: average intensity does not work need replacing for something else
# 

# z.t3 <- compute_de(input_seurat, "G190M2", 1, "G190M2", 2, 1)
# res <- left_join(z.t3$lp,z.t3$sr)
# z.t3$pdf

# if (DEBUG){
#   z.t1 <- compute_de(input_seurat, "AGGR", 1, "AGGR", 2, 1)
#   z.t2 <- compute_de(input_seurat, "AGGR", 1, "G190M2", 0, 1)
#   z.t3 <- compute_de(input_seurat, "G190M2", 1, "G190M2", 2, 1)
#   z.t4 <- compute_de(input_seurat, "G190M2", 1, "G183M1", 1, 1)
# }
# 
# View(z.t1$tmp)
# summary(sd(z.t2$tmp$adjusted_p_value))
# #zscore <- (z.t1$tmp$adjusted_p_value-median(z.t1$tmp$adjusted_p_value))/sd(z.t1$tmp$adjusted_p_value)
# zscore <- z.t2$tmp$adjusted_p_value
# zscore <- zscore[zscore != 1.0]
# hist(zscore, breaks = 100)

# z.t2 <- compute_de(input_seurat, "AGGR", 1, "G190M2", 0, 1)
# z.t1$pdf
# z.tmp <- z.t2$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# View(z.tmp)  
#
#
# z.t4 <- compute_de(input_seurat, "G190M2", 1, "G183M1", 1, 1)
# z.tmp <- z.t4$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# 
# View(z.t4$tmp)
#
# 
# View(z.t2$segex_output)
# z.t1 <- compute_de(input_seurat, "AGGR", 0, "AGGR", 2, 1)
# z.tmp <- z.t1$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# View(z.t1$tmp)
# View(z.tmp)


print("Saving Segex TSVs and prepare list of pdfs to save")
pdf_list <- pmap(de_config_validated, function(SAMPLE_ID_1, SAMPLE_ID_2, CLUSTER_ID_1, CLUSTER_ID_2,    ix) {
  
  res_list <- compute_de(input_seurat, SAMPLE_ID_1, CLUSTER_ID_1, SAMPLE_ID_2, CLUSTER_ID_2, ix)
  if (!DEBUG)  write_tsv(res_list$segex_output, res_list$segex_filename, col_names = T)

  res_list$pdf
})


print("Saving PDFs")
pdf_list %>% 
  map(function(p){plot_grid(wrap_plots(p))}) %>% 
  marrangeGrob(nrow = 3, ncol = 1) %>%
  ggsave(filename = "dotplots_by_comparision.pdf", width = 15.50, height = 20)





