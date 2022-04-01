#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Intronic/Exonic plot')
    
    p <- add_argument(p, '--intronic_h5_filtered', default="filtered_feature_bc_matrix.h5",
                      help='filtered h5 file counted in introns')
    
    p <- add_argument(p, '--intronic_h5_raw', default="raw_feature_bc_matrix.h5",
                      help='raw h5 file counted in introns')
    
    p <- add_argument(p, '--exonic_h5_filtered', default="filtered_feature_bc_matrix.h5",
                      help='filtered h5 file counted in introns')
    
    p <- add_argument(p, '--exonic_h5_raw', default="raw_feature_bc_matrix.h5",
                      help='raw h5 file counted in introns')
    
    p <- add_argument(p, '--genebody_h5_filtered', default="filtered_feature_bc_matrix.h5",
                      help='filtered h5 file counted in genebody')
    
    p <- add_argument(p, '--genebody_h5_raw', default="raw_feature_bc_matrix.h5",
                      help='raw h5 file counted in genebody')
    
    p <- add_argument(p, '--intronic_withmono_h5_filtered', default="filtered_feature_bc_matrix.h5",
                      help='filtered h5 file counted in introns + mono exonic genes')
    
    p <- add_argument(p, '--intronic_withmono_h5_raw', default="raw_feature_bc_matrix.h5",
                      help='raw h5 file counted in introns + mono exonic genes')
    
    p <- add_argument(p, '--output_rds', default="default.rds", help='output filename for rds file')

    p <- add_argument(p, '--output_cellbarcodes', default="SAMPLEID_union_cellbarcodes.csv", help='output filename for rds file')
    
    p <- add_argument(p,'--debug', default = 'FALSE', help='add to RDS all raw files')
    
    return(parse_args(p))
}

argv <- ParseArguments()

print(argv)

DEBUG <- as.logical(argv$debug)

library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)

## by default we prefer to use select and filter from dplyr package
select <- dplyr::select
filter <- dplyr::filter

## read cellranger results from filtered dataset for intronic counting (without monoexonic genes)
read_h5_matrix <- function(matrix_path, matrix_name) {
    print(paste0("Loading: ", matrix_name," - " ,matrix_path))
    mx <- Read10X_h5(matrix_path)
    mx <- CreateSeuratObject(mx)
    mx@meta.data$samples <- matrix_name
    mx@meta.data$CB <- rownames(mx@meta.data)
    mx
}

intronic_filtered <- read_h5_matrix(argv$intronic_h5_filtered, "intronic_filtered")
intronic_raw <- read_h5_matrix(argv$intronic_h5_raw, "intronic_raw")

## add to raw data marker that this cell presented in filtered data (CB_intronic) PROBABLY IS NOT REQUIRED
intronic_raw@meta.data <- left_join(intronic_raw@meta.data, 
                                    intronic_filtered@meta.data %>% select(CB, CB_intronic = samples), by=c("CB")) 

exonic_filtered <- read_h5_matrix(argv$exonic_h5_filtered, "exonic_filtered")
exonic_raw <- read_h5_matrix(argv$exonic_h5_raw, "exonic_raw")

## add to raw data marker that this cell presented in filtered data (CB_exonic) PROBABLY IS NOT REQUIRED
exonic_raw@meta.data <- left_join(exonic_raw@meta.data, 
                                  exonic_filtered@meta.data %>% select(CB, CB_exonic = samples), by=c("CB")) 

withmono_filtered <- read_h5_matrix(argv$intronic_withmono_h5_filtered, "withmono_filtered")
withmono_raw <- read_h5_matrix(argv$intronic_withmono_h5_raw, "withmono_raw")
genebody_filtered <- read_h5_matrix(argv$genebody_h5_filtered, "genebody_filtered")
genebody_raw <- read_h5_matrix(argv$genebody_h5_raw, "genebody_raw")

## full join of the intronic and exonic raw matrices, select intronic and exonic UMIs
full_raw<- full_join(intronic_raw@meta.data, exonic_raw@meta.data, by=c("CB")) %>% 
    select(CB, intronic_umi = nCount_RNA.x, exonic_umi = nCount_RNA.y)

## get information about MT contamination from intronic-mono and genebody
withmono_raw[["percent.mt_withmono"]] <- PercentageFeatureSet(withmono_raw, pattern = "^(MT|Mt|mt)-")
withmono_raw_percent <- withmono_raw@meta.data %>% select(CB, percent.mt_withmono)
genebody_raw[["percent.mt_genebody"]] <- PercentageFeatureSet(genebody_raw, pattern = "^(MT|Mt|mt)-")
genebody_raw_percent <- genebody_raw@meta.data %>% select(CB, percent.mt_genebody)
exonic_raw[["percent.mt_exonic"]] <- PercentageFeatureSet(exonic_raw, pattern = "^(MT|Mt|mt)-")
exonic_raw_percent <- exonic_raw@meta.data %>% select(CB, percent.mt_exonic)
intronic_raw[["percent.mt_intronic"]] <- PercentageFeatureSet(intronic_raw, pattern = "^(MT|Mt|mt)-")
intronic_raw_percent <- intronic_raw@meta.data %>% select(CB, percent.mt_intronic)

## full join of all 4 filtered matrices obtained from 4 GTFs
full_filtered <- full_join(intronic_filtered@meta.data %>% select(CB,samples1 = samples), 
                           exonic_filtered@meta.data %>% select(CB,samples2 = samples), by=c("CB")) %>% 
    full_join(., withmono_filtered@meta.data %>% select(CB,samples3 = samples), by=c("CB")) %>% 
    full_join(., genebody_filtered@meta.data %>% select(CB,samples4 = samples), by=c("CB")) %>% 
    rowwise() %>% 
    mutate(summ=list(c(samples1,samples2,samples3,samples4)),
           summ2 = length(summ[!is.na(summ)])) %>% 
    ungroup()

# attach information about number of UMIs from raw matrix (make it log-scale) and MT percent infromation
only_filtered <- full_filtered %>% 
    left_join(., full_raw %>% select(CB, intronic_umi, exonic_umi), by=c("CB")) %>%
    left_join(., withmono_raw_percent, by=c("CB")) %>% 
    left_join(., genebody_raw_percent, by= c("CB")) %>% 
    left_join(., exonic_raw_percent, by=c("CB")) %>% 
    left_join(., intronic_raw_percent, by= c("CB"))

## downstream raw can be as withmono_raw, exonic_raw, intronic_raw or genebody_raw
only_filtered_4gtfs <- only_filtered %>% 
    mutate(rownames = CB) %>% 
    tibble::column_to_rownames(var = "rownames") %>% 
    select(CB)

only_filtered_4gtfs %>%
    select(Barcode = CB) %>%
    write_csv(argv$output_cellbarcodes, col_names = T)

if (!DEBUG) {
    ## we can use genebody/withmono as downstream files
    ## before exporting keep only union of CB from 4 gtfs files
    genebody_raw <- genebody_raw[,only_filtered_4gtfs$CB]
    withmono_raw <- withmono_raw[,only_filtered_4gtfs$CB]

    save(only_filtered,
         genebody_raw,
         withmono_raw,
         file = argv$output_rds)
} else {
    save(genebody_raw,
         intronic_raw,
         exonic_raw,
         withmono_raw,
         only_filtered, 
         downstream_seurat,
         downstream_seurat_name,
         file = argv$output_rds)    
}

## !!! withmono_raw and genebody_raw in rds contains only union of 4 GTF
## based cellbarcodes (not all of them)

# remove big objects which is not required anymore
# we keep withmono_raw because we use it for further analysis
## rm(genebody_raw, genebody_filtered,
##    intronic_raw, intronic_filtered,
##    exonic_raw, exonic_filtered,
##    withmono_raw, withmono_filtered,
##    genebody_raw_percent,
##    withmono_raw_percent)



