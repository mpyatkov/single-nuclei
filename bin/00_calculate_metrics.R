#!/usr/bin/env Rscript

## combine information from all metrics_summary for each of 4 gtf files
## args[1] path to pipeline SC output directory
## args[2] GEX or ATAC
## args[3] output file name

library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Calculate metrics by *summary.csv files')
  p <- add_argument(p,'--summary_path', default="./", help="path to directory with summary files (usually directory with raw output from cellranger count)")
  p <- add_argument(p,'--dataset_type', default="GEX", help="dataset_type could be 'GEX' or 'ATAC', because cellranger produces different summary files for different datasets")
  p <- add_argument(p, '--output', help='output file name', default="summary.csv")
  return(parse_args(p))
}

argv <- ParseArguments()
print(argv)

library(tidyverse)

## args = commandArgs(trailingOnly=TRUE)
## print(paste("Input dir: ", argv$summary_path))
## print(paste0("ATAC or GEX(default): ", argv$dataset_type)
## print(paste0("Output file: ", argv$output))

startpath <- argv$summary_path #"./SC_LIVER_TEST/output"
                               #startpath <- "./G190_SAMPLES/output/"

if (!dir.exists(file.path(startpath))) {
    print(paste0("Directory '",startpath,"' does not exist"))
    q("no",1)
}


pattern<-ifelse(argv$dataset_type == "ATAC", "summary.csv$", "metrics_summary.csv$")

metrics_files <- list.files(path=startpath, pattern = pattern, recursive = T) %>% 
    as_tibble() %>%
    mutate(basenames = dirname(value) %>% basename(),
           value = paste0(startpath,"/", value))  

print(metrics_files)
## q("no", 1)

## some columns for GEX dataset should be properly converted to percents
convert_pct_num <- function(x){
    tmp <- x %>% str_replace(.,"%","") %>% as.numeric()
    round(tmp/100, 5)
}

## table with nsamples * 4 rows
res <- metrics_files %>% 
    pmap_dfr(function(value, basenames){
        tmp <- read_csv(value,col_names = T) %>% 
            mutate(sample_id = basenames)
        
        if (argv$dataset_type == "GEX") {
            ## for GEX convert some columns from percents to real values
            tmp<- tmp %>% 
            rowwise() %>% 
            mutate(across(c(5:17), ~convert_pct_num(.)))
        }
        
        tmp
    }) %>% 
    relocate(sample_id, .before = everything())

res %>% write_csv(argv$output)


print("DONE")
