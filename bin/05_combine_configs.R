#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Combine configs')
    p <- add_argument(p, '--samples', default="samples.csv",
                      help='file with initial samples configuratation (sample_id, batch, etc columns )')
    p <- add_argument(p, '--modified_h5', default="h5_modified_list.csv",
                      help='files with modified molecule_info.h5 files (added extra cell barcodes) ')
    p <- add_argument(p, '--output', default="aggregated.csv",
                      help='sample id')
    
    return(parse_args(p))
}

argv <- ParseArguments()

library(dplyr)
library(readr)

samples <- read_csv(argv$samples, col_names = T) %>%
    select(-fastq_prefix,-fastqdir,-chemistry)

modified_h5 <- read_csv(argv$modified_h5, col_names = F) %>%
    select(sample_id = X1, molecule_h5 = X2)

full_join(samples, modified_h5, by="sample_id") %>%
    select(sample_id, molecule_h5, everything()) %>%
    distinct() %>%
    write_csv(argv$output)

