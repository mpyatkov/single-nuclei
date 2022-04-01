#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Extract the user-defined cellbarcodes for all samples from aggregated.csv ')
    p <- add_argument(p, '--aggregation', help="path to aggregation.csv",
                      default='aggregation.csv')
    p <- add_argument(p, '--userdefined_meta', help='raw h5 file counted in introns',
                      default="sample_paths.csv")
    p <- add_argument(p, '--output', default="aggr-user-defined-barcodes.csv",
                      help='output file')
    return(parse_args(p))
}

argv <- ParseArguments()

library(dplyr)
library(readr)
library(tidyr)

aggr_filename <- argv$aggregation
user_defined_metafiles_list <- argv$userdefined_meta
#aggr_filename <- "aggregation.csv"

# read aggr file
aggregation <- read_csv(aggr_filename, col_names = T) %>% 
    mutate(ord = row_number())

# read file with paths for user-defined-meta files
user_files <- read_csv(user_defined_metafiles_list, col_names = F) %>% 
    pull(X1)

# create output file with user-defined barcodes only and correct 'dashes' (-1,-2,-3,...)
# TODO it will be require probably to add additional columns with meta information to this table
output <- bind_rows(lapply(user_files, read_csv, col_names = T, show_col_types = FALSE)) %>% 
    separate(CB, c("CB", "tmp"), sep = "-") %>% 
    select(CB, sample_id) %>% 
    left_join(aggregation, by = "sample_id") %>% 
                  #select(sample_id, ord), by = "sample_id") %>%
    filter(!is.na(ord)) %>% 
    arrange(ord) %>% 
    mutate(CB = paste0(CB, "-",ord)) %>% 
    select(-molecule_h5,-ord) %>% 
    write_csv(argv$output, col_names = T)


