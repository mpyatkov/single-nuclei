#!/usr/bin/env Rscript

## Convert samples.xlsx to samples.csv for Nextflow pipeline
## Usage: Rscript convert_samples_xlsx.R samples.xlsx > samples.csv

library(readxl)
library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Convert samples.xlsx to samples.csv')
  p <- add_argument(p, 'input', help='Input XLSX file path')
  p <- add_argument(p, '--output', help='Output CSV file path (default: stdout)', default='-')
  return(parse_args(p))
}

argv <- ParseArguments()

# Read XLSX file
tryCatch({
  data <- read_excel(argv$input, sheet = 1)
}, error = function(e) {
  stop(paste("Error reading XLSX file:", e$message))
})

# Validate required columns
required_cols <- c('sample_id', 'chemistry', 'condition', 'path_to_r1')
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse=", ")))
}

# Write CSV
if (argv$output == '-') {
  write.csv(data, stdout(), row.names = FALSE)
} else {
  write.csv(data, argv$output, row.names = FALSE, quote = FALSE)
  cat(sprintf("Successfully converted %s to %s\n", argv$input, argv$output))
}
