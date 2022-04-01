#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Inplace modifying molecule_info.h5 files')

    p <- add_argument(p, '--molecule_info_h5', default="molecule_info.h5",
                      help='molecule_info.h5 for appending cells')

    p <- add_argument(p, '--cellbarcodes', default="SAMPLEID_union_cellbarcodes.csv",
                      help='specific set of cellbarcodes we would like to use')
    
    return(parse_args(p))
}

argv <- ParseArguments()

library(hdf5r)
library(tidyverse)

## this function sets cell barcodes (CB) presented in meta_data as filtered
## by appending additional CB which was before in raw matrix 
replace_h5_cb_dataset <- function(h5file, meta_data) {
    library(rhdf5)
    ix <- h5read(h5file, "/barcode_idx")
    barcodes <- h5read(h5file, "/barcodes")
    passed <- h5read(h5file, "/barcode_info/pass_filter")
    
    barcodes_indexes <- tibble(CB = paste0(barcodes,"-1"), ix = seq_along(barcodes)-1)
    
    ## new barcodes
    new_passed <- barcodes_indexes %>%  
        inner_join(., meta_data %>% select(CB)) %>% 
        arrange(ix) %>% 
        select(ix) %>% 
        mutate(ix = as.integer(ix),
               x = as.integer(0), 
               y= as.integer(0)) %>% 
        as.matrix(., rownames=FALSE) %>% 
        t(.) %>% 
        as.array(new_passed, dim = dim(new_passed))
    
    ## delete dataset from h5 file
    h5delete(file = h5file, name = "/barcode_info/pass_filter")
    
    # create dataset and add it to file; useful links:
    # https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
    # https://support.hdfgroup.org/HDF5/doc1.6/PredefDTypes.html
    d = h5createDataset(file = h5file, dataset = "/barcode_info/pass_filter", dims = dim(new_passed), H5type="H5T_NATIVE_UINT64")
    h5write(new_passed, h5file, name = "/barcode_info/pass_filter")
    
    ## remove "rhdf5-NA.OK" attribute
    ## https://githubmemory.com/repo/grimbough/rhdf5/issues/84
    ## short way: h5deleteAttribute(h5File, name = "A", attribute = "rhdf5-NA.OK") # for new versions of rhdf5
    ## long way: 
    fid <- H5Fopen(h5file)
    did <- H5Dopen(fid, name = "/barcode_info/pass_filter")
    H5Adelete(did, "rhdf5-NA.OK")
    H5Dclose(did)
    H5Fclose(fid)
}

## read cellbarcodes
cellbarcodes<-read_csv(argv$cellbarcodes, col_names=T) %>%
    select(CB = 1)

##### Modifiyng molecule_info.h5 #####
## make a copy of the file in case we need to compare the original and modified files
##file.copy(argv$molecule_info_h5, paste0(argv$molecule_info_h5,".orig"))

## replace all CB in the original file to new one
replace_h5_cb_dataset(argv$molecule_info_h5, cellbarcodes)
