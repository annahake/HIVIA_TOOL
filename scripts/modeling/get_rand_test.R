#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to retrieve the random HLA samples given a specific patient id. 
#
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", "dplyr")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="RDS object with the random data", metavar="FILE"), 
  make_option(c("-p", "--patient_id"), type="character", default=NULL, help="Character denoting the patient_id (for LOOCV)", metavar="STR"), 
  make_option(c("-o", "--output"), type="character", default=NULL, help="File where the result is stored", metavar="FILE"), 
  make_option(c("-t", "--test_df"), type="character", default=NULL, help="RDS object with the corresponding true data", metavar="FILE"), 
  make_option(c("--cv"), type="character", default=NULL, help = "Either 'cv' or 'loocv'. ")
  # TODO: Add comma if more than one option is used. 
)

# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# check if requied arguments are not empty
stopifnot(!is.null(opt$cv))
stopifnot(!is.null(opt$output))
stopifnot(!is.null(opt$input))

rand_df<-readRDS(opt$input)

switch(opt$cv, 
  cv={
    test_df<-readRDS(opt$test_df)
    stopifnot("patient_id" %in% colnames(rand_df))
    stopifnot("patient_id" %in% colnames(test_df))
    # select only the random samples which are in the current fold
    saveRDS(rand_df %>% filter(patient_id %in% test_df$patient_id), opt$output)
  }, 
  loocv={
    stopifnot(!is.null(opt$patient_id))
    patient_id_ <-opt$patient_id
    stopifnot("patient_id" %in% colnames(rand_df))
    saveRDS(rand_df %>% filter(patient_id == patient_id_), opt$output)
  },
  {stop("Error: unknown 'cv' option provided")}

)

