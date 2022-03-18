#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to subtype HIV-1 sequences with the COMET Tool 
# https://comet.lih.lu/
#
# Code developed by Anna Hake
# Date: 2019-05-06
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, help="fasta file with HIV-1 env nucleotide sequences", metavar='FILE'), 
    make_option(c("-o", "--output_file"), type="character", default=NULL, help="csv filename where the result will be stored", metavar="FILE"), 
    make_option(c('-a', '--api_source'), type='character', default=NULL, help="filename with R source code for the comet api", 
metavar='FILE')
)

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# ------------------------------------------------------------------------------
# FUNCTIONS

source(opt$api_source)


# ------------------------------------------------------------------------------
# MAIN


# subtype sequences
# TODO add get function
download_comet_subtypes(opt$input_file, opt$output_file)


