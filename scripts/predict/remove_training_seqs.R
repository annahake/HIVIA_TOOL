#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to
#
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", 'seqinr')

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--input"), default=NULL, help="Fasta file with aligned amino-acid test sequences to training sequences", metavar="FILE"), 
   make_option(c("-s", "--select"), default=NULL, help="Fasta file with test sequences", metavar="FILE"), 
     make_option(c("-o", "--output"), default=NULL, help="Fasta file with only aligned amino-acid test sequences", metavar="FILE")
)
# 
#------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# read in alignment
alignment <- read.fasta(opt$input)
# read in test seqs only fasta file
test_seqs<- read.fasta(opt$select)
# select only sequences from alignment which are also in test_seqs
alignment_new<- alignment[names(alignment) %in% names(test_seqs)]
# write reduced alignment file to fasta
write.fasta(alignment_new, file.out = opt$output, names = names(alignment_new))
