#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to remove the reference sequence from the nt_codon_align file. 
#
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", "seqinr")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="The nt codon alignment file", metavar="FILE"), 
  make_option(c("-r", "--refseq_name"), type="character", default=NULL, help="The name of the reference sequence", metavar="STRING"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="The nt codon alignment file without the reference sequence.", metavar="FILE")
)
# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# read in alignment_cleaned
alignment<-read.fasta(opt$input, as.string=TRUE, seqtype="DNA")
# skip the first sequence (should be named CONSENSUS_C)
stopifnot(names(alignment)[1]==opt$refseq_name)
new_alignment<-alignment[-1]
# store alignment
write.fasta(new_alignment, opt$output, names=names(new_alignment))
