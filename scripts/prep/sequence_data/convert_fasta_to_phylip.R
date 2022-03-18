#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to convert fasta sequences into phylip format. 
#
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", "ape")

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="A fasta file.", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="A phylip file.", metavar="FILE"), 
  make_option(c("-f", "--format"), type="character", default="interleaved", help="Sequece format of phylip. Either 'interleaved' or 'sequential'", metavar="string") 
)
# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-optparse::OptionParser(option_list=option_list)
# read arguments
opt<-optparse::parse_args(opt_parser)
# read in fasta
fasta<-ape::read.dna(opt$input, format="fasta")
# write phylip format
ape::write.dna(fasta, file=opt$output, format=opt$format)
