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
  make_option(c("-i", "--input"), default=NULL, help="Fasta file with codon-aligned amino acid test sequences to reference sequence C", metavar="FILE"), 
  make_option(c("-o", "--output"), default=NULL, help="Fasta file with codon-aligned amino acid test sequences without gaps", metavar="FILE")  
)
# -----------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# read input seqs
seqs<-read,fasta(opt$input, as.string=TRUE, seqtype ='AA')
# remove all gaps from sequences
seqs_new<-lapply(seqs, function(str){gsub('-', '', str)})
names(seqs_new)<-names(seqs)
# write sequences to file
write.fasta(seqs_new, opt$output, names = names(seqs_new))
