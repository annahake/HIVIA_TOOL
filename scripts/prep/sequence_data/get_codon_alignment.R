#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to
#
# Code developed by Anna Hake
# Date:
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
  make_option(c("-i", "--input"), type="character", default=NULL, help="The file of the dna alignmnet or the alignment in text form", metavar="FILE/TEXT"), 
#
  make_option(c("-f", "--frame"), type="character", default='program', help="A character denoting the starting reading frame to chose. The frame is based on the first (master) sequence in the alignment. Possible options are '1', '2','3', 'all', 'program'."), 
#
  make_option(c("-c", "--compensating_num"), type="integer", default=5, help="an integer denoting how many codons are allowed for frameshift compensation. A compensating mutation is an insertion or deletion (indel) that can be combined with another nearby indel to produce a multiple of 3 and thus preserve an intact reading frame."), 
#
  make_option(c("-t", "--input_type"), type="character", default="file", help="A character denoting the type of alignment provided. Options are 'file' and 'text' ", metavar="string"),
#
  make_option(c("--output_format"), type="character", default="fasta", help="A character denoting the format of the output alignment. Options are 'fasta', 'table', 'mase', 'pretty', 'oal', 'rsf', 'gde', 'pir', 'gcg', 'clust', 'phylips', 'phylipi', 'msf', 'megas', 'megai', 'nexuss', 'nexusi', 'stockholm', 'slx'.",metavar="string"),
#
  make_option(c("--output_seqs_type"), type="character", default="all", help="A chacter denoting the type of the sequences in the alignment. Options are 'DNA', 'AA', or 'all'. If 'all' is selected, two alignments are downloaded.",metavar="string"), 
#
  make_option(c("--output_dir"), type="character", default=NULL, help="A character denoting the output directory. Note, that if the option 'all' is selected somewhere, several alignments will be downloaded and thus a directory is expected.",metavar="FILE/DIR"), 
#
  make_option(c("--api_source"), type="character", default=NULL, help="A file with the api functions",metavar="FILE"), 
#
  make_option(c("-n", "--protein_name"), type="character", default=NULL, help="A character denoting the protein name.",metavar="FILE")
)

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# ------------------------------------------------------------------------------
# FUNCTIONS

## Write down your functions here. If more than five, outsource them in a separate file and use source
source(opt$api_source)
 #------------------------------------------------------------------------------
# MAIN
download_codon_align(dna_alignment=opt$input, frame=opt$frame, compensating_num=opt$compensating_num, input_type=opt$input_type, output_format=opt$output_format, output_seqs_type=opt$output_seqs_type, id=opt$protein_name, output_dir=opt$output_dir)

