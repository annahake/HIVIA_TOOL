#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to
#
# Code developed by Anna Hake
# Date:
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", 'rentrez', 'seqinr')

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--id"), type="character", default=NULL, help="identifier of the reference sequence", metavar="string"), 
  make_option(c("-s", "--start"), type="integer", default=NULL, help="start position of the gene", metavar="int"),
  make_option(c("-e", "--end"), type="integer", default=NULL, help="end position of the gene", metavar="int"),
  make_option(c("-r", "--ref_seq_name"), type="character", default=NULL, help="name of the reference sequence", metavar="string"),
  make_option(c("-n", "--gene_name"), type="character", default=NULL, help="name of the gene", metavar="string"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output filename", metavar="FILE"),
  make_option(c("-d", '--database'), type="character", default=NULL, help="database", metavar="string")
)

# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# fetch the reference sequence
ref_seq<-rentrez::entrez_fetch(db=opt$database, id=opt$id, rettype='fasta')

# read the text object to fasta sequence
seq<-seqinr::read.fasta(textConnection(ref_seq))
switch(opt$database, 
nucleotide={gene<-seq[[1]][opt$start:opt$end];name<-paste(names(seq[1]), opt$gene_name, sep='_')},
protein={gene<-seq[[1]]; name<-names(seq[1])}
)
# extract sequence
# save fasta sequence to file
write.fasta(list(gene), names=name, file.out=opt$output)
