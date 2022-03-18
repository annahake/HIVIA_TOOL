#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to select sequences for the HIVIA analyses according their 
# clinical parameters
# Code developed by Anna Hake
# Date: 2020-02-04
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", 'seqinr', "dplyr", "tidyr")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=T))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-s", "--seq_file"), type='character', default=NULL, metavar="FILE"), 
  make_option(c("-c", "--clin_file"), type="character", default=NULL, help="a csv file with the clinical parameter", metavar="FILE"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, help="a fasta file", metavar="FILE"), 
  make_option(c('--cd4_cutoff'), type="integer", default="500", help="CD4 cutoff for patient selection", metavar="INT")
)

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)

# TODO change to fasta file (before alignment takes place)
#seq<-read.csv(opt$seq_file, stringsAsFactors=FALSE)
seq<-read.fasta(opt$seq_file)
#print(length(names(seq)))
clin_df<-read.csv(opt$clin_file, stringsAsFactors=FALSE)
pid_list <- clin_df %>% filter(onART =="NO", cd4_count<=opt$cd4_cutoff) %>% select(patient_id) %>% pull
#seq <- seq %>% filter(X %in% pid_list) 
seq<- seq[names(seq) %in% pid_list]
#print(length(names(seq)))
#write.csv(seq, row.names=FALSE, file=opt$output_file)

write.fasta(seq, names=names(seq), file.out = opt$output_file)
