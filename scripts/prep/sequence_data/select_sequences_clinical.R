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
required_packages<-list("optparse", 'seqinr', 'dplyr')

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Fasta file with DNA sequences", metavar="FILE"), 
  make_option(c("-c", "--clinical"), type="character", default=NULL, help="CSV with clinical data", metavar="FILE"), 
  make_option(c("-o", "--output"), type="character", default=NULL, help="Fasta file with selected sequences", metavar="FILE"),
  make_option(c("--leftout"), type="character", default=NULL, help="Fasta file with deselected sequences", metavar="FILE")
)
#  ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# read in sequence file
seqs<-read.fasta(opt$input, seqtype = 'DNA')
# read in clinical file
clin<-read.csv(opt$clinical, stringsAsFactors = FALSE)
## assert that patient_id and cd4_count are existing columns in the clinical file
stopifnot('cd4_count' %in% colnames(clin))
stopifnot('patient_id' %in% colnames(clin))
stopifnot('onART' %in% colnames(clin))
## ensure that clinical information exists for every sequence
stopifnot(all(names(seqs) %in% clin$patient_id))
selected_ids <- clin %>% filter(cd4_count <= 500, onART == 'NO') %>% select(patient_id) %>% pull
seqs_selected<-seqs[names(seqs) %in% selected_ids]
write.fasta(seqs_selected, names = names(seqs_selected), file.out = opt$output)
leftout_ids <-clin %>% filter(cd4_count > 500, onART == 'NO') %>% select(patient_id) %>% pull
seqs_leftout<-seqs[names(seqs) %in% leftout_ids]
write.fasta(seqs_leftout, names = names(seqs_leftout), file.out = opt$leftout)
traceback()
