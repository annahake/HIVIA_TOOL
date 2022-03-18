#!/usr/bin/env Rscript
# 
# This is code to retrieve a specific sequence from the LANL align tool to be 
# used as reference sequence. 
# 
# Code developed by Anna Hake
# Date: 2019-05-02                                                              
# load libraries and functions
library("optparse")
library("seqinr")

# read arguments
option_list = list(
  make_option(c("-a", "--alignment_file"), type="character", default=NULL, 
              help="output file name for alignment", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="reference_sequence.fasta", 
              help="output file name for reference sequence [default= %default]", metavar="character"), 
  make_option(c("-t", "--alignment_type"), type="character", default=NULL, 
              help="alignment type", metavar="character"), 
  make_option(c("-s", "--suborganism"), type="character", default="HIV1", 
              help="suborganism [default= %default]", metavar="character"), 
  make_option(c("--subtype"), type="character", default="All", 
              help="subtype [default= %default]", metavar="character"),
  make_option(c("-d", "--region_definition"), type="character", default="predefined", 
              help="region definition to use [default= %default]", metavar="character"),
  make_option(c("-r", "--region"), type="character", default="reference_sequence.fasta", 
              help="output file name for reference sequence [default= %default]", metavar="character"), 
  make_option(c("-b", "--basetype"), type="character", default=NULL, 
              help="basetype [default= %default]", metavar="character"),
  make_option(c("-y", "--year"), type="character", default=NULL, 
              help="year of alignment [default= %default]", metavar="character"),
  make_option(c("-f", "--format"), type="character", default="fasta", 
              help="format of alignment [default= %default]", metavar="character"), 
  make_option(c("-q", "--query"), type="character", default=NULL, 
              help="Query token to retrieve reference sequence", metavar="character"), 
  make_option(c("--api_source"), type="character", 
              help="filename of api code", metavar="character")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# load api functions
source(opt$api_source)

# get alignment
if(!is.null(opt$alignment_file)){
  # download alignment
  download_alignment(alignment_type=opt$alignment_type, suborganism=opt$suborganism, region_definition=opt$region_definition, region=opt$region, subtype=opt$subtype, basetype=opt$basetype, year=opt$year, format=opt$format, output_file=opt$alignment_file)
  switch(opt$basetype, 
         "DNA"={seqtype<-"DNA"}, 
         "PRO"={seqtype<-"AA"}
         )
  align<-read.fasta(opt$alignment_file, seqtype = seqtype, as.string=TRUE)
} else {
  align<-get_alignment(alignment_type=opt$alignment_type, suborganism=opt$suborganism, region_definition=opt$region_definition, region=opt$region, subtype=opt$subtype, basetype=opt$basetype, year=opt$year, format=opt$format)
}


# get index of consensus C sequence
index<-which(names(align)==opt$query)

# store consensus C sequence
reference_seq<-align[[index]][[1]]

# remove gaps?
#reference_seq_wo_gaps<-gsub("-", "", reference_seq)

# write reference sequence to fasta
write.fasta(sequences=reference_seq, names=opt$query, as.string=TRUE, file.out=opt$out)

# save full alignment
write.fasta(sequences=align, names=names(align), as.string=TRUE, file.out=opt$alignment_file)



