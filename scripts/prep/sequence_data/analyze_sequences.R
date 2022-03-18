#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to perform some descriptive statistics on sequences. 
#
# Code developed by Anna Hake
# Date: 2019-05-08
# Update:
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
  make_option(c("-i", "--input_file"), type="character", default=NULL, help="a fasta file with the sequences to analyze", metavar="FILE"), 
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="a directory where all results are stored", metavar="FILE"), 
)
# ------------------------------------------------------------------------------
# FUNCTIONS

## Write down your functions here. If more than five, outsource them in a separate file and use source
## source(".R")

#' Short title.
#'
#' @description Description.
#'
#' @param arg_name description of arg_name
#' @param arg2
#'
#' @return Description
#' 
#' @examples
#'
#' @export

#function_name<-function(arg1){
#  # assertive tests for each argument
#  
#  # missing values/coercion handling
#  
#  # warnings
#  
#  # code
#  

#}

# Roxygen documentation
#' Get the outliers based on the p-th percentile z-scores.
#'
#' @description Description.
#'
#' @param samples a vector of values.
#' @param prob the probability denoting the p-th percentile 
#'
#' @return The indices of the samples that exceed the p-th percentile based on the z-scores. 
#' 
#' @examples
#'
#' @export

get_char_freq_outliers<-function(samples, prob){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  is_outlier<-outliers::scores(samples, prob=prob)
  outlier_index<-which(is_outlier)
  return(outlier_index)

}

#' Retrieve outliers beyond p-th percentile based on z-scores
#'
#' @description Description.
#'
#' @param df data frame with the character frequencies per sample
#' @param prob The probabilty denoting the p-th percentile to look at
#'
#' @return A list with the indices of the outliers for each character.
#' 
#' @examples
#'
#' @export

get_char_freqs_outliers<-function(df, prob){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  outliers<-lapply(df, get_char_freq_outliers, prob=prob)
	return(outliers)
}

#' Computes the character frequencies in each sequence.
#'
#'
#' @param sequences character vector with the sequences
#' @param kmer integer denoting the size of the kmer to count for
#'
#' @return A data frame with all occurring characters and their occurrence 
#' frequency in each sequence
#' 
#' @examples
#'
#' @export

get_char_freqs<-function(sequences, kmer){
  # assertive tests for each argument

  # TODO: sequences should be a character vector
  # TODO: kmer should be an integer
  # TODO: is there an upper range for kmer?
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  char_freqs_list<-lapply(sequences, function(seq){
    # assert that the sequence is in character form
    char_seq<-s2c(seq)
    # compute length
    char_freqs<-data.frame(rbind(count(char_seq,kmer, alphabet=amb())))
  })
  char_freqs_df<-do.call(rbind,char_freqs_list)
  return(char_freqs_df)
}
#' Transform sequence list to matrix.
#'
#' @description Description.
#'
#' @param sequences list with sequences

#' @return A data frame with the sequences and their names.
#' 
#' @examples
#'
#' @export

transform2mat<-function(sequences){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  
  # get name of sequences
  seqs_names<-names(sequences)
  # convert sequences to string 
  sequences<-sapply(sequences, c2s)
  # create data frame
  seqs_mat<-data.frame(ID=seqs_names, sequence=sequences, stringsAsFactors=FALSE)
  # return data frame
  return(seqs_mat)
}



#' Get the length of sequences.
#'
#' @description Counts the number of nucleotides/amino acids in each sequence.
#'
#' @param sequences a character vector with the sequences
#'
#' @return A vector with the length of each sequence. 
#' 
#' @examples
#'
#' @export

get_sequences_length<-function(sequences){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  # count characters
  seqs_length<-sapply(sequences, function(seq){
    # assert that the sequence is in character form
    char_seq<-s2c(seq)
    # compute length
    char_seq_length<-length(char_seq)
  })
  return(seqs_length)
}

#' Plots the sequences distribution.
#'
#'
#' @param df data frame with the length of the sequences as column. 
#'
#' @return plot object
#' 
#' @examples
#'
#' @export

plot_sequence_distribution<-function(df){
  # assertive tests for each argument

  # missing values/coercion handling

  # warnings
  
  # code
  
  p<-ggplot(data=df, aes(df$length)) + geom_histogram()
  return(p)
}
# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-optparse::OptionParser(option_list=option_list)
# read arguments
opt<-optparse::parse_args(opt_parser)
# read in sequences
seqs<-seqinr::read.fasta(opt$input_file)
# transform seqs to matrix object
seqs_df<-transform2mat(seqs)
# check length
seqs_df$length<-get_sequences_length(seqs_df$sequence)
# get length outlier
length_outlier<-which(outliers::scores(seqs_df$length, prob=opt$outlier_cutoff))
# save plot
#pdf(opt$output_sequence_length)
## plot length histogram
#p<-plot_sequence_distribution(seqs_df)
#dev.off()

# check character distributions
char_freqs_df<-get_char_freqs(seqs_df$sequence, kmer=1)
# get outlier according to z score
char_freqs_outliers<-get_char_freqs_outliers(char_freqs_df)
write.csv(char_freqs_df, row.names=F, opt$output_char_freqs_df)
# check dimer frequencies
# dimer_freqs_df<-get_char_freqs(seqs_df$sequence, kmer=2)
# write.csv(dimer_freqs_df, row.names=F, opt$output_dimer_freqs_df)




