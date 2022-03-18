#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to clean the amino acid alignment returned by the codon align 
# tool. Cleaning comprises the handling of stop signs, frameshift signs etc. 
# Code developed by Anna Hake
# Date: 2019-05-13
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", "seqinr", "xtable", "dplyr")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="a fasta file with an AA alignment", metavar="FILE"), 
  make_option(c("-r", "--refseq"), type="character", default=NULL, help="if given, the first sequence is treated as reference sequence, dropped from the alignment and stored under this filename", metavar="FILE"),
  make_option(c("-a", "--alignment"), type="character", default=NULL, help="a fasta file name for the cleaned AA alignment", metavar="FILE"), 
  make_option(c('-d', '--dist_matrix'), type='character', default=NULL, help="a csv file name for the distance matrix", metavar="FILE"),
  make_option(c("-s", "--seq_freq_table"), type="character", default=NULL, help="csv file name for the sequence frequency table", metavar="FILE"),
  make_option(c("-p", "--pos_freq_table"), type="character", default=NULL, help="csv file name for the position frequency table", metavar="FILE"), 
  make_option(c("-l", "--log"), type="character", default=NULL, help="log file", metavar="FILE")#,
  #make_option(c( "--phylo"), type="character", default=NULL, help="fasta file for the cleaned codon alignment", metavar="FILE") 
)
# ------------------------------------------------------------------------------
# FUNCTIONS

## Write down your functions here. If more than five, outsource them in a separate file and use source
## source(".R")

# Roxygen documentation
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


#' Count the frequency of every symbol in each sequence of the aa alignment.
#'
#' @description Description.
#'
#' @param alignment amino acid alignment
#'
#' @return a data frame with the frequencies for each sequence
#' 
#' @examples
#'
#' @export

count_freqs<-function(alignment){
  # assertive tests for each argument
  if(class(alignment)!="list"){
    stop("input alignment has to be of class 'list'.")
  }
  # missing values/coercion handling
  
  # warnings
  
  # code
  
  ## define alphabet
  aa_alphabet<-s2c("abcdefghiklmnpqrstvwxyz-*#")
  
  ## 
  freqs_list<-lapply(alignment, function(seq){
    oligomer<-1
    seqinr::count(seq, oligomer, alphabet=aa_alphabet)
  })
  freqs_df<-do.call(rbind, freqs_list)
  return(freqs_df)

}

#' Replace given characters with substitutes.
#'
#' @description Description.
#'
#' @param alignment a list of sequences of same length representing the alignment
#' @param char_vec a character vector of characters that have to be substituted
#' @param sub_vec a character vector of characters to replace the the char_vec\
#'
#' @return The alignment file with the replacements
#' 
#' @examples
#'
#' @export

substitute_chars<-function(alignment, char_vec, sub_vec){
  # assertive tests for each argument
  if (class(alignment)!="list"){
    stop("alignment has to be a list of sequences")
  }
  if (class(char_vec)!="character"){
    stop("char_vec has to be a vector of characters")
  }
  if (class(sub_vec)!="character" && !is.na(sub_vec)){
    stop("sub_vec has to be a vector of characters")
  }
  if (length(char_vec)!=length(sub_vec)){
    stop("replacement vector sub_vec has different length than input vector")
  }
  # missing values/coercion handling
  
  # warnings
  
  # code
  
  ## iterate over every sequence
  new_alignment<-lapply(alignment, function(seq_){
    # change all occurrences of the character list with the replacements in the sequence
    new_seq_<-seq_
    for (i in seq(length(char_vec))){
      new_seq_<-gsub(pattern=char_vec[i], replacement=sub_vec[i], x=new_seq_)
    }
    # last sequence is the most changed sequence
    return(new_seq_)
  })
  #
  return(new_alignment)
}
#' Checks if character is a single char or not.
#'
#'
#' @param char character
#'
#' @return A boolean denoting if character is a single char or not. 
#' 
#' @examples
#'
#' @export

is.single_character<-function(char){
  # assertive tests for each argument
  if (class(char)!='character'){
  	stop('Input to "is.single_character" must be character.')
  }
  # code
  length(unlist(strsplit(char, split=''))==1)
}

#' Converts gaps to NA.
#'
#' @description Converts leading and trailing gaps to NA. Removes entire gap columns.
#'
#' @param alignment list of sequences of the same length/
#'
#' @return alignment with leading and trailing gaps converted to NA and no gap columns. 
#' 
#' @examples
#'
#' @export

convert_gaps<-function(char_seq){
  # assertive tests for each argument
  if (class(char_seq)!='SeqFastadna'){
  	stop('Input to convert_gaps has to be a character vector')
  }
  if (!all(sapply(char_seq, is.single_character))){
  	stop('Input to convert_gaps has to be a character vector with single chars not strings.')
  }
  # missing values/coercion handling
  
  # warnings
  
  # code
  # definition leading gaps=sequences start with gap
  not_gaps_ind<-which(char_seq!='-')
  leading_end<-not_gaps_ind[1]-1
  trailing_start<-not_gaps_ind[length(not_gaps_ind)] + 1 
  if(leading_end>=1){
  	char_seq[1:leading_end]<-NA
  }
  if(trailing_start<=length(char_seq)){
  	char_seq[trailing_start:length(char_seq)] <-NA
  }
	return(char_seq)
}


# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# open log connection
sink(opt$log)
# read in alignment_file as fasta (no benefit in reading in as alignment)
alignment<-seqinr::read.fasta(opt$input)

# drop first sequence from alignment if this is the reference sequence and store it 
if (!is.null(opt$refseq)){
  ref_seq_aligned<-alignment[[1]]
  seqinr::write.fasta(list(ref_seq_aligned), file.out=opt$refseq,names=names(alignment)[1])
  alignment<-alignment[-1]
}

if (!is.null(opt$seq_freq_table)){
  # char frequency table
  char_freqs_df<-count_freqs(alignment)
  # if interested in percentage divide by number of positions
  write.csv(char_freqs_df, file=opt$seq_freq_table)
}

if (!is.null(opt$pos_freq_table)){
  # get column wise char frequency table -> create df where rows are not patients but alignment pos
  # TODO: needs improvement
  position_list<-pos_list<-lapply(as.data.frame(do.call(rbind, alignment)), function(x)x)

  # char frequency table
  char_freqs_df<-count_freqs(position_list)
  write.csv(char_freqs_df, file=opt$pos_freq_table)
}
frameshifts<-lapply(alignment, function(seq) {which(seq=="#")})
real_frameshifts<-sapply(frameshifts, function(pos){ any(pos!=1) })
real_frameshift_ids<-names(frameshifts)[real_frameshifts]
# if any real frameshifts
if(length(names(frameshifts)[real_frameshifts])>0){
  # create table with id and frameshift positions
  vec<-sapply(frameshifts[real_frameshifts], function(inds){paste(inds, collapse=",", sep="")})
  names(vec)<-real_frameshift_ids
  # print table to log (sink)
  print("Following sequences contain frameshifts.")
  xtable(tbl_df(vec), type="latex")
}


stop_signs<-sapply(alignment, function(seq){which(seq=="*") })
# look for stop signs that are not at the end of the sequence
unusual_stop_signs<-names(stop_signs)[sapply(alignment, function(seq){grepl("\\*(?!(\\-)*$)", c2s(seq), perl=TRUE) })]
# if any unusual stop signs are present:
if(length(unusual_stop_signs)>0){
  # collapse all indices to string
  vec<-sapply(stop_signs[unusual_stop_signs], function(inds){paste(inds, collapse=",", sep="")})
  # set names
  names(vec)<-unusual_stop_signs
  # print table
  print("Following sequences contain stop signs inside the sequence.")
  xtable(tbl_df(vec), type="latex")
}

# Exclude frameshifts and stop signs
# TODO: store in one excluded file?
#excluded_ids <- union(unusual_stop_signs, real_frameshift_ids)
#new_alignment<-alignment[-excluded_ids]
## convert some characters to create a valid alignment for a phylogentic tree
#phylo_alignment<-substitute_chars(alignment, c("#", "\\*", "x"), c("-","-", "-"))
#write.fasta(phylo_alignment, file.out=opt$phylo, names=names(alignment))

# substitute some characters
alignment_converted<-substitute_chars(alignment, c("#", "\\*", "x"), c("-","-", NA))
# convert leading and trailing gaps to NA
# CHANGED 20190712: gaps mean variant not present no matter if in sequence or outside the sequence. 
#alignment_cleaned<-lapply(alignment_converted, convert_gaps)
alignment_cleaned <- alignment_converted
#alignment mat
alignment_cleaned_df<-do.call(rbind, alignment_cleaned)

# save alignment
write.csv(alignment_cleaned_df, file=opt$alignment)
# Note, gap columns are not removed since then the mapping would be wrong wrt to reference sequence
final_alignment<-seqinr::as.alignment(nb=length(alignment_cleaned),nam=names(alignment_cleaned), seq=alignment_cleaned, com=NA)
# Not run! segmentation fault using dist mat
#if (!is.null(opt$dist_matrix)){
#  dist_alignment<-seqinr::dist.alignment(final_alignment, matrix="similarity")
#  write.csv(dist_alignment, file=opt$dist_alignment)
#}

# close log connection
sink()
