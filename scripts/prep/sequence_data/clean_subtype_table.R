#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to clean the names of the output of the COMET Tool. 
# Currently the COMET Tool trims any special character from the fasta file as 
# '=', '%' etc. Very specific code for the HIVIA project. not generalized. 
# Code developed by Anna Hake
# Date: 2019-05-17
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
  make_option(c("-s", "--subtype_table"), type="character", default=NULL, help="file with the subtypes from the COMET Tool", metavar="FILE")# TODO: Add comma if more than one option is used. 
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


#' Splits the character names to the subject ids. 
#'
#'
#' @param names vector of characters
#'
#' @return vector with same lenth as input containing the subject ids.
#' 
#' @examples
#'
#' @export

convert2subjectID<-function(name){
  # assertive tests for each argument
  
  # a single character of the following shape "subject_id[ID]|protein_name[protein_name]|consensus_cutoff[consensus_cutoff]|batch_id[batch_id]"
  # throw error if not in this form. 
  
  # missing values/coercion handling
  
  # warnings
  
  # code 	
  subject_id_items<-unlist(strsplit(name, split="\\|"))[1]
  subject_id<-unlist(strsplit(subject_id_items, "subject_id"))[2]
  return(subject_id)
}

#' Splits the character names to the subject ids. 
#'
#'
#' @param names vector of characters
#'
#' @return vector with same lenth as input containing the subject ids.
#' 
#' @examples
#'
#' @export

convert2subjectIDs<-function(names){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code 	
  new_names<-sapply(names, convert2subjectID)
  return(new_names)
}

# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# read in subtype table
subtype_df<-read.csv(opt$subtype_table, sep="\t", stringsAsFactors=F)
# generate new_names
#new_names<-convert2subjectIDs(subtype_df$name)
# set names to new names
#subtype_df$name<-new_names
# substitute commas within columns with /
subtype_df<-lapply(subtype_df, gsub, pattern=',', replacement='|')
# store table
write.csv(data.frame(subtype_df), gsub('.tsv', '_cleaned.csv',opt$subtype_table), row.names=FALSE)
