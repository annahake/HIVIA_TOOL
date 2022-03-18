#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to select the complete cases of the modelling relevant features
# and outcome from the full data.
#
#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
required_packages<-list("optparse", "dplyr")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
# ------------------------------------------------------------------------------
option_list=list(
  make_option(c("-i", "--input"), type="character", default=NULL,help="RDS object with the SAAV data", metavar="FILE"),
  make_option(c("-f", "--fixed"), type="character", default=NULL, help="a character vector with fixed variables", metavar="STRING"),
  make_option(c("-r", "--random"), type="character", default=NULL, help="a character vector with random variables", metavar="STRING"),
  make_option(c("-y", "--outcome"), type="character", default=NULL, help="a character vector with outcome variables", metavar="STRING"),  
  make_option(c("--hla"), type="character", default=NULL, help="a character vector with hla alleles to consider", metavar="STRING"),  
  make_option(c("-o", "--output"), type = "character", default=NULL, help="Output filename", metavar="FILE"), 
  #make_option(c("-l", "--log"), type = "character", default=NULL, help="Log filename", metavar="FILE"),
  make_option(c("--id"), type = "character", default=NULL, help="The SAAV id column", metavar="FILE")
)



# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------
#' Converts list arguments ("X,X,X") to vector of characters.
#'
#'
#' @param char character representing the list argument
#'
#' @return vector of characters where each element represents element of the 
# argument list
#' 
#' @examples
#'
#' @export

convert_arglist <-function(char){
	print(char)
  stopifnot(!is.null(char))
  stopifnot(!is.na(char))
  stopifnot(is.character(char))
  char_vec<-unlist(strsplit(char, split=""))
  if ("," %in% char_vec){
    args<-unlist(strsplit(char, split="\\,"))
  }
  else {
    if(char=="NULL"){
      args<-NULL
    } else {
      args <-char
    }
  }
  return(args)

}





# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

## SETTING OPTION PARSER

opt_parser<-OptionParser(option_list=option_list)

# read arguments
opt<-parse_args(opt_parser)

## OPEN LOG CONNECTION
#if (!is.null(opt$log)){
#tryCatch(sink(opt$log), error=function(e){
#    e$message <- paste0("Error: sink connection could not be opened with the following error msg: \n", e$message)
#    stop(e)
#})
#}


## READ IN ARGUMENTS
print(opt)
df<-readRDS(opt$input)
print(opt$fixed)
variables<-c(convert_arglist(opt$fixed), convert_arglist(opt$random), convert_arglist(opt$outcome), convert_arglist(opt$hla), opt$id)
print(variables)
print(variables[!is.null(variables)])
## CHECK IF VARIABLES ARE PRESENT IN THE DF
valid<-sapply(variables[!is.null(variables)], function(var, df){var %in% colnames(df)}, df=df)
stopifnot(all(valid))

## SELECT VARIABLES OF INTEREST

df_selected<-df %>% select(one_of(variables[!is.null(variables)]))

## SELECT ONLY COMPLETE CASES

df_complete<-df_selected %>% filter(complete.cases(.))

## SAVE DF
saveRDS(df_complete, opt$output)

## CLOSE LOG CONNECTION
#if(!is.null(opt$log)){
#  sink()
#}

