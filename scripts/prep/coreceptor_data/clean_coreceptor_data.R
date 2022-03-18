#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to
#
# Code developed by Anna Hake
# Date:
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", "dplyr")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list <- list(
  make_option(c("-i", '--input'), action="store", 
	help="Input file with the FPR percentages from g2p_corecepor tool."),  
	  make_option(c("-l", "--log"), type="character", default=NULL, help="A log file.", metavar="FILE"),
  make_option(c("-o", "--output"), action="store", 
	help="Output file where the coreceptor usage is stored together with the input.")
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

getCoreceptorFromFPR<-function(fpr_perc, guideline){
  if(is.na(fpr_perc)){
	return(NA)
  }
  switch(guideline, 
	German={
	  if (fpr_perc < 5){
		return('X4')
	  } else if (fpr_perc >= 15){
		return('R5')
	  }
	  else{
		return(NA)
	  }
	},
	European={
	  if (fpr_perc < 20){
		return('X4')
	  } else {
		return('R5')
	  }
	}, 
	stop('Error: Unknown guideline provided.')
  )
}
# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# read g2p coreceptor output
g2p_output<-read.csv(opt$input, stringsAsFactors=FALSE)
# rename ID column 
g2p_output<-g2p_output %>% rename("patient_id"="header")
# create coreceptor status according to german guidelines
g2p_output<-g2p_output %>% rowwise %>% mutate(coreceptor_GER = getCoreceptorFromFPR(FPR, guideline="German"))
# create coreceptor status according to european guidelines
g2p_output<-g2p_output %>% rowwise %>% mutate(coreceptor_EU = getCoreceptorFromFPR(FPR, guideline="European"))
# extract only the coreceptor information and store it
coreceptor<-g2p_output %>% select(patient_id, coreceptor_EU, coreceptor_GER)
write.csv(coreceptor, opt$output, row.names=FALSE)


