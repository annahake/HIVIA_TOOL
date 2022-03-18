#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to convert the hla alleles to binary matrix in the test set.
#
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", 'brms', 'dplyr', 'tidyr')

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", '--input'), type="character", default=NULL, help="An rds file containing the test data set", metavar="FILE"), 
  make_option(c("--hla"), type="character", default=NULL, help="A list of hlas", metavar="STRING"),
  make_option(c("-o", '--output'), type="character", default=NULL, help="Test set with binarized hla alleles", metavar="FILE")
  # TODO: Add comma if more than one option is used. 
)
# ------------------------------------------------------------------------------
# FUNCTIONS
# TODO same functions as in convert_train_hla2binary.R. Consider having a single function libraries for functions that are often used (wrt updates/debugging)


#' Checks if at least one entry in binary vector is set to 1.
#'
#'
#' @param binary a numeric vector containing only 1,0s, and NAs. 
#
#' @return Returns 1, if at least one entry in binary vector is set to 1, else 0.
#' 
#'
#' @export

any.binary<-function(binary){
  # assertive tests for each argument
  ## binary is only allowed to contain 0,1, and NA
  stopifnot(all(binary %in% c(0,1,NA)))
  # missing values/coercion handling
  
  # warnings
  
  # code
  # TODO: alternative is that sum >= 1 (na.rm = TRUE)
  return(as.numeric(any(as.logical(binary))))
}


#' Converts list arguments ("X,X,X") to vector of characters.
#'
#'
#' @param char character representing the list argument
#'
#' @return vector of characters where each element represents element of the 
# argument list
#' 
#'
#' @export

convert_arglist <-function(char){
  # check that character is not null or na
  stopifnot(!is.null(char))
  stopifnot(!is.na(char))
  # ensure that input is a character
  stopifnot(is.character(char))
  # split the string into single characters
  char_vec<-unlist(strsplit(char, split=""))
  #TODO trim trailing spaces
  # TODO use grepl on char instead
  # if there is a comma split according to the comma
  if ("," %in% char_vec){
    args<-unlist(strsplit(char, split="\\,"))
  }
  else {
    args <-char
  }
  return(args)

}

#------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
df<-readRDS(opt$input)
hla_names<-convert_arglist(opt$hla)

select_vars<-c(hla_names, 'patient_id')
# ensure that all column names are in data frame
stopifnot(all(sapply(select_vars, function(var,df){var %in% colnames(df)}, df=df)))
# select variables
df_hla<- df %>% dplyr::select(one_of(select_vars))

# convert to long format
df_long <- df_hla %>% tidyr::gather(key="HLA", "alleleID", -patient_id)
## BINARIZE IDs
# create dummy variables for alleles
df_binary_long<-df_long %>% 
	# create dummy variables
	model.matrix(~ alleleID -1, data=.)  %>% 
	# convert matrix to data.frame
	as.data.frame() %>%
	# combine with original df
	bind_cols(df_long%>% filter(complete.cases(.)), .)


df_binary<-df_binary_long %>% 
	# drop HLA and alleleID
	select(-HLA, -alleleID) %>% 
	# for each patient_id currently there is a binary vector for all 6 different HLA genes (2 alleles each)
	# Note, due to the dummy variable (wide format) representation there is no need anymore to have the 12 entry vector per patient
	group_by(patient_id) %>% 
	# each distinct HLA allele can only appear 0,1, or 2 times per patient and only in one of the 6 genes
	# consequently, the binary vector per patient can be shrunk
	# as haplotype information is not considered separately, the 12 entry patient vector can be combined to a single value: 
	# presence or abscence of the allele: 
	summarize_at(vars(-group_cols(), starts_with("alleleID")),list(any.binary)) %>% 
	# convert to factor
	#mutate_at(vars(starts_with('alleleID')), list(factor)) %>%
	# drop the prefix
	rename_at(vars(starts_with('alleleID')), list(~gsub('alleleID', '', .)))
	
df_binary <- df_binary %>% mutate_at(vars(-patient_id), list(factor))

# deselect all hlas which are not of interest
final_df<-df %>% dplyr::select(-one_of(hla_names)) %>% dplyr::left_join(x=.,y= df_binary, by="patient_id")
final_df<-tbl_df(final_df)
saveRDS(final_df, opt$output)




