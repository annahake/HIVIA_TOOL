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
required_packages<-list("optparse", "dplyr",'lubridate')

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c('--hla'), type="character", default=NULL, help="A csv file with the hla allele occurrence.", metavar="FILE"),# TODO: Add comma if more than one option is used.,
  make_option(c("-c", '--clinical'), type="character", default=NULL, help="A csv file with the clinical data", metavar="FILE"), 
  make_option(c('-o', '--output'), type="character", default=NULL, help="File where the randomized dataset is stored.", metavar="FILE"),
    make_option(c("--cd4_cutoff"), default=NULL, help="CD4 cutoff to distinguish chronics from early patients", metavar="INT"), 
  make_option(c("--rand_nr"), default=NULL, help="The number of random perturbations to perform", metavar="INT"),
  make_option(c("-s", "--seed"), default=NULL, help="Seed for the randomization", metavar="INT")
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
perturb<-function(vec){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  return(vec[sample(length(vec))])
  
}

# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
rand_nr<-as.numeric(opt$rand_nr)
hla_df<-read.csv(opt$hla, stringsAsFactors=FALSE)
hla_names<-hla_df %>% select(-patient_id) %>% colnames %>% unlist
clinical_df<-read.csv(opt$clinical, stringsAsFactors=FALSE)

# TODO: selection of samples should be in separate script?



#not used anymore since perturbation has to be on final df
#hla_df <-hla_df %>% mutate_at(.vars=vars(-patient_id), .funs=perturb)

df_selected <- left_join(hla_df, clinical_df, by= 'patient_id') 
df_selected <- df_selected %>% dplyr::filter(onART == "NO", cd4_count <= as.numeric(opt$cd4_cutoff))
df_selected <- df_selected %>% mutate_at(vars(ends_with('date')), as_datetime)

# format all integer
df_selected<-df_selected %>% mutate_at(vars(age, years_of_infection), as.integer)

# format all double
df_selected<-df_selected %>% mutate_at(vars(log_vl, cd4_count), as.double)
# format all factors
df_selected<-df_selected %>% mutate_at(vars(-vl, -age, -years_of_infection, -log_vl, -cd4_count, -ends_with('date')),factor)
df_selected<-tbl_df(df_selected)

# create perturbed data sets
set.seed(opt$seed)
perturb_dfs_list<-list()
# make number of random perturbation a parameter
for (i in 1:rand_nr){
	#add run_id?
	perturb_dfs_list[[i]] <-df_selected %>% mutate_at(.vars = vars(one_of(hla_names)), .funs=perturb)
}

perturb_df <- bind_rows(perturb_dfs_list, .id='rand_id')



saveRDS(perturb_df, opt$output)
