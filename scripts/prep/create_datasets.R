#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to
#
# Code developed by Anna Hake
# Date:
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", 'dplyr', 'lubridate')

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c('--hla'), type="character", default=NULL, help="A csv file with the hla allele occurrence.", metavar="FILE"), 
  make_option(c("-s", "--saav"), type="character", default=NULL, help="A csv file with binary columns denoting the single amino acid variants", metavar="FILE"),
  make_option(c('--coreceptor'), type="character", default=NULL, help="A csv file with the predicted coreceptor usage from the g2p coreceptor tool", metavar="FILE"), 
  make_option(c('--clinical'), type="character", default=NULL, help="A csv file with the clinical data", metavar="FILE"), 
  #make_option(c('-o', '--all_output'), type="character", default=NULL, help="Filename of the combined dataset", metavar="FILE"),
  make_option(c('-d', '--output_dir'), type="character", default=NULL, help="Directory where the datasets per SAAV should be stored.", metavar="DIR"),
  #make_option(c('--log'), type="character", default=NULL, help="A file for logging", metavar="FILE"), 
  make_option(c("--cd4_cutoff"), default=NULL, help="CD4 cutoff to distinguish chronics from early patients", metavar="INT")#, 
  #make_option(c("--rand"), default=NULL, help="Either the HLA features (alleles) are randomized ('rand') or not ('norand').", metavar="STRING"),
  #make_option(c("--seed"), default=NULL, help="Seed for the randomization", metavar="INT")
  # TODO: Add comma if more than one option is used. 
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
#sink(opt$log)
#df_names<-c('coreceptor', 'saav', 'clinical','hla')
#df_list<-lapply(df_names, function(name){name<-read.csv(opt[[name]], stringsAsFactors=FALSE)})
coreceptor_df<-read.csv(opt$coreceptor, stringsAsFactors=FALSE)
hla_df<-read.csv(opt$hla, stringsAsFactors=FALSE)
clinical_df<-read.csv(opt$clinical, stringsAsFactors=FALSE)
saav_df<-read.csv(opt$saav, stringsAsFactors=FALSE)

## PERTURB HLA FEATURES IF ASKED
#if(opt$rand=="rand"){
#  set.seed(opt$seed)
#  hla_df <-hla_df %>% mutate_at(.vars=vars(-patient_id), .funs=perturb)
#}
#join only those that are in sequences (correct subtype, 1 stop codon in the end, 1 frame shift in the beginning)
df_selected <- left_join(saav_df, coreceptor_df, by='patient_id') %>% left_join(hla_df, by= 'patient_id') %>% left_join(clinical_df, by='patient_id')
#
# filter those that are treatment naiive: Currently report which are not treatment naiive but not discarded!
#print('Study participants that are not marked as offART:')
#df_selected %>% dplyr::filter(saav_id==df_selected$saav_id[1], onART!='NO' | is.na(onART)) %>% select(patient_id,cd4_count,  vl, onART)
df_selected <- df_selected %>% dplyr::filter(onART == "NO", cd4_count <= as.numeric(opt$cd4_cutoff))
## FORMATTING
# format all dates
df_selected <- df_selected %>% mutate_at(vars(ends_with('date')), as_datetime)

# format all integer
df_selected<-df_selected %>% mutate_at(vars(age, years_of_infection), as.integer)

# format all double
df_selected<-df_selected %>% mutate_at(vars(log_vl, cd4_count), as.double)
# format all factors
df_selected<-df_selected %>% mutate_at(vars(-vl, -age, -years_of_infection, -log_vl, -cd4_count, -ends_with('date')),factor)
#TODO: delete all NAs in y
df_selected<-df_selected %>% filter(!is.na(y))
df_selected<-tbl_df(df_selected)
# save combined data for all saav
#saveRDS(df_selected, opt$all_output)
#print("after saving combined output")

# save df per SAAV
df_list<-df_selected %>% group_split(saav_id)
void<-lapply(df_list, function(df){
	#check if the saav is still frequent,  
	# TODO general selection should be before creating data set and selecting SAAVs
	# frequent cutoff should be passed not hardcoded
	if (sum(as.numeric(df$y)) >= 0.01 * length(df)){
		saav_id<-unique(df$saav_id)
  	#file.path for all options should work as well
  	output_dir<-file.path(opt$output_dir,saav_id, "full_data")
  	if (!dir.exists(output_dir)){
  	  dir.create(output_dir, recursive=TRUE)
  	}
  	saveRDS(df, file.path(output_dir, 'df.rds'))
	}

  })
#sink()


