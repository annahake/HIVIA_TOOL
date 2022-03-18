#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to convert the hla alleles in the training data set for the 
# HIVIA project. It converts the alleles to 4 digits, creates dummy variables for 
# each allele and selects only the frequent alleles.
#
# Code developed by Anna Hake
# Date: 2019-05-20 
# Update: 2019-08-27
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", "dplyr", "tidyr")

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="An rds file with the training data for a protein site", metavar="FILE"), 
  make_option(c("-c","--freq_cutoff"), type="double", default=0.01, help="Prevelance cutoff to determine frequent HLA allele.", metavar="NUM"), 
  make_option(c("-o", "--output"), type="character", default=NULL, help="An rds file with the processed binary hla alleles", metavar="FILE"), 
  make_option(c("-l", "--log"), type="character", default=NULL, help="A txt file for logging.", metavar="FILE"),  
  make_option(c("--hla"), type="character", default=NULL, help="hla genes to consider given by the config", metavar="NUM")
)
# ------------------------------------------------------------------------------
# FUNCTIONS

#' Perturbs the entries of a vector using the built-in sample function.
#'
#'
#' @param vec Vector to be perturbed
#'
#' @return Perturbed vector. 
#' 
#'
#' @export

perturb<-function(vec){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  return(vec[sample(length(vec))])
  
}

#' Checks if an HLA allele for a specific HLA gene occur frequently given a frequency cutoff and its occurrences.
#'
#'
#' @param allele numeric vector (binary := 0,1,NAs)
#' @param freq_cutoff a double
#'
#' @return A boolean denoting if this is a frequent allele. 
#' 
#'
#' @export

is_freq_allele<-function(allele, freq_cutoff){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  
  # how many samples in total 
  population_size<-length(which(!is.na(allele)))
  # determine cutoff
  cutoff<-freq_cutoff * population_size

  is_frequent<-sum(allele, na.rm=TRUE)>=cutoff
  return(is_frequent)
}


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
# ------------------------------------------------------------------------------
# MAIN

## READ ARGUMENTS

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
## open log connection
#sink(opt$log)
df<-readRDS(opt$input)
hla_names<-convert_arglist(opt$hla)

# subset the hla data (and the patient id)
df_hla<- df %>% dplyr::select(one_of(hla_names), patient_id)

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
	
frequent_alleles<-df_binary%>% 
  dplyr::select(-patient_id) %>%
  dplyr::select_if(list(~is_freq_allele(.,freq_cutoff=opt$freq_cutoff)))

# factor the new allele cols
frequent_alleles <- frequent_alleles %>% mutate_all(factor)

# append df to patient_id 
frequent_df<- dplyr::bind_cols(patient_id=df_binary$patient_id, frequent_alleles)

final_df<-df %>% dplyr::select(-one_of(hla_names)) %>% dplyr::left_join(x=.,y= frequent_df, by="patient_id")

## SAVE NEW TRAINING DATA
final_df<-tbl_df(final_df)
saveRDS(final_df, opt$output)

## LOGGING
# which alleles are not expressed
#have_N_suffix<-grepl(df_long$alleleID, pattern="N$")
#if (any(have_N_suffix)){
#  print("Following patient_ids have non-expressing alleles:")
#  print(df_long[have_N_suffix, ])
#}
#log_str<-paste0('A total of ', ncol(df_binary)-1, ' different alleles were found in the study population\n', ncol(frequent_alleles), ' of them are frequent.\n')
#print(log_str)

#sink()

