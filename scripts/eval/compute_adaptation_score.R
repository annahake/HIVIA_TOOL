#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to compute the adaptation of a virus to its host HLA profile by 
# multiplying the likelihood odds for each variant site that an amino acid 
# occurred there more likely with HLA pressure (full model) than without (null model).
# Given the predictions of the full and null model for each variant site, 
# the adaptation score is computed as described in the paper. 
#-------------------------------------------------------------------------------
required_packages<-list("optparse", 'dplyr', 'tidyr')

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c('-i', '--input'), type="character", default=NULL, help="A string with all prediction files", metavar="FILE"), 
  make_option(c("--output_pred"), type="character", default=NULL, help="The output file for the aggregated prediction file", metavar="STRING"), 
  make_option(c("--output_adapt"), type="character", default=NULL, help="The output file for the computed adaptation score", metavar="STRING"), 
  make_option(c("--outcome_type"), type="character", default=NULL, help="Either categorical or bernoulli", metavar="STRING"), 
  make_option(c("--rand"), type="character", action= "store_true", default=NULL, help="Predictions contain random samples")
)
# ------------------------------------------------------------------------------
#FUNCTIONS
# ------------------------------------------------------------------------------

#' For a specific variant site, convert the true amino acid to 'OTHER' 
#' if the amino acid is not one of the outcomes of the prediction models. 
#'
#' @param y character or character vector (random predictions) with the true amino acids at this variant site.
#' @param pred_saav character vector with the possible outcomes of the prediction model at this variant site.
#'
#' Returns the converted amino acids at the variant site. 
convert_y <-function(y, pred_saav){
  # TODO same length
  # TODO all y the same
  
  # if the true amino acid is one of the possible predictions outcome or if 'OTHER' is not a possible prediction outcome
  if (unique(y) %in% pred_saav  || !("OTHER" %in% pred_saav)){
    # return the true amino acid label
    return(y)
  } else {
    # otherwise convert the amino acid to OTHER
    return(rep("OTHER", length(y)))
  }
}


#' Computes the adaptation score.
#' 
#' 
#' @param df Data frame with the combined predictions for each site and each patient.
#' @param rand Boolean denoting if random predictions exist
#' @param outcome_type Character denoting the outcome type of the prediction model. Either bernoulli (single amino acid variant) or categorical (variant site). 
#' 
#' @return A data frame with the (transformed) adaptation score for each patient. 
 
compute_adaptation_score <-function(df, rand, outcome_type){
  #TODO: assertion for all column names
  # for categorical outcomes
  if(outcome_type == "categorical"){
    # if random predictions exist
    if(!is.null(opt$rand)){
      adapt_df<-df %>%
      # remove any existing groupings as a precaution 
      ungroup() %>% 
      # convert prediction outcome from wide format to long format. 
      # column pred contains the predicted likelihoods
      # column pred_saav contains the corresponding amino acid at this position
      # There are many possible amino acids and corresponding likelihood predictions for the same site
      gather(key="pred_saav", value = "pred", starts_with("P.Y")) %>%
      # remove dots
      mutate(pred_saav=gsub("\\.$", "", gsub("P\\.Y\\.\\.\\.", "", pred_saav))) %>%
      # substitute remaining dot for dash
      mutate(pred_saav = gsub("\\.", "-", pred_saav)) %>%
      # for each patient, variant site, random id (random and not random predictions), and model (full and null)
      group_by(patient_id, saav_id, rand_id,model) %>%
      # keep only those where there is true information about the amino acid at the variant site and there is a prediction
      filter(!is.na(y), !is.na(pred)) %>%
      # convert the true outcome if necessary to OTHER
      mutate(y = convert_y(y, pred_saav)) %>% 
      # ungroup
      ungroup %>% 
      # keep only those predictions where we have predictions for the true amino acid
      filter(y == pred_saav) %>% 
      # convert long format to wide s.t. we have new columns full and null_clin_rand
      tidyr::spread(model, pred) %>% 
      # for each patient, for each sample (random and not random), for each site
	    group_by(patient_id, rand_id, saav_id) %>% 
#	    summarize(odd = ifelse(abs(full -null_clin_rand)>0.1, full/null_clin_rand,1)) %>%
      # compute the odds between the likelihoods of the full and the null model
      summarize(odd = max(full, 0.0001)/max(null_clin_rand, 0.0001)) %>%  
      # ungroup
	    ungroup %>%
	    # group by sample -> iterate over all variant sites
	    group_by(patient_id, rand_id) %>%
	    # keep only odds which are not NA
	    # TODO: check and warn if there are NA odds
	    filter(!is.na(odd)) %>%
	    # multipliy the odds over all variant sites for a sample
	    summarize(adapt = prod(odd)) %>%
	    # transform the adaptation score
	    mutate(adapt_trans = atan(log(adapt))*2/pi) %>%
	    # sort in decreasing order
	    arrange(desc(adapt))
    } else {
      # if there are no random predictions 
      adapt_df<-df %>% 
      # ungroup
      ungroup() %>%
      # convert prediction outcome from wide format to long format. 
      # column pred contains the predicted likelihoods
      # column pred_saav contains the corresponding amino acid at this position
      # There are many possible amino acids and corresponding likelihood predictions for the same site
      tidyr::gather(key="pred_saav", value = "pred", starts_with("P.Y")) %>%
      # remove dots
      mutate(pred_saav=gsub("\\.$", "", gsub("P\\.Y\\.\\.\\.", "", pred_saav))) %>%
      # substitute remaining dot for dash
      mutate(pred_saav = gsub("\\.", "-", pred_saav)) %>%
      # for each patient, variant site, and model (full and null)
      # This is the difference to the rand== TRUE case, no iteration over random predictions
      group_by(patient_id, saav_id, model) %>%
      # keep only those where there is true information about the amino acid at the variant site and there is a prediction
      filter(!is.na(y), !is.na(pred)) %>%
      # convert the true outcome if necessary to OTHER
      mutate(y = convert_y(y, pred_saav)) %>% 
      ungroup %>% 
      # keep only those predictions where we have predictions for the true amino acid
      filter(y == pred_saav) %>% 
      # convert long format to wide s.t. we have new columns full and null_clin_rand
      tidyr::spread(model, pred) %>% 
      # for each patient and for each site
	    group_by(patient_id, saav_id) %>% 
#	    summarize(odd = ifelse(abs(full -null_clin_rand)>0.1, full/null_clin_rand,1)) %>%
      # compute the odds between the likelihoods of the full and the null model
      summarize(odd = max(full, 0.0001)/max(null_clin_rand, 0.0001)) %>% 
	    ungroup %>%
	    # for each patient -> iterate over all sites
	    group_by(patient_id) %>%
	    # keep only odds which are not NA
	    # TODO: check and warn if there are NA odds
	    filter(!is.na(odd)) %>%
	    # multipliy the odds over all variant sites for a sample
	    summarize(adapt = prod(odd)) %>%
	    # transform the adaptation score
	    mutate(adapt_trans = atan(log(adapt))*2/pi) %>%
	    # sort in decreasing order
	    arrange(desc(adapt))
    }
  } else if (outcome_type=="bernoulli"){
    if(!is.null(rand)){
        adapt_df<-df %>% 
        # ungroup
        ungroup() %>% 
        # select only those predictions where the true amino acid was predicted
	      filter(y == 1) %>%
	      # select only relevant columns
	      select(patient_id, model, Estimate, rand_id, saav_id) %>% 
	      # convert long format to wide format by 
	      # having a column for the estimates of the full and the null model
        tidyr::spread(model, Estimate) %>%
        # for each patient, and each single amino acid variant (saav)  
	      group_by(patient_id, rand_id, saav_id) %>% 
	      # compute odds
	      summarize(odd = max(full, 0.0001)/max(null_clin_rand, 0.0001)) %>% 
	      # ungroup
	      ungroup %>%
	      # for each patient (random and non-random samples), iterate over all SAAVs
	      group_by(patient_id, rand_id) %>%
	      # keep only the non NA odds
	      # TODO keep track of NA odds and raise warnings
	      filter(!is.na(odd)) %>%
	      # compute adaptation score by multiplying the odds
	      summarize(adapt = prod(odd)) %>%
	      # compute transformed adaptation score
	      mutate(adapt_trans = atan(log(adapt))*2/pi) %>%
	      # sort in decreasing order
	      arrange(desc(adapt))
    } else {
        adapt_df<-df %>% 
        # ungroup
        ungroup() %>% 
        # select only those predictions where the true amino acid was predicted
	      filter(y == 1) %>%
	      # select only relevant columns
	      select(patient_id, model, Estimate, saav_id) %>% 
	      # convert long format to wide format by 
	      # having a column for the estimates of the full and the null model
        tidyr::spread(model, Estimate) %>%
        # for each patient, and each single amino acid variant (saav)  
	      group_by(patient_id, saav_id) %>% 
	      # compute odds
	      summarize(odd = max(full, 0.0001)/max(null_clin_rand, 0.0001)) %>% 
	      # ungroup
	      ungroup %>%
	      # for each patient, iterate over all SAAVs
	      group_by(patient_id) %>%
	      # keep only the non NA odds
	      # TODO keep track of NA odds and raise warnings
	      filter(!is.na(odd)) %>%
	      # compute adaptation score by multiplying the odds
	      summarize(adapt = prod(odd)) %>%
	      # compute transformed adaptation score
	      mutate(adapt_trans = atan(log(adapt))*2/pi) %>%
	      # sort in decreasing order
	      arrange(desc(adapt))
    }
  
  } else {
    stop("Error unknown outcome type provided. See -h for help")
  
  
  }

  return(adapt_df)





}
# -------
# MAIN
# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# initialize prediction data frame
pred_df<-NULL
# open connection to file with all the separate prediction file names
con = file(opt$input, "r")
while ( TRUE ) {
  # one line is one file name
  fname = readLines(con, n = 1)
  # if there are only empty lines left/no line 
  if ( length(fname) == 0 ) {
    # break from while loop
    break
  }
  # read data frame from file name
  df <-read.csv(fname, stringsAsFactors=FALSE)
  # bind current data frame to existing data
  pred_df<-bind_rows(pred_df, df)
}
# close the connection
close(con)
# save the prediction per variant sites
saveRDS(pred_df, opt$output_pred)


#TODO check that column names exist (model, pred, patient_id, saav_id), what about protein...per protein or together?
if("saav" %in% colnames(pred_df)){
  pred_df <- pred_df %>% rename(saav_id = saav)
}
switch(opt$outcome_type, 
  bernoulli = {adapt_score = compute_adaptation_score(pred_df, outcome_type = 'bernoulli', rand = opt$rand)}, 
  categorical = {adapt_score = compute_adaptation_score(pred_df, outcome_type = 'categorical', rand = opt$rand)}, 
  stop('Provided outcome type not supported')
)



write.csv(adapt_score, opt$output_adapt, row.names = FALSE)
