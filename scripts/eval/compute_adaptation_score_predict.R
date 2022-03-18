#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
required_packages<-list("optparse", 'dplyr', 'tidyr', 'seqinr')

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c('-i', '--input'), type="character", default=NULL, help="A csv file with all predictions", metavar="FILE"), 
  make_option(c("-s", '--seq'), type="character", default=NULL, help='CSV file with the sequence information', metavar="FILE"), 
  make_option(c('-o', "--output"), type="character", default=NULL, help="The output file for the computed adaptation score", metavar="FILE"), 
  make_option(c('-r', '--ref_map'), type = 'character', default = NULL, help = 'Mapping file from aligned Consenus C reference sequence and the HXB2 reference sequence', metavar = 'FILE'),
  make_option(c("--outcome_type"), type="character", default=NULL, help="Either categorical or bernoulli", metavar="STRING")
)
# ------------------------------------------------------------------------------
#FUNCTIONS
# ------------------------------------------------------------------------------

get_ref_pos<-function(saav_ids, ref_map, patient_ids){
	# get the mapping index of current saav_id
	# for each saav_id of each patient get the reference index
	# TODO rename from saav_id to seq_position
	ref_pos_vec<-sapply(seq(saav_ids), function(row_ind, ref_map){
		# get current saav id of current patient
		saav_id<-saav_ids[row_ind]
		# get the corresponding row in the map file
		ind<-which(as.character(ref_map[,3]) == as.character(saav_id))
		# if there is no mapping at all, this is an insertion at the end wrt reference sequence 
		# TODO this is an actual error in the mapping (* is aligned to gap, but discarded for the reference mapping)
		if (length(ind) == 0){
			print(saav_id)
			print(patient_ids[row_ind])
			count<-0
			while(length(ind) == 0) {
				count = count + 1
				# check if there is a mapping for the precedent position
				saav_id<-as.numeric(saav_id) -1
				ind<-which(as.character(ref_map[,3]) == as.character(saav_id))
				print(count)
			}
			ref_pos <-paste0(ref_map[ind, 2], '_', count)
			print(ref_pos)
		} else {
			# TODO in case that position is an insertion wrt reference, the map will be a gap
			ref_pos<-ref_map[ind, 2]
		}
		return(ref_pos)
	}, ref_map = ref_map)
	return(ref_pos_vec)
	

}



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
# read in predictions
pred_df<-read.csv(opt$input)
pred_df <- pred_df %>% mutate(patient_id = as.character(patient_id))
# read in sequences
seqs_df <- read.fasta(opt$seq, as.string = FALSE)
# read in map file
ref_map<-read.csv(opt$ref_map, stringsAsFactors=FALSE, skip=1)
# trim spaces #TODO where are the spaces, what is stored at column three?
ref_map[,3]<-gsub(" ", "", ref_map[,3])

# format sequences in long format
seqs_df<-do.call(rbind, seqs_df)
seqs_df<- as_tibble(data.frame(patient_id = rownames(seqs_df), seqs_df))

seqs_df_long <- seqs_df %>%
	gather(saav_id, y, -patient_id) %>% 
	mutate(saav_id = gsub('X','', saav_id)) %>% 
	mutate(saav_id = get_ref_pos(saav_id, ref_map, patient_id)) %>%
	# patient id from sequence file has to match patient id from hla file
	mutate(patient_id = gsub('_Capsid', '', patient_id))

pred_df$saav_id<-as.character(pred_df$saav_id)
comb_df <- left_join(pred_df, seqs_df_long, by = c('patient_id', 'saav_id'))

# TODO currently only categorical is supported, for bernoulli combining sequence and predictions have to be expanded
switch(opt$outcome_type, 
  #bernoulli = {adapt_score = compute_adaptation_score(comb_df, outcome_type = 'bernoulli', rand = NULL)}, 
  categorical = {adapt_score = compute_adaptation_score(comb_df, outcome_type = 'categorical', rand = NULL)}, 
  stop('Provided outcome type not supported')
)


write.csv(adapt_score, opt$output, row.names = FALSE)
