#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to predict the amino acid at a specific protein site given 
# the HLA allele information and the fitted Bayesian GLMM.
# In addition the visual and theoretical MCMC diagnostics, as well as the model 
# are stored 
#
#-------------------------------------------------------------------------------
# LIBRARIES./
required_packages<-list("optparse", 'brms', 'dplyr')

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-f", "--fit"), type="character", default=NULL, help="RDS file with the brms fit", metavar="FILE"), 
  make_option(c("-t", "--test"), type="character", default=NULL, help="RDS file with the test data", metavar="FILE"),
  make_option(c("-o", '--output'), type="character", default=NULL, help="CSV file with the prediction", metavar="FILE"), 
  make_option(c("-r", '--test_rand'), type = "character", default=NULL, help="RDS file with the randomized test data", metavar="FILE"), 
  make_option(c("--model"), type="character", default=NULL, help="Name of the model (full or null)", metavar="STRING"), 
  make_option(c("--model_version"), type="character", default=NULL, help="Name of the model-version (hlaboth_clin, hla1, or hla2)", metavar="STRING"), 
  make_option(c("-p", "--protein"), type="character", default=NULL, help="Name of the HIV protein used (gag_p24 or env)", metavar="STRING"),
  make_option(c("-s", "--site"), type="character", default=NULL, help='Site in the protein', metavar="STRING"), 
  make_option(c("--seed"), type="character", default=NULL, help='random seed', metavar="STRING")
  # TODO: Add comma if more than one option is used. 
)
# ------------------------------------------------------------------------------
# FUNCTIONS


#' Format test data. 
#'
#' @description Checks if all required features are present. If missing, HLA allele features are set to zero. Checks if levels of categorical features exist in model fit. Samples with new levels can't be predicted and are removed. 
#'
#' @param df_ Data frame with the test data
#' @param fit Brms fit
#'
#' @return Data frame with corrected format. 
#' 
#'
#' @export


fit_format <- function(df_, fit){

	### get colnames of the training data
	fit_vars<-colnames(fit$data)
	### remove the outcome variable, only feature variable are of interest
	y_ind<-which(fit_vars=='y')
	fit_vars<-fit_vars[-y_ind]
	# get colnames of the test data
	test_vars<-colnames(df_)

	### check which features in the test set are not appearing in the model
	#omit_vars<-setdiff(test_vars, fit_vars)
	### delete those
	#df_ <- df_ %>% select(-one_of(omit_vars))

	### check which features are missing in the test set
	add_vars<-setdiff(fit_vars, test_vars)
	#print(add_vars)
	# if the variables are hla features create zero vectors (not present in the test data)
	if(length(add_vars)>0){
		for(var in add_vars){
			if(grepl("\\w\\_\\d{2,3}\\_\\d{2,3}", var)){
		  	#print(var)
				df_<-df_ %>% mutate(!!var := factor(rep(0,nrow(df_)), levels=c('0', '1')))
			} else {
				#print(var)
				warning(paste('variable:', var, 'present in fit but not in test set.'))
			}
		}
	}

	## FROM HERE ON FIT AND DF_TEST HAVE SAME FEATURES/COLNAMES
	#stopifnot(length(setdiff(fit_vars, colnames(df_)))==0)


	## remove any samples that have levels which are not present in the fit
	## doesn't work if only one sample

	## for each categorical feature, 
	fac_var<-fit$data %>% select_if(is.factor) %>% colnames
	for (var in fac_var){
		# which is not the outcome variable, the patient id, or the random id
  	if(!(var %in% c('y', "patient_id", "rand_id") )){
  		# get the levels of the variable in the fit
  		fit_levels<-levels(fit$data[, var])
  		# get the levels of the variables in the test data
  		test_levels<-levels(as.factor(df_[[var]]))
  		# get the levels that are in the test data but not in the fit data
  		omit_levels<-setdiff(test_levels, fit_levels)
  		# create data set where the categorical variable has same levels as in fit
  		df_ <- df_ %>% 
  			filter(!(get(var) %in% omit_levels)) %>%
				mutate(!!var:=factor(get(var)))
  	}
	}
	# factor all categorical variables in test data
	df_<- df_ %>% mutate_if(sapply(df_, is.factor), factor)
	# return test data
	return(df_)
}


#' Scale numerical data according to fit data.
#'
#' @description Scale numerical features as the age variable with the same parameters as in the fit.
#'
#' @param df_ Data frame with the test data
#' @param fit Brms fit
#'
#' @return Test data with normalized numerical features with same center and scale as in the fit data.
#' 
#'
#' @export

scale_<-function(df, fit){
  num_col<-sapply(fit$data, is.numeric)
  if(any(num_col)){
    num_names<-colnames(fit$data)[num_col]
    for (n in num_names){
      df[, n] <- scale(df[, n], center = attr(fit$data[,n], "scaled:center"), scale = attr(fit$data[,n],"scaled:scale"))
    }
  }

  return(df)
}

# ------------------------------------------------------------------------------
# MAIN
# silence warnings if desired. 
# Set warnings back to old setting at the end of the script.
oldw <- getOption("warn")
options(warn = -1)

# READ IN PARAMETERS

## set option parser
opt_parser<-OptionParser(option_list=option_list)
## read arguments
opt<-parse_args(opt_parser)
## read in fit
fit<-readRDS(opt$fit)
## read in test
df_test<-readRDS(opt$test)
## set seed (This is important as otherwise the predictions might differ slightly)
set.seed(opt$seed)
# PREDICT TEST DATA
## format test data frame to fit the model train data
df_test_formatted <- fit_format(df_test,fit)
df_test_scaled<-scale_(df_test_formatted, fit)
#if ('y' %in% colnames(df_test_scaled)){
#	df_test_scaled<-df_test_scaled %>% select(-y)
#}

# if random test data is provided
if(!is.null(opt$test_rand)){
	# perform the same steps as with the normal test data
	## read in randomized test data
	df_test_rand <- readRDS(opt$test_rand)
	df_test_rand_formatted <- fit_format(df_test_rand, fit)
	df_test_rand_scaled<- scale_(df_test_rand_formatted, fit)
}

#TODO: check for empty data frames 
if ('y' %in% colnames(df_test_scaled)){
	preds<-predict(fit, df_test_scaled %>% select(-y), allow_new_levels = TRUE)
	preds_df<-data.frame(protein = opt$protein, model = opt$model, model_version = opt$model_version, patient_id= df_test_scaled$patient_id, preds, saav_id = opt$site,y = df_test_scaled$y)
} else {
	preds<-predict(fit, df_test_scaled, allow_new_levels = TRUE)
	preds_df<-data.frame(protein = opt$protein, model = opt$model, model_version = opt$model_version, patient_id= df_test_scaled$patient_id, preds, saav_id = opt$site)
}

# prepend the true label and the patient id
#preds_df<-data.frame(protein = opt$protein, model = opt$model, model_version = opt$model_version, patient_id= df_test_scaled$patient_id, preds, saav_id = opt$site,y = df_test_scaled$y)
#if ('y' %in% colnames(df_test_scaled) && 'saav_id' %in% colnames(df_test_scaled)){
#	preds_df<-data.frame(preds_df, saav_id = unique(df_test_scaled$saav_id),y = df_test_scaled$y )
#}

if(!is.null(opt$site)){
	preds_df<-data.frame(preds_df, site = opt$site)
}

# if test data is provided
if(!is.null(opt$test_rand)){
	# add the outcome variable
	df_test_rand_full<- left_join(df_test_rand_scaled, df_test %>% select(patient_id, y), by = "patient_id")
	# get predictions
	rand_preds<-predict(fit, df_test_rand_scaled, allow_new_levels=TRUE)
	# add further variables of interest
	rand_preds_df <-data.frame(protein = opt$protein, model = opt$model, model_version = opt$model_version, saav_id = unique(df_test$saav_id), patient_id= df_test_rand_scaled$patient_id, y = df_test_rand_full$y, rand_preds,rand_id = df_test_rand_scaled$rand_id)
print('after rand')
	# combine random predictions with true predictions
	preds_df <-bind_rows(preds_df,rand_preds_df)
}



if(!dir.exists(dirname(opt$output))){
	dir.create(dirname(opt$output))
}

# save result
write.csv(preds_df, opt$output, row.names=FALSE)

# STORE MODEL CHARACTERISTICS
#if(all(rhat(fit)<0.99)){
#	warning('rhat smaller than 0.99')
#}
#model_info=list(summary=summary(fit),logp=log_posterior(fit), nuts_params=nuts_params(fit), rhat=rhat(fit), neff_ratio=neff_ratio(fit))
## save model summary
#saveRDS(summary(fit), opt$info_output)

# STORE GGMCMC PLOT/too much storage space, not informative
#S<-ggs(fit)
#ggmcmc(S, file=opt$mcmc_plot)
# check MCMC diagnostic
traceback()
options(warn = oldw)
