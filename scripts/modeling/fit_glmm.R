#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 
#-------------------------------------------------------------------------------
# This is R code to fit a bernoulli or categorical GLMM model (Bayesian and 
# non-Bayesian) to predict the HLA adaptation of a particular site in a protein 
# of HIV-1 while incorporating the phylogenetic relatedness of the viral samples 
# as group-level effect. 
#
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", "brms", "ape", "MCMCglmm", "dplyr")

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-i", "--data"), type="character", default=NULL, help="RDS object where the data per protein site is stored", metavar="FILE"),
  make_option(c("-t", "--tree"), type="character", default=NULL, help="Tree object representing the viral phylogeny of the samples", metavar="FILE"), 
  make_option(c('--model'), type = 'character', default=NULL, help = 'Determines if the null (wihtout hla information) or the full model is built. One of the following: "null_clin_rand" or "full".'), 
  make_option(c("--method"), type='character', default=NULL, help ="Determines the glmm package to use. Currently only the brms package is supported."), 
  make_option(c("-p", "--prior"), type='character', default=NULL, help ="Determines the prior for the coefficients."), 
  make_option(c("--adapt_delta"),type='numeric', default=NULL, help="Parameter to control the MCMC sampling."), 
  make_option(c("-s", "--seed"), type = "integer", default=100, help="seed", metavar="INT"), 
  make_option(c("-o", "--output"), type = "character", default=NULL, help="Model output filename", metavar="FILE"), 
  #make_option(c("-l", "--log"), type = "character", default=NULL, help="Log filename", metavar="FILE"),
  make_option(c("-f", "--fixed"), type = "character", default=NULL, help="fixed variables to include", metavar="STRING"), 
  make_option(c("-r", "--random"), type = "character", default=NULL, help="random variables to include", metavar="STRING"), 
  make_option(c("-y", "--outcome"), type = "character", default=NULL, help="outcome variable", metavar="STRING"), 
  make_option(c("--coeffs"), type = "character", default=NULL, help="File name where the coefficients of the model are stored", metavar="STRING"), 
  make_option(c('--outcome_type'),type = "character", default=NULL, help = 'Either categorical or bernoulli', metavar = 'STRING' ), 
  make_option(c("--fitted"),type = "character", default=NULL, help = 'File name where the predictions for the training samples are stored.', metavar = 'STRING') ,
  make_option(c("--saav"), type = "character", default=NULL, help = 'Protein site (categorical fit) or singla amino acid variant position (bernoulli fit)', metavar = 'STRING'), 
  make_option(c("--protein"), type = "character", default=NULL, help = 'HIV protein', metavar = 'STRING'), 
  make_option(c("--model_version"), type = "character", default=NULL, help = 'Defining which HLA features are taken into account (hlaboth_clin, hla1, hla2, hla1_clin, hlaboth_clin, ...)', metavar = 'STRING'), 
  make_option(c('--ref_map'), type = "character", default=NULL, help = 'Table with the HXB2 reference amino acids per position', metavar = 'STRING'), 
  make_option(c('--par_ratio'), type = "character", default=NULL, help = 'Ratio of the expected number of non-zero coefficients to the expected number of zero coefficients in the horseshoe prior', metavar = 'STRING')
)
# ------------------------------------------------------------------------------
# FUNCTIONS

#' Calls and fits a glmm model.
#'
#'
#' @param df Data frame with the input data.
#' @param phylo Tree object with the viral phylogeny of the samples in the data frame.
#' @param model Character denoting if the null model (wihtout HLA information) or the full model is built. 
#' @param method Character denoting the glmm method/R package to use.
#' @param fixed Character vector with column names of the data frame that should enter the model as fixed variable.
#' @param random Character with the column name of the random component. Currently only one random-intercept model is possible. 
#' @param outcome Character with the column name of the outcome variable. 
#' @param prior Character denoting the prior for the coefficients of the Bayesian GLMM (method = brms). 
#' @param seed Integer for initializing the random pseudogenerator. 
#' @param adapt_delta Parameter for the MCMC sampling for the brms method. 
#' @param sample_prior Indicate if samples from priors should be drawn additionally to the posterior sample. Set to no. 
#' @param outcome_type Character denoting if categorical or bernoulli fit should be performed.
#' @param par_ratio Parameter for the horseshoe prior for the brms method. 
#'
#' @return The fitted model. 
#'
#' @export

build_glmm<-function(df, phylo, model, method, fixed, random, outcome, prior=NULL, seed, adapt_delta, sample_prior, outcome_type, par_ratio = NULL){
  # assertive tests for each argument
  valid_models<-c("null_clin_rand",  "full")
  if (!model %in% valid_models){
    stop(paste0("Unsupported 'model' type provided. Model must be one of the following: \n", paste(valid_models, collapse="\n")))
  }
  valid_methods <- c("MCMCglmm", "brms")
  if (!method %in% valid_methods){
    stop(paste0("Unsupported 'method' option provided. Method must be one of the following: \n", paste(valid_methods, collapse="\n")))
  }
  valid_priors = c(NULL, "weak1","weak2", "horseshoe", "NULL")
  if (!prior %in% valid_priors | is.null(prior)){
    stop(paste0("Unsupported 'prior' option provided (", prior, ")", ". Priors have to be one of the following:\n", paste(valid_priors, collapse="\n")))
  }  
  # TODO: 
  # assert that dimension of tree/resulting covariance matrix is the same as random component. 
	# Assert that ids of tree match random component

  # missing values/coercion handling
  # warnings
  
  # code

  ## get phylogenetic information 
  phylo_info <- get_phylogeny(tree=phylo, method=method)
  ## set priors if option prior is set
  if (!is.null(fixed) && !is.null(prior) && prior != "NULL"){
		priors <- setup_priors(prior_type=prior, par_ratio= par_ratio, method=method)
  } else {
    priors <- NULL
  }
  ## build formula
  formula <-build_formula(method=method, model=model, fixed=fixed, random =random, outcome =outcome)
  ## call glmm method
  switch(method, 
    MCMCglmm={
      set.seed(seed)
      rand_formula<-as.formula(paste0('~', random))
      # TODO change to random variable to random input
      mod<-MCMCglmm(fixed=formula, random=rand_formula, data=df, family="categorical", prior=priors, ginverse=phylo_info, verbose=FALSE)
    },
    brms={
      # note that brm formula consists of both fixed and random part
      # TODO change cov_ranef input to random variable
      #print(df[,outcome])
      set.seed(seed)
      mod<-brm(formula=formula, cov_ranef=list(patient_id = phylo_info), data=df, family=outcome_type, prior=priors, seed=seed, control = list(adapt_delta = adapt_delta), sample_prior=sample_prior, verbose=TRUE)
    }, 
    {stop("unsupported 'method' option provided")}
  )
  return(mod)
}




#' Setup formula.
#'
#' @description Set up of model formula based on the chosen method and model, given the fixed, random and outcome variable. 

#' @param method The package method name. Either brms or MCMCglmm. Note, that currently only brms is used.
#' @param model The type of model. Either null_clin_rand (not HLA pressure) or full (with HLA pressure). 
#' @param fixed A vector with the fixed predictors in the model. 
#' @param random The name of the group-level effect variable. Currently, only one group-level effect is modelled as intercept only model. 
#' @param outcome The name of the outcome variable. 
#'
#' @return A formula object. 
#' 
#' @examples
#'
#' @export
build_formula<-function(method, model, fixed, random, outcome){
  # assertive tests for each argument
  # if model = single, then hla has to be not NULL
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  # select subset df

# ensure all factor variables are factors 
# build formula

#	print(paste0('method: ', method))
#	print(paste0('model: ', model))
#	print(paste0('fixed: ', fixed))
#	print(paste0('random: ',random))
#	print(paste0('outcome: ',outcome))
# TODO: switch not necessary since only brms method is used and tested
  switch(method, 
#    MCMCglmm = {
#      switch(model, 
#        null = {
#          formula= as.formula(paste0(outcome,  " ~. "))
#        }, 
#        {
#          formula=as.formula(paste0(outcome, " ~ ", paste(fixed, collapse="+")))
#        }
##        contr_rand_single = {
##          stopifnot(!is.null(hla), "Single HLA allele not provided. Model option single requires a single hla allele.")
##          formula=as.formula(paste("y ~ age", "coarse_race", "sex", "years_of_infection" , hla , sep="+"))
##        }, 
##        full = {
##          hla_column_names<-df %>% select(matches("\\w\\.\\d")) %>% colnames
##          formula=as.formula(paste("y ~ age", "coarse_race", "sex", "years_of_infection" ,  paste(hla_column_names, collapse="+"), sep="+"))
##        }, 
##        {stop("unknown model provided in build_formula")}
#      
#      
#      )
#    }, 
    brms = {
    	#TODO no switch necessary anymore since same formula structure. 
      switch(model, 
        null_clin_rand={
        	# TODO: check if this can be the case
        	if(is.null(fixed)){
        		formula=as.formula(paste0(outcome,  " ~ ", "(1|", random, ")"))
        	} else {
        		formula=as.formula(paste0(outcome, " ~ ", paste(fixed, collapse="+"), "+(1|", random, ")"))
        	}
        },
        full={
          formula=as.formula(paste0(outcome, " ~ ", paste(fixed, collapse="+"), "+(1|", random, ")"))
        }, 
        {stop("unknown model provided in build_formula")}
      )
    }, 
    {stop("unknown method provided in build_formula")}
  )
	#print(formula)
  return(formula)
}

#' Set up priors.
#'
#' @description Set up the prior for the coefficients of the model. Currently, the same prior is set on all coefficients. In the brms package this can be extended to set prior for different population-level effects as well as group-level effects. 
#'
#' @param prior_type name of the prior
#' @param method package used for the modelling, either MCMCglmm or brms
#' @param par_ratio additional parameter for the horseshoe prior in the brms package
#'
#' @return A prior object with the corresponding settings.
#' 
#' @examples
#'
#' @export
# TODO: method is currently set to brms. Add option to pass more parameter without specifying beforehand, for example for the horseshoe prior. 
setup_priors<-function(prior_type,method, par_ratio){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  switch(prior_type,
  weak1={
    switch(method, 
#      MCMCglmm={
#        # R is for the residual variance
#        # G is for the random effect variance
#        ## V is the covariance matrix and nu the degree of belief for the inverse-Wishart, mean vector=alpha_mu, covariance matrix for working parameters = alpha.V)
#        # B is for the fixed effects (mu=expected value and V=covariance matrix, representing the strength of the belief
#        # G_V is set to 1 since only one random effect
#        # R_V is fixed in the binomial setting
#        # priors correspond to an inverse-Gamma distribution with shape and scale parameters equal to 0.01
#        # TODO: check if V of G is set to 1 or to number of levels of random components
#        priors=list(R=list(V=0.5,nu=0.02, fix=TRUE),G=list(G1=list(V=diag(n_rand_levels),nu=0.02, alpha.mu=rep(0, n_rand_levels), alpha.V=diag(n_rand_levels*25^2))) )
      #}, 
      # flat prior
      brms={
        priors=c(set_prior('normal(0,1)'))
      }, 
      {stop("unknown method provided in setup_priors function")}
    )
  }, 
  weak2={
    switch(method, 
      brms={
        priors=c(set_prior('normal(0,10)'))
      },
      {stop("prior only implemented for the brms package")}
    
    )
  },
  horseshoe={
    switch(method,
      brms={
		if(!is.null(par_ratio) | par_ratio=="NULL"){
			horseshoe_call =paste0('horseshoe(df = 1, par_ratio = ', par_ratio, ')')
		} else {
			horseshoe_call ='horseshoe(df = 1)'
		}
      	
        priors=c(set_prior(horseshoe_call, class="b"))
      }, 
      {stop("unknown prior for given method provided in setup_priors function")}
    )
  }, 
  {stop("unknown prior_type provided in setup_priors function.")}
  
  )
  return(priors)
}

#' Get phylogenetic information according to method specification.
#'
#' @description Description.
#'
#' @param tree 
#' @param method
#'
#' @return The phylogenetic information as taken by the given glmm method.
#' 
#' @examples
#'
#' @export

get_phylogeny<-function(tree, method){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  switch(method, 
  MCMCglmm = {
    phylo_info<-MCMCglmm::inverseA(tree,nodes="TIPS", scale=TRUE)$Ainv
#    TODO: check if this is an alternative or outdated
#    V <- ape::vcv(tree)
#    V <- V/max(V)
#    detV <- exp(determinant(V)$modulus[1])
#    V <- V/detV^(1/n)
#    phylo_info <- Matrix(solve(V),sparse=T)
  }, 
  brms={
    # compute variance-covariance matrix A
    phylo_info <- ape::vcv.phylo(tree)
  }, 
  {stop("unknown glmm method provided in get_phylogeny function")}
  
  )
  return(phylo_info)

}


#' Standardize numeric variables to have mean 0 and sd 0.5.
#'
#'
#' @param v Numeric vector. 
#'
#' @return Numeric vector with standardized entries to have mean 0 and 0.5
#' 
#' @examples
#'
#' @export

my_standardize<-function(v){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  nvar<-(v - mean(v, na.rm=TRUE))/(2 *sd(v, na.rm=TRUE))
  return(nvar)
}


##' Computes the number of parameters/coefficents given the fixed and random component. 
##'n.
##'
##' @param fixed Character vector with column names in the data frame that will enter model as fixed variables.
##' @param random Character that denotes the column name of the random component. 
##' @param df Data frame with the input data.
##'
##' @return Integer denoting the number of parameters of the model.
##' 
##'
##' @export

#compute_n_par<-function(fixed, random, df){
#  # assertive tests for each argument
#  
#  # missing values/coercion handling
#  
#  # warnings
#  
#  # code
#  pars<-sapply(df[fixed], function(col){
#    if(is.numeric(col) | is.integer(col)){
#      return(1)
#    }
#    else if (is.factor(col)){
#      return(length(levels(col))-1)
#    }
#  })
#  n_par<-sum(pars)+length(random)+1
#  return(n_par)
#}




#' Prepares the data for fitting a glmm.  
#'
#' @description Calls the build_glmm method after ensuring that all format requirements are fulfilled. 
#'
#' @param df Data frame with the input data.
#' @param phylo Tree object with the viral phylogeny of the samples in the data frame.
#' @param method Character denoting the glmm method/R package to use.
#' @param model Character denoting if the null model (wihtout HLA information) or the full model is built. 
#' @param prior Character denoting the prior for the coefficients of the Bayesian GLMM (method = brms). 
#' @param seed Integer for initializing the random pseudogenerator.
#' @param fixed Character vector with column names of the data frame that should enter the model as fixed variable.
#' @param random Character with the column name of the random component. Currently only one random-intercept model is possible. 
#' @param outcome Character with the column name of the outcome variable. 
#' @param adapt_delta Parameter for the MCMC sampling for the brms method. 
#' @param outcome_type Character denoting if categorical or bernoulli fit should be performed.
#' @param par_ratio Parameter for the horseshoe prior for the brms method.
#'
#' @return A glmm fit
#' 
#'
#' @export

build_model<-function(df, phylo, method, model, prior, seed, fixed, random, outcome, adapt_delta, outcome_type, par_ratio){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code

  # subset the tree accordingly

  # TODO: Note, that it is hard coded that phylo depends on the patient_id
  # though this could be changed to the random variable if set
  phylo <-keep.tip(phylo, tip=as.character(df$patient_id))
  # Check if all factor variables have still contrast
  if(!is.null(fixed)){
    is_factor<-sapply(fixed, function(v,df){is.factor(df[[v]])}, df=df)
    has_contrast<-sapply(fixed[is_factor], function(fac, df){
      # check if factor variable has still contrasts/ alternatively with unique
      return(length(levels(factor(df[[fac]]))) >=2) 
    }, df=df)
    # remove non-contrast variables from model
    if(!all(has_contrast)){
    # log
    print("Warning: following categorical input variable removed from model since no contrast left:")
    print(fixed[is_factor][!has_contrast])
    # delete from model
    del_inds<-sapply(fixed[is_factor][!has_contrast], function(f){which(fixed==f)})
    fixed<-fixed[-del_inds]
    }
    
    # scale non-categorical features 
    is_integer<-sapply(fixed, function(v,df){is.integer(df[[v]])}, df=df)
    integer_features<-fixed[is_integer]
    # standardize all numeric feature, 
    if(length(integer_features)>0){
      df<-df %>% mutate_at(vars(matches(paste0("(", paste0(integer_features, collapse="|"), 
")"))), function(x) {scale(x, center=T, scale=2*sd(x))})
    }
  }
  # Call model fit
  fit<-build_glmm(sample_prior="no", adapt_delta=adapt_delta, outcome=outcome, df=df, phylo=phylo, method=method, model=model, prior=prior, seed=seed, fixed=fixed, random=random, outcome_type = outcome_type, par_ratio = par_ratio)
#  models<-list("mod_prior"=mod_prior, "mod_fit"=mod_fit)
  return(fit)


}
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
# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# open log connection
sink(opt$log)
#if(class(opened_log)="try_error"){
#  stop("sink connection could not be opened")
#}
#tryCatch(sink(opt$log), error=function(e){
#    e$message <- paste0("Error: sink connection could not be opened with the following error msg: \n", e$message)
#    stop(e)
#})

## READ IN ARGUMENTS

# load phylogeny, only needed for random component
phylo<-ape::read.tree(opt$tree)
# load data 
df <-readRDS(opt$data)
stopifnot(!is.null(df$saav_id))
stopifnot(length(unique(df$saav_id))==1)
fixed<-convert_arglist(opt$fixed)
# random component
random<-convert_arglist(opt$random)
outcome<-convert_arglist(opt$outcome)
#
seed<-as.numeric(opt$seed)

#change order of outcome levels s.t. the first level is the amino acid from the reference protein
ref_map<-read.csv(opt$ref_map, skip = 2, header = FALSE)
ref_level <-ref_map %>% filter(as.numeric(V2) ==as.numeric(opt$saav)) %>% select(V1) %>% pull
df<- df %>% mutate(!!outcome := recode(get(outcome), '-' = 'x',.default = levels(get(outcome)))) %>% 
mutate(!!outcome := relevel( get(outcome), tolower(ref_level))) %>%
mutate(!!outcome := factor(get(outcome)))




## SET ATTRIBUTES

print(opt$par_ratio)
switch(opt$model, 
  full = {
    # get all hla allele column names
    hla_column_names<-df %>% select(matches("\\w\\_\\d{2,3}\\_\\d{2,3}")) %>% colnames
    fixed=c(fixed, hla_column_names)
    mod<-build_model(df=df, phylo=phylo, method=opt$method, prior=opt$prior, seed=seed, fixed=fixed, random=random, adapt_delta=opt$adapt_delta,outcome=outcome, model=opt$model, outcome_type= opt$outcome_type, par_ratio =opt$par_ratio)
  }, 
  # clincal and random component but without hla
  null_clin_rand ={
    mod<-build_model(df=df, phylo=phylo, method=opt$method, prior=opt$prior, seed=seed, fixed=fixed, random=random, adapt_delta=opt$adapt_delta,outcome=outcome, model=opt$model, outcome_type= opt$outcome_type,par_ratio =opt$par_ratio)
  }, 
  {
   stop("Error: not implemented model option provided. See help -h.")
  
  }
  )
  
# save model(s)
if(!dir.exists(dirname(opt$output))){
  dir.create(dirname(opt$output), recursive =TRUE)
}
saveRDS(mod, opt$output)
mod_summary <- rbind(summary(mod)$random[[1]], summary(mod)$fixed)

if(!is.null(opt$coeffs)){
	write.csv(mod_summary, opt$coeffs)
}

if(!is.null(opt$fitted)){
  if(!all(is.null(c(opt$protein, opt$model_version)))){
    df<-data.frame(protein=opt$protein, model_version = opt$model_version, saav = opt$saav, model = opt$model, outcome_type = opt$outcome_type, patient_id = mod$data$patient_id, predict(mod), y = mod$data$y)
    write.csv(df, opt$fitted, row.names = FALSE)
  
  }

}

traceback()
# close log
sink()
