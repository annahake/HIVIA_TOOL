#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to generate k stratified CV folds and save the corresponding 
# test and train data given a data frame, a seed, and a group variable for the 
# stratification 
#
# Code developed by Anna Hake
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", 'dplyr')

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c('-i', "--input"), type="character", default=NULL, help="An rds file with the data frame.", metavar="FILE"), 
  make_option(c('-d', '--output_dir'), type="character", default=NULL, help="Directory name of the output", metavar="DIR"),
  make_option(c('-k', "--kfold"), type="integer", default=NULL, help="Number of split folds", metavar="FILE"),
  make_option(c('-s', '--seed'), type="integer", default=NULL, help="Seed for sampling", metavar="DIR"),
  make_option(c('--group'),type="character", default=NULL, help="variable name wrt perform stratification", metavar="STRING"),
  make_option(c('--log'), type="character", default=NULL, help="A file for logging", metavar="FILE")
  # TODO: Add comma if more than one option is used. 
)
# ------------------------------------------------------------------------------
# FUNCTIONS

#' Create stratified cross-validation folds.
#'
#' @description Given a data frame, the number of folds, a random seed and a group variable according to stratification should be performed, the function returns a vector with the corresponding assigned stratified fold indices. This function is adopted from the caret package function createFolds.
#'
#' @param data data frame with the data 
#' @param kfold integer denoting the number of folds in the cross-validation
#' @param seed seed for repeatability
#' @param group column name according to which the stratification should be performed
#'
#' @return A vector of length of the rows of the data frame where each entry denotes the assignment to one of the k folds. 
#' 
#' @export

create_stratifiedCV_folds<-function(data, kfold, seed, group){

  # assertive tests for each argument

  ## amount of folds has to be an integer 
  stopifnot(is.integer(kfold))
  ## seed has to be an integer
  stopifnot(is.integer(seed))
  ## data has to be not empty
  stopifnot(!is.null(data))
  ## group has to be a column name in data
  stopifnot(group %in% colnames(data))
  
  # set some attributes
  ## size of data samples
  N =  nrow(data)
  ## number of folds
  k = kfold
  ## outcome variable according to which stratification should be performed
  y = as.factor(unlist(data[, group]))
  ## set seed
  set.seed(seed)
  
  ## k has to be in the range k >=1 <=N
  stopifnot(k >=1)
  stopifnot(k <= length(y))
  ## group variable should have a least two levels
  stopifnot(length(levels(y))>=2)
  
  
  #TODO missing values/coercion handling
  # if y_i=NA, set folds_ind_i to NA
  
  # warnings
  
  # code
  
  ## number of samples per fold
  fold_size<-floor(N/k)
  
  ## initialize fold_inds vector
  fold_inds <- rep(NA, length(y))
  ## get the number of samples per class
  class_size<-table(y)

  
  if (k < N){
    ## for each class 
    for (i in 1:length(class_size)){
      ## samples of this class per fold
      min_class_size<-class_size[i]%/% k
      # if at least one sample of this class per fold     
      if (min_class_size > 0){
        remaining<-class_size[i] %% k
        # assign the minimum amount of samples to each class
        assign_vector<-rep(1:k, min_class_size)
        # assign the remaining samples to each class_size
        if (remaining >0){
          assign_vector<-c(assign_vector, sample(1:k, remaining))
        }
        fold_inds[which(y == names(class_size)[i])] <- sample(assign_vector)
     } else {
      # less class samples than folds, distribute randomly
      fold_inds[which(y == names(class_size)[i])] <- sample(1:k, size=class_size[i])
     }
    }
  } else {
  # k == N, that is one sample per fold, randomly assign each sample to a fold
    fold_inds<-seq(along.with=data$group)
    fold_inds<-sample(fold_inds)
  }
  return(fold_inds)

}
# TODO Is there a reason to store all kfold of splits into train and test sets or is the fold indices vector sufficient. 
#' 
#'
#' @description Save stratified split folds to given directory as rds files.
#'
#' @param data data frame with the data
#' @param folds vector of size of the rows of the data with assigned fold indices
#' @param dir output directory where split folds should be stored
#'
#' @return None.
#' 
#'
#' @export

save_folds<-function(data, folds, dir){
  # assertive tests for each argument
  stopifnot(length(folds)==nrow(data))
  # missing values/coercion handling
  folds<-factor(as.character(folds))
  # warnings
  
  # code
  
  ## for each fold, divide data into 
  for (k in levels(folds)){
    ## train data and 
    k_train <-  data %>% filter(folds!=k)
    ## test data
    k_test  <-  data %>% filter(folds==k)
    ## create subdirectory
    output_dir<-file.path(dir, paste0("fold_", k))
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    ## save train and test data
    saveRDS(k_train, file=file.path(output_dir, "train.rds"))
    saveRDS(k_test, file=file.path(output_dir, "test.rds"))
  }
}
  

# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
## read in data frame
df<-readRDS(opt$input)
## read in k fold
k<-as.integer(opt$kfold)
## read in seed
s<-as.integer(opt$seed)
## read in column vector according to which stratification should be performed
group<-as.character(opt$group)
# group has to be a variable in the df
stopifnot(group %in% colnames(df))

## generate folds
kfolds<-create_stratifiedCV_folds(data=df, kfold=k, seed=s, group=group)
# TODO: currently the train and test data is stored, storing the indices directly would be the lighter, more efficient version. 
# save the train and test folds
save_folds(data=df, folds=kfolds, dir=opt$output_dir)
