#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to make a single amino acid variant table out of an aa alignment.
#
# Code developed by Anna Hake
# Date: 2019-05-14
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", "tidyr", "dplyr")

## load libraries
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-a", "--alignment"), type="character", default=NULL, help="A csv file with the cleaned aa alignment", metavar="FILE"),
  make_option(c("-r", "--ref_map"), type="character", default=NULL, help="A map file with the reference to HXB2", metavar="FILE"), 
  make_option(c("-d", "--degree"), type="double", default=0.01, help="The degree to which a variant has to be prevalent in the population"), 
  make_option(c("-o", "--output"), type="character", default=NULL, help="A csv file where the SAAV will be stored."), 
  make_option(c("--outcome_type"), type="character", default= NULL, help="Option for multiclass vs binary model")
)
# ------------------------------------------------------------------------------
# FUNCTIONS


#TODO find better word for degree!

is_SAAV <-function(col,degree){
	# assert that degree is a double in the range 0 to 1
	stopifnot(is.double(degree))
	stopifnot(degree >=0 & degree <= 1)
	# assert that column is not empty
	stopifnot(!is.null(col))
	stopifnot(length(col)>0)
	stopifnot(!all(is.na(col)))
	
	# compute the absolute cutoff to define frequent
	cutoff <- degree * length(col[!is.na(col)])
	# check which AAs are frequent
	are_frequent_aas <- as.data.frame(table(col))$Freq>=cutoff
	# saav is defined if there are at least 2 frequent AAs
	is_saav<-sum(as.numeric(are_frequent_aas))>=2
	return(is_saav)
}

is_frequent_aa <-function(aa, col, degree){
	stopifnot(aa %in% col)
	cutoff <- degree * length(col[!is.na(col)])
	freq_tbl<-as.data.frame(table(col))
	return(freq_tbl %>% filter(col == aa) %>% select(Freq) %>% pull >= cutoff)
}

convert2SAAV_categorical <- function(col, degree){
	stopifnot(is.vector(col))
	saav_col<-sapply(col, convert_AA_categorical, col = col, degree = degree)
	saav_col<-as.factor(saav_col)
	return(saav_col)
}
convert_AA_categorical<-function(aa,col, degree){
	if(is.na(aa)){
		return(NA)
	} else if (is_frequent_aa(aa, col, degree)){
		return(aa)
	} else {
		return('OTHER')
	}


}
make_SAAV_df_categorical<-function(alignment, degree=opt$degree){
 							# remove the patient_id
 							# TODO check that only AAs columns are considered automatically instead of hardcoding
 saav_df <- alignment[,-1] %>% 
 							# select only SAAVs
 							select_if(~sum(!is.na(.)) > 0) %>%
 							select_if(function(x) is_SAAV(x,degree=degree)) %>%
 							# convert to multiclass vector: non-frequent AA are converted to OTHER
 							mutate_all(function(x) {convert2SAAV_categorical(col=x, degree=degree)})
 saav_df<-bind_cols(patient_id = alignment[,1], saav_df )
 return(saav_df)



}

#' Checks if two binary columns are complements.
#'
#' @description Description.
#'
#' @param col_pair integer vector with the indices of the two columns to be tested.
#' @param df data frame where the two columns are stored
#' 
#' @return Boolean denoting if the two binary columns are complements or not. 
#' 
#' @examples
#'
#' @export

is.complement<-function(col_pair, df){
  # assertive tests for each argument
  
  ## indices are within column range
  valid_indices<-sapply(col_pair, function(col_ind){col_ind<=ncol(df)})
  stopifnot(all(valid_indices))
  ## indices are integers
  stopifnot(all(is.integer(col_pair)))
  ## vectors are binary vectors
  valid_entries<-sapply(unlist(df), function(entry){entry %in% c(NA,0,1)})
  stopifnot(all(valid_entries))
  ## col_inds has length 2
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  ## definition complement: union returns a vector full of ones 
  nas<-lapply(col_pair, function(ind){which(is.na(df[, ind]))})
  stopifnot(nas[[1]]==nas[[2]])
  union<-df[-nas[[1]], col_pair[1]] + df[-nas[[1]],col_pair[2]]
  return(all(union==1))
}

#' Removes complements from a data frame of binary columns.
#'
#'
#' @param df data frame with binary columns
#'
#' @return Data frame with binary columns where no complements exist. 
#' 
#' @examples
#'
#' @export

remove_complements<-function(df){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  # generate all possible pairwise combinations
  combs<-combn(seq(ncol(df)),2)
  # for each pairwise combination
  are_complement<-sapply(as.data.frame(combs), is.complement, df)
  if (any(are_complement)){
    complements_comb<-which(are_complement)
    complements_ind<-sapply(complements_comb, function(comb){
    	ones_1<-sum(df[, combs[1,comb]], na.rm=TRUE)
    	ones_2<-sum(df[,combs[2,comb]], na.rm=TRUE)
    	comb_row<-ifelse(ones_1 > ones_2, 1,2)
    	return(combs[comb_row, comb])
    })
    new_df<-as.data.frame(df[,-complements_ind])
    colnames(new_df)<-colnames(df)[-complements_ind]
  } else {
  	new_df<-df
  }

  return(new_df)
}

test_remove_complements<-function(){
	df<-data.frame(c(1,1,0,0,0, NA), c(0,0,1,1,1,NA))
	res_df<-remove_complements(df)
	identical(res_df,df[,-2])
}
#' Computes the position of the alignment to the reference.
#'
#' @description Description.
#'
#' @param alignment_length integer denoting the length of the aligned sequences. 
#' @param ref_map data frame with the reference sequence, its position and the mapped position. 
#'
#' @return A character vector with reference positions. 
#' 
#' @examples
#'
#' @export

get_ref_pos<-function(alignment_length, ref_map){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  ## iterate over all columns
  ## initialize colnames
  new_colnames<-rep(NA, alignment_length)

  for (i in 1:alignment_length){
    ## get the index of the corresponding reference position in the map
    ref_pos_ind<-which(ref_map[,3]==as.character(i))
    ## if there is no corresponding reference position, this is an insertion.
    ## TODO check if this is the correct interpretation. 
    ## or if there is alway the position of the current sequence, but potentially a gap in the ref map [,2]
    if(length(ref_pos_ind)==0){
    ## by induction the precedent cases have been fixed
    # TODO: BUG: if first element is missing, no induction base
      last_pos_entry<-unlist(strsplit(new_colnames[i-1], split="_"))
      ## no insertion present
      if(length(last_pos_entry)==1){
        # insert first insertion
        new_colnames[i]<-paste(last_pos_entry[1], 1, sep="_")
      } else {
        new_colnames[i]<-paste(last_pos_entry[1], as.numeric(last_pos_entry[2])+1, sep="_")
      }
    } else {
      new_colnames[i]<-as.character(ref_map[ref_pos_ind,2])
    }
  }
  return(new_colnames)
}

#' Converts an alignment column into a SAAV column.
#'
#' @description Description.
#'
#' @param col character vector representing the alignment column
#' @param degree double denoting the prevalence cutoff of the variant in the population
#'
#' @return a list of corresponding SAAVs
#' 
#' @examples
#'
#' @export

convert2SAAV_binary<-function(col, degree){
  # assertive tests for each argument

  # missing values/coercion handling
  
  # warnings
  
  # code
  
  # if NA column
  if (all(is.na(col))){
    # return empty list/no variants present
    return(NULL)
  }

      freq_cutoff<-ceiling(degree*length(col[!is.na(col)]))

  # get frequencies of all non NA characters in the columns
	# NA's are not counted 
	all_variants<-as.data.frame(table(col))

  

  # which variants are frequent
  freq_variants_ind<-which(all_variants$Freq>=freq_cutoff)
  # if less than two frequent variants
  if (length(freq_variants_ind)<2){
  	# no variant present, only conservative or no frequent character
    return(NULL)
  }
  # for each frequent variant create SAAV binary column
  SAAV_list<-lapply(all_variants[freq_variants_ind,1], function(var){
    # initialize binary column with zeros
    SAAV_col<-rep(0, length(col))
    # identify occurrences of the variant
    var_inds<-which(col==var)
    # mark the variant occurrences with ones
    SAAV_col[var_inds]<-1
    # identify occurrences of NA
    if(any(is.na(col))){
    	nas <-which(is.na(col))
      SAAV_col[nas]<-NA
    }

    # return the binary vector
    return(SAAV_col)
  })
  # coerce to data frame
  SAAV_df<-as.data.frame(SAAV_list)
  # no null dfs up to here
  
  # set names to the variants
  colnames(SAAV_df)<-all_variants[freq_variants_ind,1]
  # check for complements
#  if (ncol(SAAV_df)>1){
#    SAAV_df<-remove_complements(SAAV_df)
#  }
  return(as.data.frame(SAAV_df)) 
}

#' Creates an SAAV table out of an alignment.
#'
#' @description Description.
#'
#' @param alignment a data frame where apart from the id column, each column is an alignment position. 
#' @param degree double denoting the prevalence cutoff of the variant in the population
#'
#' @return Data frame where each column is a valid SAAV. 
#' 
#' @examples
#'
#' @export

make_SAAV_df_binary<-function(alignment, degree){
  # assertive tests for each argument
  if(!is.double(degree)){
    stop("degree must be a double. See --help.")
  }
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  # convert each column to SAAVs apart from the id column
  SAAV_list<-lapply(alignment[,-1], convert2SAAV_binary, degree)
  # set names
  names(SAAV_list)<-colnames(alignment[,-1])
  no_SAAVs<-which(sapply(SAAV_list, is.null))
  # remove those where no variant was present
  #TODO what if there are no no_SAAVs - bug
  SAAV_list<-SAAV_list[-no_SAAVs]
  new_names<-lapply(seq(SAAV_list), function(pos){
  	pos_name<-names(SAAV_list)[pos]
  	SAAV_names<-colnames(SAAV_list[[pos]])
  	new_SAAV_names<-paste(paste('P', pos_name,sep=''), SAAV_names, sep='_')
  	colnames(SAAV_list[[pos]])<<-new_SAAV_names
  	return(new_SAAV_names)
  })
  # concatenate all SAAVs
  SAAV_df<-as.data.frame(SAAV_list)
  # set colnames correctly
  colnames(SAAV_df)<-unlist(new_names)
  # make table
  SAAV_table<-data.frame(alignment[,1],SAAV_df)
  colnames(SAAV_table)<-c('patient_id', colnames(SAAV_df))
  return(SAAV_table)
}


# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
# read in alignment
alignment<-read.csv(opt$alignment, stringsAsFactors=FALSE)
# read in map file
ref_map<-read.csv(opt$ref_map, stringsAsFactors=FALSE, skip=1)
# trim spaces #TODO where are the spaces, what is stored at column three?
ref_map[,3]<-gsub(" ", "", ref_map[,3])
# change the column name of the patient id
colnames(alignment)[1]<-"patient_id"
# change position names of alignment according map file
colnames(alignment)[-1]<-get_ref_pos(ncol(alignment[,-1]), ref_map)
# create SAAV table
if (opt$outcome_type=="bernoulli"){
  SAAV_df<-make_SAAV_df_binary(alignment, degree=opt$degree)
} else if (opt$outcome_type =="categorical") {
  SAAV_df<-make_SAAV_df_categorical(alignment, degree=opt$degree)
} else {
  stop("Unknown outcome_type option provided")
}

# TODO:log how many SNVS, check that null are only those that should be null
print(dim(SAAV_df))
# convert to long table format
SAAV_df_long<-tidyr::gather(SAAV_df, "saav_id", "y", -patient_id)
# save SAAV_df
write.csv(SAAV_df_long, opt$output, row.names=FALSE)

