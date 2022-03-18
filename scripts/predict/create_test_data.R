#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to
#
# Code developed by Anna Hake
# Date:
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", 'dplyr', 'tidyr')

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
  make_option(c("-s", "--seqs"), type="character", default=NULL, help="A csv file with the aligned cleaned sequences", metavar="FILE"),
  make_option(c('--coreceptor'), type="character", default=NULL, help="A csv file with the predicted coreceptor usage from the g2p coreceptor tool", metavar="FILE"), 
  make_option(c('--clinical'), type="character", default=NULL, help="A csv file with the clinical data", metavar="FILE"), 
  make_option(c('--saav_list'), type="character", default=NULL, help="A string with the final saav_list", metavar="FILE"), 
    make_option(c('-d', '--output_dir'), type="character", default=NULL, help="Path to the output directory", metavar="STRING"), 
  make_option(c('--output_suffix'), type="character", default=NULL, help="Suffix of the output", metavar="STRING"),
  make_option(c("--cd4_cutoff"), default=NULL, help="CD4 cutoff to distinguish chronics from early patients", metavar="INT"),
  make_option(c("-r", "--ref_map"), type="character", default=NULL, help="A map file with the reference to HXB2", metavar="FILE")
  )
  
# ----
# FUNCTIONS

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


#----
# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)

saav_list<-unlist(strsplit(opt$saav, split = ','))


coreceptor_df<-read.csv(opt$coreceptor, stringsAsFactors=FALSE)
hla_df<-read.csv(opt$hla, stringsAsFactors=FALSE)
clinical_df<-read.csv(opt$clinical, stringsAsFactors=FALSE)

alignment<-read.csv(opt$seqs, stringsAsFactors=FALSE)
# read in map file
ref_map<-read.csv(opt$ref_map, stringsAsFactors=FALSE, skip=1)
# trim spaces #TODO where are the spaces, what is stored at column three?
ref_map[,3]<-gsub(" ", "", ref_map[,3])
# change the column name of the patient id
colnames(alignment)[1]<-"patient_id"
# change position names of alignment according map file
colnames(alignment)[-1]<-get_ref_pos(ncol(alignment[,-1]), ref_map)
alignment<-alignment %>% select(one_of(saav_list), patient_id)
alignment_long <- alignment %>% gather(key = saav_id,value = y, -patient_id) 

df_selected <- left_join(alignment_long, coreceptor_df, by='patient_id') %>% left_join(hla_df, by= 'patient_id') %>% left_join(clinical_df, by='patient_id')

df_selected <- df_selected %>% dplyr::filter(onART == "NO", cd4_count > as.numeric(opt$cd4_cutoff))

df_list <- df_selected %>% group_split(saav_id)

lapply(df_list, function(df){
	saav_id<-unique(df$saav_id)
	fname <- file.path(opt$output_dir, saav_id, opt$output_suffix)
	if (!dir.exists(dirname(fname))){
  	dir.create(dirname(fname), recursive=TRUE)
	}
	saveRDS(df, fname)
})

warnings()
traceback()
