#!/usr/bin/env Rscript
#------------------------------------------------------------------------------- 
# This is R code to select sequences for the HIVIA analyses according to their 
# subtype (only subtype C)
#
# Code developed by Anna Hake
# Date: 2019-05-07
# Update:
#-------------------------------------------------------------------------------
# LIBRARIES
required_packages<-list("optparse", 'seqinr', "dplyr", "tidyr", "xtable")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=T))
})

# ------------------------------------------------------------------------------
# ARGUMENTS HANDLING
option_list=list(
  make_option(c("-f", "--fasta_file"), type='character', default=NULL, help="a fasta file with the nt sequences", metavar="FILE"), 
make_option(c("-s", "--subtype_file"), type="character", default=NULL, help="a tsv file with the subtypes of the nt sequences", metavar="FILE"),
make_option(c("-o", "--output_file"), type="character", default=NULL, help="a fasta file", metavar="FILE"),
make_option(c("-l", "--log_file"), type="character", default=NULL, help="a csv file", metavar="FILE"), 
make_option(c('-t', "--nmer_tbl"), type="character", default=NULL, help="a tex file", metavar="FILE")
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

#' Maps the subject id to the corresponding original name.
#'
#' 
#' @param id character denoting subject id
#' @param full_names character vector with the original sequence names
#'
#' @return The matched orginal sequence name. 
#' 
#' @examples
#'
#' @export

map_id2full_name<-function(id, full_names){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code
  
  # get the index of the match in the original sequence set
  match_index<-grep(id, full_names)
  # get the matched complete original name 
  matched_name<-full_names[match_index]
  return(matched_name)

}
#' Maps the subject IDs to the full sequence names.
#'
#'
#' @param ids character vector with the subject ids
#' @param full_names character vector with the original sequence names
#'
#' @return The matched orginal sequence names. 
#' 
#' @examples
#'
#' @export

get_full_names<-function(ids, full_names){
  # assertive tests for each argument
 
  # is character
  if (!is.character(ids) | !is.character(full_names)){
  	stop('Inputs have to be character vectors.')
  }
  # have same length
  if (length(ids)!= length(full_names)){
  	stop('Matching not possible: different number of names provided.')
  }
  # missing values/coercion handling
  
  # warnings
  
  # code 	
  # map name for every sequence
  mapped_names<-sapply(ids, map_id2full_name, full_names)
  return(mapped_names)
  

}
#' Merge rows according to a specific column.
#'
#' @description Description.
#'
#' @param filter_col character denoting the filter column
#' @param merge_cols character vector denoting the columns to be merged. 
#' @param df1 data frame one
#' @param df2 data frame two
#'
#' @return A single data frame with merged rows. 
#' 
#' @examples
#'
#' @export

merge_cells<-function(filter_col, merge_cols, df1, df2){
  # assertive tests for each argument
  
  # missing values/coercion handling
  
  # warnings
  
  # code 
  in_df2<-sapply(df1[,filter_col], function(val){val%in% df2[,filter_col]})
  in_df1<-sapply(df2[,filter_col], function(val){val%in% df1[,filter_col]})
  if(!all(in_df2)){
  	new_df<-rbind(df1[!in_df2, c(filter_col, merge_cols)], df2[!in_df1, c(filter_col, merge_cols)])
  }
  if(any(in_df2)){
  	merged_rows_list<-lapply(which(in_df2), function(i){
  	df2_ind<-which(df1[i,filter_col]==df2[,filter_col])
  	merged_row<-sapply(merge_cols, function(col){
  		paste(df1[i, col], df2[df2_ind, col], sep='|')
  	})
  	merged_row<-c('IDs'=df1[i, filter_col], merged_row)
  })
  merged_df<-do.call(rbind, merged_rows_list)
  new_df<-rbind(new_df, merged_df)
  }

  return(new_df)
}
# ------------------------------------------------------------------------------
# MAIN

# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)
subtypes<-read.csv(opt$subtype_file, stringsAsFactors=F)
sequences<-read.fasta(opt$fasta_file, seqtype='DNA', as.string=TRUE)
## correct the names returned from comet
#subtypes$name<-map_comet2orig_names(unlist(subtypes$name), names(sequences))
## store the corrected subtype information 
#write.table(subtypes, opt$subtype_file, sep='\t', row.names=F)
# Identify subtype c sequences
# C_seqs<-subtypes %>% filter(subtype=='C') %>% select(name)
C_seqs_ind<-which(subtypes$subtype=='C')
#subtypes$name<-get_full_names(subtypes$name, names(sequences))
C_seqs_names<-subtypes$name[C_seqs_ind]
not_C_seqs_names<-subtypes$name[-C_seqs_ind]
if(!is.null(opt$nmer_tbl) & !(is.null(opt$log_file))){
	# check for long n-nmers
# create tbl 
seqs_df<-tbl_df(sequences) %>% 
  # convert to long format
  gather(key="patient_id", value="sequence") 
# check for n-mers of size i
for (i in 1:10){
  seqs_df <- seqs_df %>% 
  # check if n-mer of size i occurs in sequence
  mutate(!!paste0("n_mer_",i) := grepl(paste0("n{", i, "}"), sequence))
}
# convert to long format
seqs_df<-seqs_df %>% gather(key="n_mer", value="n_mer_occurrence", -patient_id, -sequence)

nmer_df<-seqs_df %>% 
  # for each n-mer size
  group_by(n_mer) %>% 
  # count the occurrences, the percentage and the corresponding patient_ids
  summarise(count=sum(n_mer_occurrence, na.rm=TRUE), perc=sum(n_mer_occurrence, na.rm=TRUE)/n(), ids=paste(patient_id[n_mer_occurrence], collapse=","))

print(xtable(nmer_df, type="latex"), file=opt$nmer_tbl)
#TODO make n-cutoff a global variable
#longNmers_ind<-sapply(sequences, function(seq){grepl("n{1,}", seq)})
longNmers_ind<-rep(FALSE, length(sequences))
longNmers_names<-names(sequences)[longNmers_ind]

if(any(longNmers_ind)){
	if(length(C_seqs_ind!=0)){
		excluded<-merge_cells(filter_col='IDs', merge_cols='Criteria', df1=data.frame('IDs'=longNmers_names, 'Criteria'='N>8'), df2=data.frame('IDs'=not_C_seqs_names, 'Criteria'='not subtype C'))
	} else {
		excluded<-data.frame('IDs'=longNmers_names, 'Criteria'='N>8')
	}
} else {
	if (length(C_seqs_ind!=0)){
		excluded<-data.frame('IDs'=not_C_seqs_names, 'Criteria'='not subtype C')
	} else {
	 excluded<-data.frame()
	
	}

}


colnames(excluded)<-c('PatientID', 'Criteria')
# store excluded sequences 
write.table(excluded, opt$log_file, sep='\t', row.names=F)
selected_names<-intersect(names(sequences)[!longNmers_ind], C_seqs_names)

# store the selected sequences
#write.fasta(sequences[selected_names], file.out=opt$output_file, names=selected_names)


}

write.fasta(sequences[names(sequences) %in% C_seqs_names], file.out=opt$output_file, names=names(sequences[names(sequences) %in% C_seqs_names]))
