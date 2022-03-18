# -----------------------------------------------------------------------------
# This is code to collect the fasta sequences for the HIVImmunoAdapt project
# given the corresponding config file. 
# Code developed by Anna Hake
# Date: 2019-04-10
#
# -----------------------------------------------------------------------------

# PACKAGE HANDLING (TODO install? better way?)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(tibble))

# -----------------------------------------------------------------------------

# ARGUMENT HANDLING
option_list <- list(
	make_option(c("-n", "--protein_name"),
	default=NULL, type="character", 
	help="The protein for which the sequences are collected."),
	make_option(c("-c", '--consensus_cutoff'), action="store",
	default=10, type="integer", 
	help="The prevelance cutoff for the consensus sequence"), 
	make_option(c("-q", '--search_token'), action="store",
	default=NULL, type="character", 
	help="The regular expression to retrieve the sequences from the path"), 
	make_option(c("-o", '--output_file'), action="store",
	default=NULL, type="character", 
	help="The file name where the fetched sequences are stored"), 
	make_option(c("-p", '--search_paths'), action="store",
	default=NULL, type="character", 
	help="comma seperated list of directories where the sequences are stored")
)
# -----------------------------------------------------------------------------

# FUNCTION DEFINITION

#' Convert fasta filename to subject id

convert_file2sequence_name<-function(filename){
	subject_id<-unlist(strsplit(basename(filename), split='\\_'))[1]
	return(subject_id)
}
#' Extracts a sequence from a fasta file matching provided pattern in its 
#' title.
#' 
#' Reads in a fasta file and extracts a sequence which contains a specific 
#' pattern in its title. 
#' @param fasta_file The file name of the fasta sequence. Compressed files 
#' are also possible. (character)
#' @param pattern A regular expression that is used as pattern to identify the 
#' sequence of interest. (character)
#'
#' @return Returns the sequence of interest (list) if pattern is found. 
#' Otherwise NA is returned. The list contains the sequence of class SeqFastadna and the attributes \code{name} and \code{Annot}
#TODO set names correctly!!!!!!
extract_sequence_from_fasta<-function(fasta_file, sequence_pattern){
	seqs<-read.fasta(fasta_file, as.string=TRUE)
	seq_index<-grep(names(seqs), pattern=paste("-", sequence_pattern, sep=""))
	if(length(seq_index)>0){
		return(seqs[seq_index])
	} else {
		return(NA)
	}
}

#' Collects all sequences from the same batch. 
#' 
#' Collects all sequences from the same batch. The batch is given by the 
#' \code{path} argument. There are several compressed fasta files in the path 
#' path directory. The fasta files of interest are identified using the 
#' \code{file_pattern} in the given path. The \code{sequence_pattern} identifies
#' the sequence of interest in a fasta file.
#' 
#' @param path The path where the fasta sequences of this batch are stored. (character)
#' @param file_pattern A regular expression that is used as pattern to identify the 
#' fasta files of interest in the path directory. (character)
#' @param sequence_pattern A regular expression that is used as pattern to identify
#' the sequence of interest in a fasta file.(character)
#'
#' @return Returns a list of all relevant sequences in that batch.
collect_batch_sequences<-function(path, file_pattern, sequence_pattern){

fasta_names_list<-grep(list.files(path=path, full.names=TRUE), pattern=file_pattern, 
                           perl=TRUE, value=TRUE)
seqs_list<-lapply(fasta_names_list, extract_sequence_from_fasta, 
                  sequence_pattern=sequence_pattern)
return(seqs_list)
}

#' Set sequence name.
#' 
#' Changes the sequence name to the subject_id|protein_id|consensus_cutoff|batch_id
#' 
#' @param fasta_entry fasta sequence with attribute name, annotation etc. 
#' @param protein_id	character denoting the corresponding protein name of the sequence. 
#' @param batch_id integer denoting the batch number. 
#'
#' @return Returns the new name of the sequence (character). 
set_sequence_name_rich<-function(fasta_entry, protein_id, batch_id){
	# split the current name at '-cutoff-'. Not working since in first batch cutoff is missing for gag. 
	sequence_name<-unlist(strsplit(attr(fasta_entry, 'name'), split="-cutoff-"))
	# first part is the subject_id
  subject_id<-sequence_name[1]
  # second part is the cutoff used
  consensus_cutoff<-sequence_name[2]
  # name tags of new name
  names_list<-tibble::lst(subject_id, protein_id, consensus_cutoff, batch_id)
  # set name of this sequence
  attr(fasta_entry, 'name')<-paste(names(names_list), names_list , sep="=", collapse="|")
  # return name
  return(attr(fasta_entry, 'name'))
}
#' Set sequence name.
#' 
#' Changes the sequence name to the subject_id
#' 
#' @param fasta_entry fasta sequence with attribute name, annotation etc. 
#'
#' @return Returns the new name of the sequence (character). 
set_sequence_name<-function(fasta_entry){
	# delete everything after the dash
	sequence_name<-gsub("(-cutoff)?-10(\\%)?$", "", attr(fasta_entry, 'name'))
	# first part is the subject_id
  subject_id<-sequence_name
  # return name
  return(subject_id)
}
#' Get sequence names for all sequences.
#' 
#' @param fasta_list list with all sequences from different batches. 
#' @param protein_id	character denoting the corresponding protein name of the sequence. 
#'
#' @return Returns the names of the sequences (character). 
get_sequence_names_rich<-function(fasta_list, protein_id){
	#iterate over all batches
	batch_seqs_names<-lapply(seq(fasta_list), function(batch_id){
		#retrieve names for each sequence in the batch
		new_seq_names<-lapply(fasta_list[[batch_id]], set_sequence_name, protein_id=protein_id, batch_id=batch_id)
	})
	return(batch_seqs_names)
}
#' Get sequence names for all sequences.
#' 
#' @param fasta_list list with all sequences from different batches. 
#'
#' @return Returns the names of the sequences (character). 
get_sequence_names<-function(fasta_list){
	#iterate over all batches
	batch_seqs_names<-lapply(seq(fasta_list), function(batch_id){
		#retrieve names for each sequence in the batch
		new_seq_names<-lapply(fasta_list[[batch_id]], set_sequence_name)
	})
	return(batch_seqs_names)
}


# -----------------------------------------------------------------------------

# MAIN
opt <- parse_args(OptionParser(option_list=option_list))
opt$search_token<-as.character(opt$search_token)
# get search paths
search_paths<-unlist(strsplit(opt$search_paths, ","))
# collect sequences for all batches
sequences_list<-lapply(search_paths, collect_batch_sequences, file_pattern=opt$search_token, sequence_pattern=opt$consensus_cutoff)
# format fasta header of each sequence to subject_id|protein_id|consensus_cutoff|batch_id
sequences_names<-unlist(get_sequence_names(sequences_list))
# format names
sequences_names<-gsub('-', '', sequences_names)
# write sequences to fasta file
write.fasta(do.call(c, sequences_list), names=sequences_names, file.out=opt$output_file, as.string=TRUE)

