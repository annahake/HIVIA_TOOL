required_packages<-list("httr", "rvest")

## load libraries
# TODO: add functions for loading libraries in loadlibraries.R/function library
loaded_libraries<-lapply(required_packages, function(package){
  # TODO check if installed, if not raise error/error message
  suppressMessages(library(package, character.only=TRUE))
})
# TODO: check out_ https://cran.r-project.org/web/packages/httr/vignettes/api-packages.html


#' Downloads alignment from download page. 
#'
#' Submits a request form to the download page https://www.hiv.lanl.gov/cgi-bin/common_code/MultiFormatDownload.cgi with given parameters and downloads the corresponding alignment. 
#' 
#' @param post_form The submit_form to download the alignment. 
#' @param output_format A character denoting in which format the alignment will be downloaded. 
#' @param output_dir A character denoting the directory to save the alignment.
#' @param id A character containing a fixed part of the output_filename.
#' @param lanl_url A character denoting the url the main service LANL. 
#'

download_alignment<-function(post_form, output_format, output_dir,id, lanl_url){
  #
  program_fname <- basename(post_form$fields$TableFile$value)
  basetype <- unlist(strsplit(program_fname, "codonalign"))[1]
  output_filename<-file.path(output_dir, paste(paste(id,basetype,"codonalign", sep="_"), output_format, sep="."))
  httr::POST(url = file.path(lanl_url, post_form$url), body = list("TableFile" = post_form$fields$TableFile$value, format=output_format), enctype=post_form$enctype, write_disk(output_filename, overwrite = TRUE))
}

#' Codon alignment using the LANL CodonAlign Tool 
#'
#' Submits a request form to the Codon Alignment Tool 'https://www.hiv.lanl.gov/content/sequence/CodonAlign/codonalign.html' with given parameters and downloads the corresponding alignment. 
#' 
#' @param dna_alignment Depending on the parameter input_type, a file with the dna alignment or directly a character vector with the sequences in fasta format. The following page contains all supported alignment formats: https://www.hiv.lanl.gov/content/sequence/HelpDocs/SEQsamples.html
#' @param frame  A character denoting the starting reading frame to chose. The frame is based on the first (master) sequence in the alignment. Possible options are '1', '2','3', 'all', 'program'. If option 'all' is selected, all possible alignments are computed and stored. Option 'program' lets the program decide on the best starting frame with respect of the longest open reading frame.
#' @param compensatingNum An integer denoting how many codons are allowed for frameshift compensation. A compensating mutation is an insertion or deletion (indel) that can be combined with another nearby indel to produce a multiple of 3 and thus preserve an intact reading frame. 
#' @param input_type A character denoting the type of alignment provided. Options are 'file' and 'text'. 
#' @param output_format A character denoting the format of the output alignment. Options are 'fasta', 'table', 'mase', 'pretty', 'oal', 'rsf', 'gde', 'pir', 'gcg', 'clust', 'phylips', 'phylipi', 'msf', 'megas', 'megai', 'nexuss', 'nexusi', 'stockholm', 'slx'.
#' @param output_seqs_type A chacter denoting the type of the sequences in the alignment. Options are 'DNA', 'AA', or 'all'. If 'all' is selected, two alignments are downloaded. 
#' @param output_dir A character denoting the output directory. Note, that if the option 'all' is selected somewhere, several alignments will be downloaded and thus a directory is expected. 
#' @param id A character containing a fixed part of the output_filename.
#' @return none, since alignments are stored directly to disk
download_codon_align<-function(dna_alignment, frame, compensating_num, input_type, output_format, output_seqs_type, output_dir, id){
  # TODO: Check all parameters, throw error if not valid parameters
  #check_parameter(dna_alignment, frame, compensatingNum, input_type, output_format, output_seqs_type, output_file)
  # TODO: set session parameter like user agen 
  
  # TODO: add log file for cookie, retrieved date, settings etc
  # set parameters of the form
  print(dna_alignment)
  query_list = list(
  'frame' = frame,
  'compensatingNum' = compensating_num
  )
  switch(input_type, 
         file={query_list[['QueryFilehandle']]<-upload_file(dna_alignment)},
         text={query_list[['seq_input']]<-dna_alignment}
         )
  
  # LANL main website
  lanl_url<-'https://www.hiv.lanl.gov'
  # Codon Align tool website
  codon_align_url<-'https://www.hiv.lanl.gov/content/sequence/CodonAlign/codonalign.html'
  # read in as html object
  codon_align_html<-read_html(codon_align_url)
  # get all submit forms on this webpage
  forms<-html_form(codon_align_html)
  
  #retrieve POST form by searching for POST in forms
  post_form_id<-grep('method = "POST"', forms)
 
  #get post form
  align_submit_form<-forms[[post_form_id]]
  
  # post submit form
  download_page <- httr::POST(url = file.path(lanl_url, align_submit_form$url), body = query_list)
  # get all html forms of the download page
  download_page_forms<-html_form(content(download_page))
  
  # find all POST forms instead of hardcoding
  post_form_ids<-grep('method = "POST"', download_page_forms)
  
  # if seq_type is 'DNA' select only 'DNA' download forms
  #grep('ntcodonalign', download_page_forms)

  # if seq_type is 'AA' select only 'AA' download forms
  # else select all post forms
  
  # download for each frame option and each seq type
  last_form<-length(download_page_forms)
  #TODO: at the moment only output_seqs_type='all' is implemented
  downloaded<-lapply(download_page_forms[1:last_form], download_alignment, output_format=output_format, output_dir=output_dir, id=id,lanl_url=lanl_url)
}

