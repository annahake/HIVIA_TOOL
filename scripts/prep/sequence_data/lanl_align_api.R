# -----------------------------------------------------------------------------
# This is a script that contains functions to access the HIV Sequence Alignments
# tool: https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html
# Code developed by Anna Hake
# Date: 2019-04-25
#
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# LOAD PACKAGES
package_list<-list('httr', 'rvest')
package_availability<-sapply(package_list, function(package){suppressPackageStartupMessages(library(package, character.only=TRUE, logical.return=TRUE))})
if(!all(package_availability)){
  errors_ind<-which(!package_availability)
  stop(paste('following packages are not available: ', paste(package_list[errors_ind],collapse=','), sep=''))
}
# -----------------------------------------------------------------------------
#TODO: user agent?
# ua<-user_agent("myrepo")

align_api<-function(domain, path){
  url<-modify_url(domain, path=path)
  # API is in text/html
  resp<-GET(url)  # TODO adding user agent
  # parse object response to html
  parsed<-content(resp)
  # check for errors
  if (status_code(resp)!=200){
    stop("HIV Sequence Alignments website is not responding")
    # TODO check if there is a message attribute in parsed if there is an error. 
  }
  
  # make the api an S3 object
  structure(
    list(
      content=parsed,
      path=path,
      domain=domain,
      url=url,
      response=resp
    ), 
    class="align_api"
  )
}

print.align_api <- function(x, ...) {
  cat("<ALIGN API ", x$url, ">\n", sep = "")
  str(x$content)
  invisible(x)
}



check_args<-function(alignment_type, suborganism, region_definition, region, subtype, basetype, year, format, output_file){


  # get all arguments without the function name
#  args_list<-as.list(match.call())[-1]
  args_list<-as.list(environment())
  # check number of arguments

  if (length(args_list)!=9){
     stop("wrong number of arguments. The following 9 arguments have to be set: alignment_type, suborganism, region_definition, region, subtype, basetype, year, format and output_file.")
  }
  # for each argument check if it is a valid option
  mapply(function(arg, arg_name){
    switch(arg_name, 
      alignment_type={
        match.arg(arg, c("ALL", "CON", "COM", "FLT", "REF", "RIP"))
      },
      suborganism={
        match.arg(arg, c("HIV1","HIV2", "SIV"))
      },
      #TODO: user defined not supported yet
      region_definition={
        match.arg(arg, c("predefined"))
      },
      region={
        match.arg(arg, c("ENV", "GAG", "GENOME", "LTR", "NEF", "POL", "REV", "TAT", "VIF", "VPR", "VPU"))
      },
      subtype={
        match.arg(arg, c("All","ALLM", "A-K"))
      },
      basetype={
        match.arg(arg, c("DNA", "PRO"))
      },
      year={
        years<-as.character(c(1997:2017))
        match.arg(arg, years)
      },
      format={
        # TODO: make dictionary (warning 'raw' is the same as 'csv')
        formats<-c("fasta", "table", "mase", "pretty", "outali", "clustal", "phylipi", "philips", "rphylipi", "rphilips", "msf", "rsf", "megai", "megas", "slx", "gde", "GDEFlat", "stockholm", "nexusi", "nexuss", "pir", "gcg", "csv")
        match.arg(arg, formats)
      },
      # TODO check if output file corresponds to output format?
      output_file={
        if (!is.character(arg)){
          stop("'output_file' must be a character.")
        } else {
          if (!nzchar(arg)){
            stop("'output_file' must be a non-empty character")
          }
        }
        arg
      },
      {stop("not a valid argument.")}
    )
  }, args_list, names(args_list))
}

#TODO: userdefined with region radio not supported
set_body<-function(alignment_type, suborganism, region_definition, region, subtype, basetype, year, format){

  # store configs in a list
  body<-list(
    # set hidden variable organism
    "ORGANISM"="HIV",
    # set variable alignment type
    "ALIGN_TYPE"=alignment_type, 
    # set variable suborganism
    "SUBORGANISM"=suborganism,
    # set variable region definition
    "PRE_USER"=region_definition,
    # set variable predefined region
    "REGION"=region,
    # set variable subtype
    "GENO_SUB"=subtype,
    # set variable basetype
    "BASETYPE"=basetype,
    #s et variable year
    "YEAR"=year, 
    # set variable format
    "FORMAT"=format, 
    # set submit type
    "submit"="Get Alignment"
  )
  # return the configs
  return(body)
}

download_alignment<-function(alignment_type, suborganism, region_definition, region, subtype, basetype, year, format, output_file){

  # check if arguments are valid
  # store function call
#  check_args_call<-match.call()
#  # change function name to check_args
#  check_args_call[[1]]<-as.name("check_args")
#  # call check_args with same args and environment as this function
#  valid_args<-eval(check_args_call, sys.frame(sys.parent()))
  
  valid_args<-check_args(alignment_type, suborganism, region_definition, region, subtype, basetype, year, format, output_file)
  # all args are valid 
  
  # open tool website
  align_tool<-align_api("https://www.hiv.lanl.gov", "/content/sequence/NEWALIGN/align.html")
  # TODO try submit_form from rvest package
  # get submit form url
  submit_url_path<-html_nodes(align_tool$content, xpath="//form[@method='post']") %>% html_attr('action')
  # set body params (TODO use match.call)
  body<-set_body(alignment_type, suborganism, region_definition, region, subtype, basetype, year, format)
  # submit form
  submit_resp<-httr::POST(url=file.path(align_tool$domain, submit_url_path), body=body)
  # get download button
  a_node_ind<-which(html_nodes(content(submit_resp),"a")%>% html_text()=="Download this alignment")
  # download url
  download_url<-html_nodes(content(submit_resp),"a")[[a_node_ind]]%>% html_attr("href") 
  # download the sequences to output file
  download_resp<-httr::GET(url=file.path(align_tool$domain, download_url))
  # convert binary output to text
  alignment_text<-rawToChar(content(download_resp))
  # write output to disk
  write(alignment_text, output_file)
}
# TODO will this work if this function is used in parallel (different (unique))naming might be a solution)
get_alignment<-function(alignment_type, suborganism, region_definition, region, subtype, basetype, year, format){
  # define tmp file
  tmp_file="tmp_align_download.fasta"
  # download alignment file
  download_alignment(alignment_type, suborganism, region_definition, region, subtype, basetype, year, format, tmp_file)
  # read in fasta file
  library(seqinr)
  # depending on the basetype set seqtype
  switch(basetype, 
         "DNA"={seqtype<- "DNA"}, 
         "PRO"={seqtype<-"AA"}
         )
  # read in fasta file
  fasta<-read.fasta(tmp_file, seqtype = seqtype, as.string=TRUE)
  # remove tmp file
  file.remove(tmp_file)
  # return fasta object
  return(fasta)
}

# main testing 


