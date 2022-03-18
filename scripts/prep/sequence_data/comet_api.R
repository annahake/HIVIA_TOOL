# -----------------------------------------------------------------------------
# This is a script that contains functions to subtype HIV-1 sequences using the
# COMET Tool "https://comet.lih.lu/"
# Code developed by Anna Hake
# Date: 2019-04-24
#
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# LOAD PACKAGES
package_list<-list('httr', 'rvest')
package_availability<-sapply(package_list, function(package){suppressPackageStartupMessages(library(package, character.only=TRUE, logical.return=TRUE))})
if(!all(package_availability)){
  errors_ind<-which(!package_availability)
  stop(paste('Error: following packages are not available: ', paste(package_list[errors_ind],collapse=','), sep=''))
}
# -----------------------------------------------------------------------------
#TODO: user agent?
# ua<-user_agent("myrepo")

comet_api<-function(domain, path){
  url<-modify_url(domain, path=path)
  # API is in text/html
  resp<-GET(url)  # TODO adding user agent
  # parse object response to html
  parsed<-content(resp)
  # check for errors
  if (status_code(resp)!=200){
    stop("Error: Comet website is not responding")
    # TODO check if there is a message attribute in parsed if there is an error. 
  }
  
  # make the api an S3 object
  structure(
    list(
      content=parsed,
      path=path,
      url=url,
      domain=domain,
      response=resp
    ), 
    class="comet_api"
  )
}

print.comet_api <- function(x, ...) {
  cat("<Comet API ", x$url, ">\n", sep = "")
  str(x$content)
  invisible(x)
}

set_body<-function(fastafile){
  # store configs in a list
  body<-list(
    # upload the fasta file to subtype
    fastafile=upload_file(fastafile), 
    # check box that the tool is not used for commercial means
    non_commercial="confirmed", 
    # set the submit type
    submit="submit fasta file"
  )
  # return the configs
  return(body)
}

post_submit_form<-function(body, api){
  # send the submit form via a post request with given configs in the body
  resp<-httr::POST(url=api$url, body=body)
}

get_download_url<-function(body, api){
  # send the submit form
  response<-post_submit_form(body=body, api=api)
  # convert response to html object
  resp_html<- content(response)
  # get the node with download url
  download_url_info<-resp_html %>% html_nodes(xpath='.//meta[@http-equiv="refresh"]') %>% html_attr("content") 
  # extract the download url
  download_url<-strsplit(download_url_info, split="URL=")[[1]][[2]]
  # return the download url
  return(download_url)
}

get_download_page<-function(body, api){
  # get the download url
  download_url<-get_download_url(body=body, api=api)
  # get the download page
  download_page<-httr::GET(url=file.path(api$domain, download_url))
  Sys.sleep(5)
  download_page<-httr::GET(url=file.path(api$domain, download_url))
  # return the download page
  return(download_page)
}

#TODO: if sequence is unassigned, this sample has 5 instead of 4 columns, choose different split might solve this problem
format_results<-function(resp){
  # split the tab delimited response table wrt to the tabs
  split_resp<-lapply(content(resp), function(line){strsplit(line, split="\t")})
  # convert the splitted list to matrix
  resp_mat<-do.call(rbind.fill, do.call(c,split_resp))
  # set colnames 
  colnames(resp_mat)<-unlist(strsplit(colnames(content(resp)), split="\t"))
  #return formatted response matrix
  return(resp_mat)
}

get_subtypes<-function(download_page, api){
  # find the Download link/url
  # get all links
  links<-html_nodes(content(download_page), "a")
  # identify the download link
  download_link_index<-links %>% grep("Download", .)
  # throw error if no download button exists
  if (length(download_link_index)!=1){
    stop('Error: no download button to click')
  }
  # extract relative download url
  download_relative_url<-links[download_link_index] %>% html_attr("href") 
  # TODO check if rather get should be used 
  # download the subtypes by sending a post request to the download url
  results<-httr::POST(url=file.path(api$domain, download_relative_url))
  # convert the tab delimited text object to table
  #subtypes_mat<-format_results(results)
  return(content(results, as='text'))
}

download_comet_subtypes<-function(fastafile, output_file){
  # set the location of the comet tool
  comet_url<-"https://comet.lih.lu/"
  comet<-comet_api(domain=comet_url, path="index.php?cat=hiv1")
  # set configuration for the tool, e.g. upload your fasta file, ticking checkboxes, submit type
  body<-set_body(fastafile)
  # corresponds to clicking the submit button
  resp<-get_download_page(body=body, api=comet)
  # corresponds to clicking the download button
  subtypes<-get_subtypes(download_page=resp, api=comet)
  # downloads the subtypes to the given output file name
  write(subtypes, output_file)  
}
#main

#TEST:
#fastafile<-'/home/anna/Schreibtisch/comet_example.fasta'
#output_file<-'/home/anna/Schreibtisch/comet_results.csv'
#download_comet_subtypes(fastafile, output_file)
