#!/usr/bin/Rscript --vanilla --slave --no-site-file

################################################################################################
# WRAPPER FOR tablemerge_script.R (TABLE MERGE)                                                #
#                                                                                              #
# Author: Melanie PETERA                                                                       #
# User: Galaxy                                                                                 #
# Original data: used with tablemerge_script.R                                                 #
# Starting date: 11-05-2015                                                                    #
# V-1: First version of wrapper 
# V-1.1: r-batch removal
#                                                                                              #
#                                                                                              #
# Input files: dataMatrix ; Metadata file                                                      #
# Output files: dataMatrix ; Metadata file                                                     #
#                                                                                              #
################################################################################################

#batch package replacement
parse_args <- function() {
  args <- commandArgs()
  start <- which(args == "--args")[1] + 1
  if (is.na(start)) {
    return(list())
  }
  seq_by2 <- seq(start, length(args), by = 2)
  result <- as.list(args[seq_by2 + 1])
  names(result) <- args[seq_by2]
  return(result)
}

args <- parse_args()

source_local <- function(...){
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	for(i in 1:length(list(...))){source(paste(base_dir, list(...)[[i]], sep="/"))}
}
#Import the different functions
source_local("tablemerge_script.R","RcheckLibrary.R","miniTools.R")


if(length(args) < 4){ stop("NOT enough argument !!!") }


tab.merge(args$dataMatrix_in, args$Metadata_in, args$metatype, args$combined_out)


#delete the parameters to avoid the passage to the next tool in .RData image
rm(args)
