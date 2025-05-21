#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 08122016_ReadFids_wrapper.R
## Manon Martin
## manon.martin@uclouvain.be

## ======================================================
## ======================================================
# Preamble
## ======================================================
## ======================================================

runExampleL <- FALSE


## ------------------------------
## Options
## ------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)
options(warn = 1)

## ------------------------------
## Libraries laoding
##------------------------------
# library(batch) 
suppressPackageStartupMessages(library(ggplot2)) # nice plots
suppressPackageStartupMessages(library(gridExtra)) # nice plots
suppressPackageStartupMessages(library(reshape2)) # data manipulation
suppressPackageStartupMessages(library(stringr)) # string of characters manipulation

# This function retrieve the raw file in the working directory
#   - if zipfile: unzip the file with its directory tree
#   - if singlefiles: set symlink with the good filename
retrieveRawfileInTheWorkingDir <- function(singlefile, zipfile) {
  if (!is.null(singlefile) && (length("singlefile") > 0)) {
    for (singlefile_sampleName in names(singlefile)) {
      singlefile_galaxyPath <- singlefile[[singlefile_sampleName]]
      if (!file.exists(singlefile_galaxyPath)) {
        error_message <- paste("Cannot access the sample:", singlefile_sampleName, "located:", singlefile_galaxyPath, ". Please, contact your administrator ... if you have one!")
        print(error_message)
        stop(error_message)
      }
      
      file.symlink(singlefile_galaxyPath, singlefile_sampleName)
    }
    directory <- "."
  }
  if (!is.null(zipfile) && (zipfile != "")) {
    if (!file.exists(zipfile)) {
      error_message <- paste("Cannot access the Zip file:", zipfile, ". Please, contact your administrator ... if you have one!")
      print(error_message)
      stop(error_message)
    }
    
    # unzip
    suppressWarnings(unzip(zipfile, unzip = "unzip"))
    
    # get the directory name
    filesInZip <- unzip(zipfile, list = TRUE)
    directories <- unique(unlist(lapply(strsplit(filesInZip$Name, "/"), function(x) x[1])))
    directories <- directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
    directory <- "."
    if (length(directories) == 1) directory <- directories
    
    cat("files_root_directory\t", directory, "\n")
  }
  return(directory)
}

# In-house function for argument parsing instead of the R batch library)
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

# R script call
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
# Import the different functions
source_local("ReadFids_script.R")
source_local("DrawFunctions.R")
## ------------------------------
## Errors ?????????????????????
## ------------------------------


## ------------------------------
## Constants
## ------------------------------
topEnvC <- environment()
flagC <- "\n"


## ------------------------------
## Script
##------------------------------
if(!runExampleL)
  argLs <- unlist(parse_args())

  ## Outputs
logOut <- argLs[["logOut"]]
nomGraphe <- argLs[["graphOut"]]

sink(logOut, append = TRUE)

print(argLs)

## ======================================================
## ======================================================
## Parameters Loading
## ======================================================
## ======================================================

## Inputs
# Path
## Bruker FIDs
fileType <- "Bruker"
zipfile <- argLs[["fidzipfile"]]

directory <- unzip(zipfile, list = F)

# filesInZip <- unzip(zipfile, list = TRUE)

# path1 <- paste(getwd(), strsplit(zipfile1[1], "/")[[1]][2], sep = "/")

path <- paste(getwd(),strsplit(directory[1], "/")[[1]][2], sep = "/")

zipfile1 <- rawFilePath[[1]]
print(zipfile1)

path1 <- paste(paste(getwd(), zipfile1, sep = "/"), "/", sep = "")
print(path1)

path <- paste(paste(getwd(), strsplit(directory, "/")[[1]][2], sep = "/"), "/", sep = "")
print(path)


  # other inputs from ReadFids
l <- argLs[["title_line"]]
subdirs <- argLs[["subdirectories"]]
dirs.names <- argLs[["dirs_names"]]


# Outputs
# dataMatrix <- argLs[["dataMatrix"]]
# sampleMetadata <- argLs[["sampleMetadata"]]
logOut <- argLs[["logOut"]]
nomGraphe <- argLs[["graphOut"]]



## Checking arguments
## -------------------
error.stock <- "\n"

if (length(error.stock) > 1) {
    stop(error.stock)
}



## ======================================================
## ======================================================
## Computation
## ======================================================
## ======================================================
sink(logOut, append = TRUE)

if (length(warnings()) > 0) { # or !is.null(warnings())
    print("something happened")
}

## Starting
cat("\nStart of 'ReadFids' Galaxy module call: ", as.character(Sys.time()), "\n\n", sep = "")

outputs <- ReadFids(path = path, l = l, subdirs = subdirs, dirs.names = dirs.names)
data_matrix <- outputs[["Fid_data"]] # Data matrix
data_sample <- outputs[["Fid_info"]] # Sample metadata

pdf(nomGraphe, onefile = TRUE, width = 13, height = 13)
title <- "Raw FID data"
DrawSignal(data_matrix,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
)
invisible(dev.off())

## ======================================================
## ======================================================
## Saving
## ======================================================
## ======================================================

# Data matrix
write.table(data_matrix, file = argLs[["dataMatrix"]], quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# Sample metadata
write.table(data_sample, file = argLs$sampleMetadata, quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# log file
# write.table(t(data.frame(argLs)), file = argLs$logOut, col.names = FALSE, quote=FALSE)

# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")

# Reproductibility
sessionInfo()
## Ending

cat("\nEnd of 'ReadFids' Galaxy module call: ", as.character(Sys.time()), sep = "")

sink()

options(stringsAsFactors = strAsFacL)

rm(list = ls())

