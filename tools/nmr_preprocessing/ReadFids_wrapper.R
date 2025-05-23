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
## ------------------------------
library(batch)
library(ggplot2)
library(gridExtra)
library(reshape2)


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
## ------------------------------
if (!runExampleL) {
    argLs <- parseCommandArgs(evaluate = FALSE)
}

sink(argLs$logOut)

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
path <- paste(getwd(), strsplit(directory[1], "/")[[1]][2], sep = "/")


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
write.table(data_matrix, file = argLs$dataMatrix, quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# Sample metadata
write.table(data_sample, file = argLs$sampleMetadata, quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# log file
# write.table(t(data.frame(argLs)), file = argLs$logOut, col.names = FALSE, quote=FALSE)

# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")

argLs

## Ending

cat("\nEnd of 'ReadFids' Galaxy module call: ", as.character(Sys.time()), sep = "")

sink()


options(stringsAsFactors = strAsFacL)

rm(list = ls())
