#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 170116_NmrPreprocessing.R
## Manon Martin and Marie Tremblay-Franco

## Preamble
runExampleL <- FALSE

## Options
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## Libraries laoding
suppressPackageStartupMessages(library(ptw))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(reshape2))

# In-house function for argument parsing
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

#Import needed functions
source_local("NmrPreprocessing_script.R")
source_local("DrawFunctions.R")

if (!runExampleL)
#  argLs <- parseCommandArgs(evaluate=FALSE)
argLs <- unlist(parse_args())

sink(argLs$logOut)

print(argLs[["logOut"]])

# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")
print(argLs)

