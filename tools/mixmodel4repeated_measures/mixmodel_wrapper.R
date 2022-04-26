#!/usr/bin/env Rscript

library(grid)
library(lme4)     ## mixed model computing
library(Matrix)
library(MASS)
library(lmerTest) ## computing pvalue and lsmeans from results of lme4 package
library(multtest) ## multiple testing
library(ggplot2)
library(gridExtra)

source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}

create_parser <- function() {
  parser <- optparse::OptionParser()

  parser <- optparse::add_option(
    parser,
    c("--information"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--dataMatrix_in"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--sampleMetadata_in"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--variableMetadata_in"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--time"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--subject"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--thrN"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--dff"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--fixfact"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--trf"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--adjC"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--diaR"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--out_graph_pdf"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--out_estim_pdf"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  parser <- optparse::add_option(
    parser,
    c("--variableMetadata_out"),
    type = "character", ## This must be changed
    action = "store",
    # default = It would be good to have default values
    help = "Not redacted yet"
  )
  return(parser)
}

args <- optparse::parse_args(create_parser())

source_local("mixmodel_script.R")
source_local("diagmfl.R")


##------------------------------
## Initializing
##------------------------------

## options
##--------

strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## constants
##----------

modNamC <- "mixmodel" ## module name

topEnvC <- environment()
flagC <- "\n"

## functions
##----------

flgF <- function(tesC,
                 envC = topEnvC,
                 txtC = NA) { ## management of warning and error messages

    tesL <- eval(parse(text = tesC), envir = envC)

    if (!tesL) {

        sink(NULL)
        stpTxtC <- ifelse(is.na(txtC),
                          paste0(tesC, " is FALSE"),
                          txtC)

        stop(stpTxtC,
             call. = FALSE)

    }

} ## flgF

## log file
##---------

sink(args$information)

cat("\nStart of the '", modNamC, "' Galaxy module call: ", format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")
cat("\nParameters used:\n\n")
print(args)
cat("\n\n")

## loading
##--------

datMN <- t(as.matrix(read.table(args$dataMatrix_in,
                                check.names = FALSE,
                                header = TRUE,
                                row.names = 1,
                                sep = "\t")))

samDF <- read.table(args$sampleMetadata_in,
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")

varDF <- read.table(args$variableMetadata_in,
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")


## checking
##---------

flgF("identical(rownames(datMN), rownames(samDF))", txtC = "Column names of the dataMatrix are not identical to the row names of the sampleMetadata; check your data with the 'Check Format' module in the 'Quality Control' section")
flgF("identical(colnames(datMN), rownames(varDF))", txtC = "Row names of the dataMatrix are not identical to the row names of the variableMetadata; check your data with the 'Check Format' module in the 'Quality Control' section")

flgF("args$time    %in% colnames(samDF)", txtC = paste0("Required time factor '", args$time, "' could not be found in the column names of the sampleMetadata"))
flgF("args$subject %in% colnames(samDF)", txtC = paste0("Required subject factor '", args$subject, "' could not be found in the column names of the sampleMetadata"))

flgF("mode(samDF[, args$time])    %in% c('character', 'numeric')", txtC = paste0("The '", args$time, "' column of the sampleMetadata should contain either number only, or character only"))
flgF("mode(samDF[, args$subject]) %in% c('character', 'numeric')", txtC = paste0("The '", args$subject, "' column of the sampleMetadata should contain either number only, or character only"))

flgF("args$adjC %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")
flgF("args$trf %in% c('none', 'log10', 'log2')")

flgF("0 <= as.numeric(args$thrN) && as.numeric(args$thrN) <= 1", txtC = "(corrected) p-value threshold must be between 0 and 1")
flgF("args$diaR %in% c('no', 'yes')")


##------------------------------
## Formating
##------------------------------

if (args$dff == "Satt") {
  dffmeth <- "Satterthwaite"
} else {
  dffmeth <- "Kenward-Roger"
}


##------------------------------
## Computation
##------------------------------


varDF <- lmixedm(datMN = datMN,
                     samDF = samDF,
                     varDF = varDF,
                     fixfact     = args$fixfact,
                     time        = args$time,
                     subject     = args$subject,
                     logtr       = args$trf,
                     pvalCutof   = args$thrN,
                     pvalcorMeth = args$adjC,
                     dffOption   = dffmeth,
                     visu        = args$diaR,
                     least.confounded = FALSE,
                     outlier.limit = 3,
                     pdfC        = args$out_graph_pdf,
                     pdfE        = args$out_estim_pdf
                     )


##------------------------------
## Ending
##------------------------------


## saving
##--------

varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                          varDF)

write.table(varDF,
            file = args$variableMetadata_out,
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

## closing
##--------

cat("\nEnd of '", modNamC, "' Galaxy module call: ",
    as.character(Sys.time()), "\n", sep = "")
cat("\nInformation about R (version, Operating System, attached or loaded packages):\n\n")
sessionInfo()

sink()

options(stringsAsFactors = strAsFacL)

rm(list = ls())
