#!/usr/bin/env Rscript

library(batch) ## parseCommandArgs

library(lme4)     ## mixed model computing
library(Matrix)
library(MASS)
library(lmerTest) ## computing pvalue and lsmeans from results of lme4 package 
library(multtest) ## multiple testing
library(ggplot2)
library(gridExtra)
library(grid)

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}

source_local("mixmodel_script.R")
source_local("diagmfl.R")

argVc <- unlist(parseCommandArgs(evaluate=FALSE))

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

    if(!tesL) {

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

sink(argVc["information"])

cat("\nStart of the '", modNamC, "' Galaxy module call: ", format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

## loading
##--------

datMN <- t(as.matrix(read.table(argVc["dataMatrix_in"],
                                check.names = FALSE,
                                header = TRUE,
                                row.names = 1,
                                sep = "\t")))

samDF <- read.table(argVc["sampleMetadata_in"],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")

varDF <- read.table(argVc["variableMetadata_in"],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")


## checking
##---------

flgF("identical(rownames(datMN), rownames(samDF))", txtC = "Column names of the dataMatrix are not identical to the row names of the sampleMetadata; check your data with the 'Check Format' module in the 'Quality Control' section")
flgF("identical(colnames(datMN), rownames(varDF))", txtC = "Row names of the dataMatrix are not identical to the row names of the variableMetadata; check your data with the 'Check Format' module in the 'Quality Control' section")

#flgF("argVc['fixfact'] %in% colnames(samDF)", txtC = paste0("Required fixed factor '"  , argVc['fixfact'], "' could not be found in the column names of the sampleMetadata"))
flgF("argVc['time']    %in% colnames(samDF)", txtC = paste0("Required time factor '"   , argVc['time']   , "' could not be found in the column names of the sampleMetadata"))
flgF("argVc['subject'] %in% colnames(samDF)", txtC = paste0("Required subject factor '", argVc['subject'], "' could not be found in the column names of the sampleMetadata"))

#flgF("mode(samDF[, argVc['fixfact']]) %in% c('character', 'numeric')", txtC = paste0("The '", argVc['fixfact'], "' column of the sampleMetadata should contain either number only, or character only"))
flgF("mode(samDF[, argVc['time']])    %in% c('character', 'numeric')", txtC = paste0("The '", argVc['time']   , "' column of the sampleMetadata should contain either number only, or character only"))
flgF("mode(samDF[, argVc['subject']]) %in% c('character', 'numeric')", txtC = paste0("The '", argVc['subject'], "' column of the sampleMetadata should contain either number only, or character only"))

flgF("argVc['adjC'] %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")
flgF("argVc['trf'] %in% c('none', 'log10', 'log2')")

flgF("0 <= as.numeric(argVc['thrN']) && as.numeric(argVc['thrN']) <= 1", txtC = "(corrected) p-value threshold must be between 0 and 1")
flgF("argVc['diaR'] %in% c('no', 'yes')")


##------------------------------
## Computation
##------------------------------


varDF <- lmixedm(datMN = datMN,
                     samDF = samDF,
                     varDF = varDF,
                     fixfact    = argVc["fixfact"],
                     time       = argVc["time"],
                     subject    = argVc["subject"],
                     logtr      = argVc['trf'], 
                     pvalCutof  = argVc['thrN'],
                     pvalcorMeth= argVc["adjC"], 
                     dffOption  = "Satterthwaite",
					      visu		= argVc["diaR"], 
                     least.confounded = FALSE,
                     outlier.limit = 3,
					      pdfC     = argVc["out_graph_pdf"],
						  pdfE	   = argVc["out_estim_pdf"]
					      )


##------------------------------
## Ending
##------------------------------


## saving
##--------

varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                          varDF)

write.table(varDF,
            file = argVc["variableMetadata_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

## closing
##--------

cat("\nEnd of '", modNamC, "' Galaxy module call: ",
    as.character(Sys.time()), "\n", sep = "")

sink()

options(stringsAsFactors = strAsFacL)

rm(list = ls())
