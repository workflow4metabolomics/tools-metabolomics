#!/usr/bin/env Rscript

library(batch) ## parseCommandArgs

# Constants
argv <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=","",argv[grep("--file=",argv)])
prog.name <- basename(script.path)

# Print help
if (length(grep('-h', argv)) >0) {
	cat("Usage:", prog.name,
	    "dataMatrix_in myDataMatrix.tsv",
	    "sampleMetadata_in mySampleData.tsv",
	    "variableMetadata_in myVariableMetadata.tsv",
	    "facC qual",
	    "tesC kruskal",
	    "adjC fdr",
	    "thrN 0.05",
	    "variableMetadata_out myVariableMetadata_out.tsv",
	    "figure figure.pdf",
	    "information information.txt",
		"\n")
	quit(status = 0)
}

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}

source_local("univariate_script.R")

argVc <- unlist(parseCommandArgs(evaluate=FALSE))

##------------------------------
## Initializing
##------------------------------

## options
##--------

strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## packages
##---------

library(PMCMR)

## constants
##----------

modNamC <- "Univariate" ## module name

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

cat("\nStart of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

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

tesC <- argVc["tesC"]

## checking
##---------

flgF("identical(rownames(datMN), rownames(samDF))", txtC = "Column names of the dataMatrix are not identical to the row names of the sampleMetadata; check your data with the 'Check Format' module in the 'Quality Control' section")
flgF("identical(colnames(datMN), rownames(varDF))", txtC = "Row names of the dataMatrix are not identical to the row names of the variableMetadata; check your data with the 'Check Format' module in the 'Quality Control' section")

flgF("argVc['facC'] %in% colnames(samDF)", txtC = paste0("Required factor of interest '", argVc['facC'], "' could not be found in the column names of the sampleMetadata"))
flgF("mode(samDF[, argVc['facC']]) %in% c('character', 'numeric')", txtC = paste0("The '", argVc['facC'], "' column of the sampleMetadata should contain either number only, or character only"))

flgF("!(tesC %in% c('ttest', 'wilcoxon')) || (mode(samDF[, argVc['facC']]) == 'character' && length(unique(samDF[, argVc['facC']])) == 2)", txtC = paste0("For 'ttest' and 'wilcoxon', the chosen factor column ('", argVc['facC'], "') of the sampleMetadata should contain characters with only two different classes"))
flgF("!(tesC %in% c('anova', 'kruskal')) || (mode(samDF[, argVc['facC']]) == 'character' && length(unique(samDF[, argVc['facC']])) > 2)", txtC = paste0("For 'anova' and 'kruskal', the chosen factor column ('", argVc['facC'], "') of the sampleMetadata should contain characters with at least three different classes"))
flgF("!(tesC %in% c('pearson', 'spearman')) || mode(samDF[, argVc['facC']]) == 'numeric'", txtC = paste0("For 'pearson' and 'spearman', the chosen factor column ('", argVc['facC'], "') of the sampleMetadata should contain numbers only"))

flgF("argVc['adjC'] %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")

flgF("0 <= as.numeric(argVc['thrN']) && as.numeric(argVc['thrN']) <= 1",
     txtC = "(corrected) p-value threshold must be between 0 and 1")


##------------------------------
## Computation
##------------------------------


varDF <- univariateF(datMN = datMN,
                     samDF = samDF,
                     varDF = varDF,
                     facC = argVc["facC"],
                     tesC = tesC,
                     adjC = argVc["adjC"],
                     thrN = as.numeric(argVc["thrN"]),
                     pdfC = argVc["figure"])


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

cat("\n\n\n============================================================================")
cat("\nAdditional information about the call:\n")
cat("\n1) Parameters:\n")
print(cbind(value = argVc))

cat("\n2) Session Info:\n")
sessioninfo <- sessionInfo()
cat(sessioninfo$R.version$version.string,"\n")
cat("Main packages:\n")
for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
cat("Other loaded packages:\n")
for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")

cat("============================================================================\n")

sink()

options(stringsAsFactors = strAsFacL)

rm(list = ls())
