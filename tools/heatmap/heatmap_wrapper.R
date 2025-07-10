#!/usr/bin/env Rscript

library(batch) ## parseCommandArgs

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}

source_local("heatmap_script.R")

argVc <- unlist(parseCommandArgs(evaluate=FALSE))


##------------------------------
## Initializing
##------------------------------

## options
##--------

strAsFacL <- options()[["stringsAsFactors"]]
options(stringsAsFactors=FALSE)

## constants
##----------

modNamC <- "Heatmap" ## module name

## log file
##---------

sink(argVc["information"])

cat("\nStart of the '", modNamC, "' module: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

## loading
##--------

proMN <- t(as.matrix(read.table(argVc["dataMatrix_in"],
                                check.names = FALSE,
                                header = TRUE,
                                row.names = 1,
                                sep = "\t")))

obsDF <- read.table(argVc["sampleMetadata_in"],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")

feaDF <- read.table(argVc["variableMetadata_in"],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")

## adding default parameter values
##--------------------------------


if(!("corMetC" %in% names(argVc)))
    argVc["corMetC"] <- "pearson"
if(!("aggMetC" %in% names(argVc)))
    argVc["aggMetC"] <- "ward"
if(!("colC" %in% names(argVc)))
    argVc["colC"] <- "blueOrangeRed"
if(!("scaL" %in% names(argVc)))
    argVc["scaL"] <- "TRUE"
if(!("cexN" %in% names(argVc)))
    argVc["cexN"] <- "0.8"

## checking
##---------

if(as.numeric(argVc["cutSamN"]) > nrow(proMN))
    stop("Number of sample clusters must be inferior to the number of samples")
if(as.numeric(argVc["cutVarN"]) > ncol(proMN))
    stop("Number of variable clusters must be inferior to the number of variables")

## printing arguments
##-------------------

cat("\nArguments used:\n\n")
argMC <- as.matrix(argVc)
colnames(argMC) <- "value"
argDatVl <- grepl("\\.dat$", argVc) ## discarding dataset file names
if(sum(argDatVl))
    argMC <- argMC[!argDatVl, , drop = FALSE]
print(argMC)


##------------------------------
## Computation
##------------------------------


heaLs <- heatmapF(proMN = proMN,
                  obsDF = obsDF,
                  feaDF = feaDF,
                  disC = argVc["disC"],
                  cutSamN = as.numeric(argVc["cutSamN"]),
                  cutVarN = as.numeric(argVc["cutVarN"]),
                  fig.pdfC = argVc["figure"],
                  corMetC = argVc["corMetC"],
                  aggMetC = argVc["aggMetC"],
                  colC = argVc["colC"],
                  scaL = as.logical(argVc["scaL"]),
                  cexN = as.numeric(argVc["cexN"]))


##------------------------------
## Ending
##------------------------------


## saving
##-------

proDF <- cbind.data.frame(dataMatrix = colnames(heaLs[["proMN"]]),
                          as.data.frame(t(heaLs[["proMN"]])))
write.table(proDF,
            file = argVc["dataMatrix_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

obsDF <- cbind.data.frame(sampleMetadata = rownames(heaLs[["obsDF"]]),
                          heaLs[["obsDF"]])
write.table(obsDF,
            file = argVc["sampleMetadata_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

feaDF <- cbind.data.frame(variableMetadata = rownames(heaLs[["feaDF"]]),
                          heaLs[["feaDF"]])
write.table(feaDF,
            file = argVc["variableMetadata_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

## Ending
##-------

cat("\nEnd of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")

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
