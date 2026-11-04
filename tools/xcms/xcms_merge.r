#!/usr/bin/env Rscript

# Import the different functions
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("xcms", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n")

args <- parseCommandArgs(evaluate = FALSE) # interpretation of arguments given in command line as an R list of objects


cat("\tXSET MERGING...\n")

mergeXDataReturn <- mergeXData(args)
xdata <- mergeXDataReturn$xdata
# singlefile <- mergeXDataReturn$singlefile
md5sumList <- mergeXDataReturn$md5sumList
chromTIC <- mergeXDataReturn$chromTIC
chromBPI <- mergeXDataReturn$chromBPI

cat("\tGET SAMPLE METADATA\n")
# Create a sampleMetada file
sampleNamesList <- getSampleMetadata(xdata = xdata, sampleMetadataOutput = "sampleMetadata.tsv")

cat("\n\n")

cat("\tMsExperiment OBJECT INFO\n")
print(MsExperiment::sampleData(xdata))
print(xdata)
cat("\n\n")
cat("\tSpectra INFO\n")
print(spectra(xdata))
cat("\n\n")

cat("\tSAVE RData\n")
# saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xdata", "md5sumList", "sampleNamesList") # zipfile, singlefile removed
save(list = objects2save[objects2save %in% ls()], file = "merged.RData")
