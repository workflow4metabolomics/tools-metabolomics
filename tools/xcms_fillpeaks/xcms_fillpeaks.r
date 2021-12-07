#!/usr/bin/env Rscript

# ----- LOG FILE -----
log_file <- file("log.txt", open = "wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("xcms", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = F, quote = F, sep = "\t")

cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

#saving the specific parameters
method <- "FillChromPeaks"

if (!is.null(args$convertRTMinute)) convertRTMinute <- args$convertRTMinute
if (!is.null(args$numDigitsMZ)) numDigitsMZ <- args$numDigitsMZ
if (!is.null(args$numDigitsRT)) numDigitsRT <- args$numDigitsRT
if (!is.null(args$intval)) intval <- args$intval
if (!is.null(args$naTOzero)) naTOzero <- args$naTOzero

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
load(args$image)
if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS 2.*")

#Verification of a group step before doing the fillpeaks job.
if (!hasFeatures(xdata)) stop("You must always do a group step after a retcor. Otherwise it won't work for the fillpeaks step")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile


cat("\n\n")

# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tCOMPUTE\n")

cat("\t\t\tFilling missing peaks using default settings\n")
# clear the arguement list to remove unexpected key/value as singlefile_galaxyPath or method ...
args <- args[names(args) %in% slotNames(do.call(paste0(method, "Param"), list()))]

fillChromPeaksParam <- do.call(paste0(method, "Param"), args)
print(fillChromPeaksParam)

# back compatibility between xcms-3.0.0 and xcms-3.5.2
xdata <- updateObject(xdata)
register(SerialParam())
xdata <- fillChromPeaks(xdata, param = fillChromPeaksParam)

if (exists("intval")) {
    getPeaklistW4M(xdata, intval, convertRTMinute, numDigitsMZ, numDigitsRT, naTOzero, "variableMetadata.tsv", "dataMatrix.tsv")
}

cat("\n\n")

# ----- EXPORT -----

cat("\tXCMSnExp OBJECT INFO\n")
print(xdata)
cat("\n\n")

cat("\txcmsSet OBJECT INFO\n")
# Get the legacy xcmsSet object
xset <- getxcmsSetObject(xdata)
print(xset)
cat("\n\n")

#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xdata", "zipfile", "singlefile", "md5sumList", "sampleNamesList") #, "chromTIC", "chromBPI", "chromTIC_adjusted", "chromBPI_adjusted")
save(list = objects2save[objects2save %in% ls()], file = "fillpeaks.RData")

cat("\n\n")


cat("\tDONE\n")
