#!/usr/bin/env Rscript
# Authors:
#   - ABiMS Team
#   - LABERCA - PARC project founding

# ----- LOG FILE -----
log_file <- file("log.txt", open = "wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

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


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE)
write.table(as.matrix(args), col.names = FALSE, quote = FALSE, sep = "\t")

cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

# saving the specific parameters
method <- "ChromPeakArea"

if (!is.null(args$convertRTMinute)) convertRTMinute <- args$convertRTMinute
if (!is.null(args$numDigitsMZ)) numDigitsMZ <- args$numDigitsMZ
if (!is.null(args$numDigitsRT)) numDigitsRT <- args$numDigitsRT
if (!is.null(args$intval)) intval <- args$intval
if (!is.null(args$naTOzero)) naTOzero <- args$naTOzero

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

load(args$image)
if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata' (MsExperiment or XcmsExperiment object)")

# Verification of a group step before doing the fillpeaks job.
if (!hasFeatures(xdata)) stop("You must always do a group step after a retcor. Otherwise it won't work for the fillpeaks step")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
# rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
# zipfile <- rawFilePath$zipfile
# singlefile <- rawFilePath$singlefile

cat("\n\n")

# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

cat("Missing values before peak filling: ")
na_before <- sum(is.na(featureValues(xdata)))
print(na_before)
cat("\n\n")

cat("\t\tCOMPUTE\n")

cat("\t\t\tFilling missing peaks using specified settings\n")

# Median value used
if (!is.null(args$mzmin)) {
    fillChromPeaksParam <- ChromPeakAreaParam(
        mzmin = median,
        mzmax = median,
        rtmin = median,
        rtmax = median,
        minMzWidthPpm = args$minMzWidthPpm
    )
    # Default parameters (quartile)
} else {
    fillChromPeaksParam <- ChromPeakAreaParam(
        minMzWidthPpm = args$minMzWidthPpm
    )
}

cat("fillChromPeaks parameters\n")
print(fillChromPeaksParam)

# XCMS 4.x - MsExperiment/XcmsExperiment objects only.
if (is(xdata, "XCMSnExp") || is(xdata, "OnDiskMSnExp")) {
    stop("\n\nERROR: The RData contains a legacy '", class(xdata)[1], "' object. This function only supports the 'MsExperiment'/'XcmsExperiment' objects produced by xcms >= 4. Please reprocess your data with a recent version of xcms.")
}
# Disable parallel processing to avoid memory known bugs of xcms
register(SerialParam())
# Fill peaks with defined paramaters
xdata <- fillChromPeaks(xdata, param = fillChromPeaksParam)
cat("\n\n")
cat("Missing values after peak filling: ")
na_after <- sum(is.na(featureValues(xdata)))
print(na_after)
cat("\n\n")

if (exists("intval")) {
    getPeaklistW4M(
        xdata,
        intval,
        convertRTMinute,
        numDigitsMZ,
        numDigitsRT,
        naTOzero,
        "variableMetadata.tsv",
        "dataMatrix.tsv"
    )
}

cat("\n\n")

# ----- EXPORT -----

cat("\tXcmsExperiment OBJECT INFO\n")
print(xdata)
cat("\n\n")

# saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xdata", "md5sumList", "sampleNamesList")
save(list = objects2save[objects2save %in% ls()], file = "fillpeaks.RData")

cat("\n\n")

cat("\tDONE\n")
