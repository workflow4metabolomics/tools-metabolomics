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

pkgs <- c("xcms", "batch", "RColorBrewer", "plotly", "htmlwidgets", "htmltools", "MsExperiment")
loadAndDisplayPackages(pkgs)
cat("\n\n")


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) # interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = FALSE, quote = FALSE, sep = "\t")

cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

# saving the specific parameters
method <- args$method

if (!is.null(args$convertRTMinute)) convertRTMinute <- args$convertRTMinute
if (!is.null(args$numDigitsMZ)) numDigitsMZ <- args$numDigitsMZ
if (!is.null(args$numDigitsRT)) numDigitsRT <- args$numDigitsRT
if (!is.null(args$intval)) intval <- args$intval
if (!is.null(args$naTOzero)) naTOzero <- args$naTOzero

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

# image is an .RData file necessary to use xset variable given by previous tools
load(args$image)
if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'from XCMS version 4")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
# rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
# zipfile <- rawFilePath$zipfile
# singlefile <- rawFilePath$singlefile

cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tCOMPUTE\n")


cat("\t\t\tPerform the correspondence\n")
args$sampleGroups <- MsExperiment::sampleData(xdata)$sample_group

# clear the argument list
args <- args[names(args) %in% slotNames(do.call(paste0(method, "Param"), list(sampleGroups = c(1, 2))))]

groupChromPeaksParam <- do.call(paste0(method, "Param"), args)
print(groupChromPeaksParam)
xdata <- groupChromPeaks(xdata, param = groupChromPeaksParam)


cat("\t\tDRAW GRAPHICS\n")
if (method == "PeakDensity") {
    getPlotChromPeakDensity(xdata, param = groupChromPeaksParam)
} else {
    cat("\t\t\tNo density plot available for method", method, "\n")
}

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

# cat("\n\n")

# ----- EXPORT -----

cat("\tXcmsExperiment OBJECT INFO\n")
print(xdata)
cat("\n\n")

# saving in .Rdata file
objects2save <- c("xdata", "md5sumList", "sampleNamesList") # , "chromTIC", "chromBPI", "chromTIC_adjusted", "chromBPI_adjusted")
save(list = objects2save[objects2save %in% ls()], file = "group.RData")

cat("\n\n")

cat("\tDONE\n")
