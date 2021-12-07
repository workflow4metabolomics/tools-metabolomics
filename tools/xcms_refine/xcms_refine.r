#!/usr/bin/env Rscript

# ----- LOG FILE -----
log_file <- file("log.txt", open  =  "wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("xcms", "batch", "RColorBrewer")
loadAndDisplayPackages(pkgs)
cat("\n\n");

# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
# interpretation of arguments given in command line as an R list of objects
args <- parseCommandArgs(evaluate = FALSE)
write.table(as.matrix(args), col.names = F, quote = F, sep = "\t")

cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

#saving the specific parameters
args_method  <- args$method
args_image   <- args$image
args_msLevel <- args$msLevel
param_args <- list()

if (args_method == "CleanPeaks") {
  param_args$maxPeakwidth <- args$maxPeakwidth
} else if (args_method == "FilterIntensity") {
  param_args$threshold <- args$threshold
  param_args$value <- args$value
  param_args$nValues <- args$nValues
} else if (args_method == "MergeNeighboringPeaks") {
  param_args$expandRt <- args$expandRt
  param_args$expandMz <- args$expandMz
  param_args$ppm      <- args$ppm
  param_args$minProp  <- args$minProp
}

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
load(args_image)
if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. Such RData as this might have been created by an old version of XMCS 2.*")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile

cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tPREPARE PARAMETERS\n\n")

if (args_method == "CleanPeaks") {
  refineChromPeaksParam <- CleanPeaksParam(maxPeakwidth = param_args$maxPeakwidth)
} else if (args_method == "FilterIntensity") {
  refineChromPeaksParam <- FilterIntensityParam(
    threshold = param_args$threshold,
    nValues = param_args$nValues,
    value = param_args$value
  )
} else if (args_method == "MergeNeighboringPeaks") {
  refineChromPeaksParam <- MergeNeighboringPeaksParam(
    expandRt = param_args$expandRt,
    expandMz = param_args$expandMz,
    ppm = param_args$ppm,
    minProp = param_args$minProp
  )
}

cat(str(refineChromPeaksParam))

cat("\n\n\t\tCOMPUTE\n")

xdata <- updateObject(xdata)

xdata <- refineChromPeaks(xdata, param = refineChromPeaksParam)

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
objects2save <- c("xdata", "zipfile", "singlefile", "md5sumList", "sampleNamesList")
save(list = objects2save[objects2save %in% ls()], file = "xcmsSet.RData")

cat("\n\n")


cat("\tDONE\n")
