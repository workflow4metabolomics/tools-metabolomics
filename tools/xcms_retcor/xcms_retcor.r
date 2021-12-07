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

pkgs <- c("xcms", "batch", "RColorBrewer")
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
method <- args$method

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
load(args$image); args$image <- NULL
if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS 2.*")

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

cat("\t\t\tAlignment/Retention Time correction\n")
# clear the arguement list to remove unexpected key/value as singlefile_galaxyPath or method ...
args <- args[names(args) %in% slotNames(do.call(paste0(method, "Param"), list()))]

adjustRtimeParam <- do.call(paste0(method, "Param"), args)
print(adjustRtimeParam)

if (hasAdjustedRtime(xdata)) {
  cat("WARNING: a retention time ajustment had already been applied to your data.\nThe function applyAdjustedRtime was processed to cumulate the ajustment")
  cat("Replace raw retention times with adjusted retention times.\n")
  xdata <- applyAdjustedRtime(xdata)
}
xdata <- adjustRtime(xdata, param = adjustRtimeParam)

cat("\n\n")


# -- TIC --
cat("\t\tDRAW GRAPHICS\n")
getPlotAdjustedRtime(xdata)

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
save(list = objects2save[objects2save %in% ls()], file = "retcor.RData")

cat("\n\n")


cat("\tDONE\n")
