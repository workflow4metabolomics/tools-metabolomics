#!/usr/bin/env Rscript

# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("CAMERA", "multtest", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");

# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")

args <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = F, quote = F, sep = "\t")

cat("\n\n");


# ----- PROCESSING INFILE -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
load(args$image); args$image <- NULL

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tARGUMENTS PROCESSING INFO\n")

# Save arguments to generate a report
if (!exists("listOFargs")) listOFargs <- list()
listOFargs[[format(Sys.time(), "%y%m%d-%H:%M:%S_annotatediff")]] <- args

# We unzip automatically the chromatograms from the zip files.
if (!exists("zipfile")) zipfile <- NULL
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
args <- rawFilePath$args
directory <- retrieveRawfileInTheWorkingDir(singlefile, zipfile)


# Because so far CAMERA isn't compatible with the new XCMSnExp object
if (exists("xdata")) {
    xset <- getxcmsSetObject(xdata)
}

cat("\n\n")


# ----- PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

results_list <- annotatediff(xset = xset, args = args, variableMetadataOutput = "variableMetadata.tsv")
xa <- results_list$xa
diffrep <- results_list$diffrep
variableMetadata <- results_list$variableMetadata

cat("\n\n")

# ----- EXPORT -----

cat("\tXSET OBJECT INFO\n")
print(xa)
cat("\n\n")

#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xa", "variableMetadata", "diffrep", "cAnnot", "listOFargs", "zipfile", "singlefile")
save(list = objects2save[objects2save %in% ls()], file = "annotatediff.RData")

cat("\n\n")

cat("\tDONE\n")
