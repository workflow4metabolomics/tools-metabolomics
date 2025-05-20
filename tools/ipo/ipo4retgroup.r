#!/usr/bin/env Rscript
# Authors Gildas Le Corguille and Yann Guitton


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

# Import the different functions
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("IPO", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n")
# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) # interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = FALSE, quote = FALSE, sep = "\t")

cat("\n\n")
# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

samplebyclass <- 2
if (!is.null(args$samplebyclass)) {
    samplebyclass <- args$samplebyclass
    args$samplebyclass <- NULL
}

cat("\n\n")

# ----- INFILE PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")
options(bitmapType = "cairo")

# image is an .RData file necessary to use xset variable given by previous tools
load(args$image)
args$image <- NULL

# Because so far CAMERA isn't compatible with the new XCMSnExp object
if (exists("xdata")) {
    xset <- getxcmsSetObject(xdata)
}

if (!exists("xset")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata' which is provided by the tool: MSnbase readMSData")


# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, NULL, args)
singlefile <- rawFilePath$singlefile
print(singlefile)
directory <- retrieveRawfileInWD(singlefile, NULL)


cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

ipo4retgroup(xset, directory, "IPO_parameters4xcmsSet.tsv", args, samplebyclass)

cat("\n\n")


cat("\tDONE\n")
