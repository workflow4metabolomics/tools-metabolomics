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
options(bitmapType = "cairo")

samplebyclass <- 2
if (!is.null(args$samplebyclass)) {
    samplebyclass <- args$samplebyclass
    args$samplebyclass <- NULL
}

# ----- INFILE PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, NULL, args)
singlefile <- rawFilePath$singlefile
directory <- retrieveRawfileInWD(singlefile, NULL)

# Check some character issues
checkXmlStructure(directory)

cat("\n\n")




# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


xset <- ipo4xcmsSet(directory, "IPO_parameters4xcmsSet.tsv", args, samplebyclass)



cat("\n\n")


# ----- EXPORT -----

cat("\tXSET OBJECT INFO\n")
print(xset)
# delete the parameters to avoid the passage to the next tool in .RData image


# saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xset", "singlefile")
save(list = objects2save[objects2save %in% ls()], file = "ipo4xcmsSet.RData")

cat("\n\n")


cat("\tDONE\n")
