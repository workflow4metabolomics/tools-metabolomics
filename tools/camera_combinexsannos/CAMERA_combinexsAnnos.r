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
load(args$image_pos)
xaP <- xa

diffrepP <- NULL
if (exists("diffrep")) diffrepP <- diffrep

load(args$image_neg)
xaN <- xa

diffrepN <- NULL
if (exists("diffrep")) diffrepN <- diffrep


cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tARGUMENTS PROCESSING INFO\n")

# Save arguments to generate a report
if (!exists("listOFargs")) listOFargs <- list()
listOFargs[[format(Sys.time(), "%y%m%d-%H:%M:%S_combinexsAnnos")]] <- args

cat("\n\n")


# ----- PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

cAnnot <- combinexsAnnos_function(
    xaP = xaP, xaN = xaN,
    diffrepP = diffrepP, diffrepN = diffrepN,
    pos = args$pos, tol = args$tol, ruleset = args$ruleset, keep_meta = args$keep_meta,
    convertRTMinute = args$convertRTMinute, numDigitsMZ = args$numDigitsMZ, numDigitsRT = args$numDigitsRT,
    variableMetadataOutput = "variableMetadata.tsv"
)

# ----- EXPORT -----

#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xa", "variableMetadata", "diffrep", "cAnnot", "listOFargs", "zipfile", "singlefile")
save(list = objects2save[objects2save %in% ls()], file = "combinexsAnnos.RData")

cat("\n\n")

cat("\tDONE\n")
