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
source_local("xcms4_lib.r")

pkgs <- c("MsExperiment", "xcms")
loadAndDisplayPackages(pkgs)
cat("\n\n")

# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) # interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = FALSE, quote = FALSE, sep = "\t")

cat("\n\n")


# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")


cat("\n\n")

# ----- INFILE PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL # no zipfile allowed in 1st dev version
rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
files <- rawFilePath$files
# print(files)

md5sumList <- list("origin" = getMd5sum(files))

cat("\n\n")

# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tCOMPUTE\n")

cat("\t\t\tCreate a metaodata data.frame\n")
s_groups <- sapply(files, function(x) tail(unlist(strsplit(dirname(x), "/")), n = 1))
s_name <- tools::file_path_sans_ext(basename(files))
pd <- data.frame(sample_name = s_name, sample_group = s_groups, stringsAsFactors = FALSE)
# print(pd)

cat("\t\t\tRead Ms Experiment\n")
raw_data <- readMsExperiment(spectraFiles = files, sampleData = pd, backend = Spectra::MsBackendMemory())
print(raw_data)

# Transform the files absolute pathways into relative pathways
raw_data@sampleData@listData$spectraOrigin <- sub(paste(getwd(), "/", sep = ""), "", raw_data@sampleData@listData$spectraOrigin)
raw_data@spectra@backend@spectraData$dataOrigin <- sub(paste(getwd(), "/", sep = ""), "", raw_data@sampleData@listData$spectraOrigin)

# # Create a sampleMetada file
# sampleNamesList <- getSampleMetadata(xdata = raw_data, sampleMetadataOutput = "sampleMetadata.tsv")
# print(sampleNamesList)

cat("\n\n")

# ----- EXPORT -----

cat("\tMSExperiment OBJECT INFO\n")
print(raw_data)
cat("\t\tmetadaData\n")
print(raw_data@sampleData)
cat("\n\n")

cat("\nSAVING\n")
# saving R data in .Rdata file to save the variables used in the present tool
# objects2save <- c("raw_data", "zipfile", "singlefile", "md5sumList", "sampleNamesList") # , "chromTIC", "chromBPI")
objects2save <- c("raw_data", "zipfile", "singlefile", "md5sumList")
save(list = objects2save[objects2save %in% ls()], file = "readmsexp.RData")


cat("\tDONE\n")
