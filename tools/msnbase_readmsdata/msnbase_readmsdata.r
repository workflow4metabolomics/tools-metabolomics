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

pkgs <- c("MSnbase", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = F, quote = F, sep = "\t")

cat("\n\n")


# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")


cat("\n\n")

# ----- INFILE PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
files <- rawFilePath$files

md5sumList <- list("origin" = getMd5sum(files))

cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tCOMPUTE\n")

cat("\t\t\tCreate a phenodata data.frame\n")
s_groups <- sapply(files, function(x) tail(unlist(strsplit(dirname(x), "/")), n = 1))
s_name <- tools::file_path_sans_ext(basename(files))
pd <- data.frame(sample_name = s_name, sample_group = s_groups, stringsAsFactors = FALSE)
print(pd)

cat("\t\t\tLoad Raw Data\n")
raw_data <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")

# Transform the files absolute pathways into relative pathways
raw_data@processingData@files <- sub(paste(getwd(), "/", sep = ""), "", raw_data@processingData@files)

# Create a sampleMetada file
sampleNamesList <- getSampleMetadata(xdata = raw_data, sampleMetadataOutput = "sampleMetadata.tsv")

cat("\n\n")

# ----- EXPORT -----

cat("\tMSnExp OBJECT INFO\n")
print(raw_data)
cat("\t\tphenoData\n")
print(raw_data@phenoData@data)
cat("\n\n")

#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("raw_data", "zipfile", "singlefile", "md5sumList", "sampleNamesList") #, "chromTIC", "chromBPI")
save(list = objects2save[objects2save %in% ls()], file = "readmsdata.RData")


cat("\tDONE\n")
