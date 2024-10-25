#!/usr/bin/env Rscript

# Load CAMERA library
library(CAMERA)
library(xcms)

source_local("lib.r")

# Retrieve command-line arguments
args <- W4MRUtils::parse_args(args = commandArgs())

print("Command-line arguments retrieved:")
print(args)

print("Argument types:")
print(sapply(args, class))


# Verify if the arguments are correct
if (!file.exists(args$image)) {
  stop("The provided RData file does not exist: ", args$image)
}

# Load the RData file
load(args$image)
args$image <- NULL

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[format(Sys.time(), "%y%m%d-%H:%M:%S_groupCorr")]] <- args

# Retrieve the raw files
if (!exists("zipfile")) zipfile <- NULL
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
args <- rawFilePath$args

# Retrieve the files
directory <- retrieveRawfileInTheWorkingDir(singlefile, zipfile)

# Ensure the xa object is loaded
if (!exists("xa")) {
  stop("The xa object was not found in the RData file.")
}

print("xa object loaded:")
print(xa)

# Apply the groupCorr function to the xsAnnotate object
print("Calling groupCorr function:")
xa <- groupCorr(xa, cor_eic_th = args$cor_eic_th, pval = args$pval, graphMethod = args$graphMethod, calcIso = args$calcIso, calcCiS = args$calcCiS, calcCaS = args$calcCaS, cor_exp_th = args$cor_exp_th, intval = args$intval)

print("Result of groupCorr function:")
print(xa)

# Extract the list of annotated peaks
peakList <- getPeaklist(xa, intval = args$intval)

# Generate group names with and without decimals for mz and rt
names_default <- groupnames(xa@xcmsSet, mzdec = 0, rtdec = 0) # Names without decimals
names_custom <- groupnames(xa@xcmsSet, mzdec = args$numDigitsMZ, rtdec = args$numDigitsRT) # Names with "x" decimals

# Calculate indices of the columns to include from peakList
# Select all columns except the last sample-specific columns
ncols <- length(colnames(peakList))
sample_cols <- length(xa@sample) # Number of samples

# Indices for the columns of interest
main_cols <- 1:(ncols - (6 + sample_cols)) # Main columns before sample columns
tail_cols <- (ncols - 2):ncols # The last 3 columns (adduct, isotope, pcgroup)

# Combine the selected columns from matgrp with the group names
variableMetadata <- cbind(
  name = names_default,
  name_custom = names_custom,
  peakList[, c(main_cols, tail_cols)]
)

if (!exists("RTinMinute")) RTinMinute <- FALSE

if (args$convertRTMinute == TRUE && RTinMinute == FALSE) {
  RTinMinute <- TRUE
  variableMetadata <- RTSecondToMinute(variableMetadata = variableMetadata, convertRTMinute = args$convertRTMinute)
}

# Save the extracted peak list as a TSV file named 'variableMetadata.tsv'
output_file_tsv <- "variableMetadata.tsv"
write.table(variableMetadata, file = output_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# Save the updated xsAnnotate object
output_file_RData <- "camera_groupCorr.RData"
objects2save <- c("xa", "variableMetadata", "listOFlistArguments", "zipfile", "singlefile", "RTinMinute")
save(list = objects2save[objects2save %in% ls()], file = output_file_RData)

cat("Output files generated:", output_file_tsv, "and", output_file_RData, "\n")
