#!/usr/bin/env Rscript

# ----- PACKAGE -----
cat("\tSESSION INFO\n")

# Import the different functions
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("CAMERA", "xcms", "multtest", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n")
# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")

args <- parseCommandArgs(evaluate = FALSE) # interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = FALSE, quote = FALSE, sep = "\t")

cat("\n\n")

print("Arguments retrieved from the command line:")
print(args)

print("Argument types:")
print(sapply(args, class))

# Check if the image file exists
if (!file.exists(args$image)) {
  stop("The RData file does not exist: ", args$image)
}

# ----- PROCESSING INFILE -----

# Load the RData file
load(args$image)
args$image <- NULL

# Save arguments for report generation
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[format(Sys.time(), "%y%m%d-%H:%M:%S_findAdducts")]] <- args

# We unzip automatically the chromatograms from the zip files.
if (!exists("zipfile")) zipfile <- NULL
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
args <- rawFilePath$args
directory <- retrieveRawfileInTheWorkingDir(singlefile, zipfile)

# Check if the 'rules' argument in 'args' is NULL
if (is.null(args$rules)) {
  # If 'args$rules' is NULL, set 'rulset' to NULL
  args$rulset <- NULL
} else {
  # Try to read the rules file with different delimiters
  delimiters <- c(";", "\t", ",") # List of possible delimiters
  success <- FALSE # Flag to check if reading was successful

  for (sep in delimiters) {
    # Attempt to read the rules file with the current separator
    args$rulset <- read.table(args$rules, header = TRUE, sep = sep)

    # Check if the number of columns is at least 4
    if (ncol(args$rulset) >= 4) {
      success <- TRUE # Mark success if the format is correct
      break # Exit the loop if the file was read successfully
    }
  }

  # If reading the rules file failed for all delimiters
  if (!success) {
    # Display an error message if the file is not well formatted
    error_message <- "The rules file appears to be improperly formatted. Accepted column separators are ;, tab, and ,."
    print(error_message)
    stop(error_message) # Stop execution with an error
  }
}

# Verify that the object xa is loaded
if (!exists("xa")) {
  stop("The object xa was not found in the RData file.")
}

print("Loaded xa object:")
print(xa)

# Apply the findAdducts function on the xsAnnotate object
print("Calling findAdducts function:")
xa <- findAdducts(xa, ppm = args$ppm, mzabs = args$mzabs, multiplier = args$multiplier, polarity = args$polarity, rules = args$rulset, max_peaks = args$max_peaks, psg_list = args$psg_list, intval = args$intval)

print("Result of findAdducts function:")
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

if (args$convertRTMinute && RTinMinute == FALSE) {
  RTinMinute <- TRUE
  variableMetadata <- RTSecondToMinute(variableMetadata = variableMetadata, convertRTMinute = args$convertRTMinute)
}

# Save the extracted peak list as a TSV file named 'variableMetadata.tsv'
output_file_tsv <- "variableMetadata.tsv"
write.table(variableMetadata, file = output_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# Save the updated xsAnnotate object
output_file_RData <- "camera_findAdducts.RData"
objects2save <- c("xa", "variableMetadata", "listOFlistArguments", "zipfile", "singlefile", "RTinMinute")
save(list = objects2save[objects2save %in% ls()], file = output_file_RData)

cat("Output files generated:", output_file_tsv, "and", output_file_RData, "\n")
