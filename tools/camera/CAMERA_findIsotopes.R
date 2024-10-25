#!/usr/bin/env Rscript

# Load necessary libraries
library(CAMERA)
library(xcms)

source_local("lib.r")

# Retrieve command-line arguments
args <- W4MRUtils::parse_args(args = commandArgs())

print("Arguments retrieved from command line:")
print(args)

print("Argument types:")
print(sapply(args, class))

# Argument validation
if (!file.exists(args$image)) {
  stop("The provided RData file does not exist: ", args$image)
}

# Load the RData file
load(args$image)
args$image <- NULL

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[format(Sys.time(), "%y%m%d-%H:%M:%S_findIsotopes")]] <- args

# Retrieve raw files
if (!exists("zipfile")) zipfile <- NULL
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
args <- rawFilePath$args

# Retrieve files
directory <- retrieveRawfileInTheWorkingDir(singlefile, zipfile)

# Verify if the xa object is loaded
if (!exists("xa")) {
  stop("The xa object was not found in the RData file.")
}

print("xa object loaded:")
print(xa)

print("calcIsotopeMatrix")
calcIsotopeMatrix <- function(maxiso = 4) {
  if (!is.numeric(maxiso)) {
    stop("Parameter maxiso is not numeric!\n")
  } else if (maxiso < 1 || maxiso > 8) {
    stop(paste(
      "Parameter maxiso must be between 1 and 8. ",
      "Otherwise, use your own IsotopeMatrix.\n"
    ), sep = "")
  }

  isotopeMatrix <- matrix(NA, 8, 4)
  colnames(isotopeMatrix) <- c("mzmin", "mzmax", "intmin", "intmax")

  isotopeMatrix[1, ] <- c(1.000, 1.0040, 1.0, 150)
  isotopeMatrix[2, ] <- c(0.997, 1.0040, 0.01, 200)
  isotopeMatrix[3, ] <- c(1.000, 1.0040, 0.001, 200)
  isotopeMatrix[4, ] <- c(1.000, 1.0040, 0.0001, 200)
  isotopeMatrix[5, ] <- c(1.000, 1.0040, 0.00001, 200)
  isotopeMatrix[6, ] <- c(1.000, 1.0040, 0.000001, 200)
  isotopeMatrix[7, ] <- c(1.000, 1.0040, 0.0000001, 200)
  isotopeMatrix[8, ] <- c(1.000, 1.0040, 0.00000001, 200)

  return(isotopeMatrix[1:maxiso, , drop = FALSE])
}

isotopeMatrix <- calcIsotopeMatrix(args$maxiso)

# Apply the findIsotopes function on the xsAnnotate object
print("Calling findIsotopes function:")
xa <- findIsotopes(xa, maxcharge = args$maxcharge, maxiso = args$maxiso, ppm = args$ppm, mzabs = args$mzabs, intval = args$intval, minfrac = args$minfrac, isotopeMatrix = isotopeMatrix, filter = args$filter)

print("Result of the findIsotopes function:")
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

# Saves the extracted peak list as a TSV file named 'variableMetadata.tsv'
output_file_tsv <- "variableMetadata.tsv"
write.table(variableMetadata, file = output_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# Save the updated xsAnnotate object
output_file_RData <- "camera_findIsotopes.RData"

objects2save <- c("xa", "variableMetadata", "listOFlistArguments", "zipfile", "singlefile", "RTinMinute")
save(list = objects2save[objects2save %in% ls()], file = output_file_RData)

cat("Output files generated:", output_file_tsv, "and", output_file_RData, "\n")
