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

# Function to convert "NA" strings to actual NA values and string lists to numeric lists
convertStringToNumeric <- function(x) {
  # Force conversion to character
  x <- as.character(x)

  if (x == "NA") {
    return(NA)
  } else if (grepl("^[0-9]+$", x)) {
    # If the string represents a single numeric value
    return(as.numeric(x))
  } else {
    # Convert string representation of a list to a numeric vector
    # Use a regular expression to split by common separators
    return(as.numeric(unlist(strsplit(x, "[,;\\s]+"))))
  }
}

# Convert only the 'sample' element in args
args$sample <- convertStringToNumeric(args$sample)

print("Argument types:")
print(sapply(args, class))

# Check if the image file exists
if (!file.exists(args$image)) {
  stop("The RData file does not exist: ", args$image)
}

# ----- PROCESSING INFILE -----

# Load the RData file (it should contain the xset object, typically an xcmsSet or XCMSnExp)
load(args$image)
args$image <- NULL

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[format(Sys.time(), "%y%m%d-%H:%M:%S_groupFWHM")]] <- args

# We unzip automatically the chromatograms from the zip files.
if (!exists("zipfile")) zipfile <- NULL
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
args <- rawFilePath$args

print(paste("singlefile :", singlefile))
if (!is.null(singlefile)) {
  directory <- retrieveRawfileInTheWorkingDir(singlefile, zipfile)
}

# If the xdata object exists, convert it to xcmsSet
if (exists("xdata")) {
  phenoData <- xdata@phenoData
  xset <- getxcmsSetObject(xdata)
}

if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS 2.*")

# Verification of a group step before doing the fillpeaks job.
if (!hasFeatures(xdata)) stop("You must always do a group step after a retcor. Otherwise it won't work for the groupFWHM step")

# Convert the xset object to xsAnnotate using CAMERA
cat("Converting xset object to xsAnnotate...\n")
xsa <- xsAnnotate(xset, sample = args$sample, nSlaves = as.numeric(args$nSlaves), polarity = args$polarity)

# Apply the groupFWHM function with the parameters
cat("Applying groupFWHM...\n")
xa <- groupFWHM(xsa, sigma = as.numeric(args$sigma), perfwhm = as.numeric(args$perfwhm), intval = args$intval)

# Extract the list of annotated peaks
peakList <- getPeaklist(xa, intval = args$intval)

# Generate group names with and without decimals for mz and rt
names_default <- groupnames(xa@xcmsSet, mzdec = 0, rtdec = 0) # Names without decimals
names_custom <- groupnames(xa@xcmsSet, mzdec = args$numDigitsMZ, rtdec = args$numDigitsRT) # Names with "x" decimals

# Calculate indices of the columns to include from peakList
# Select all columns except the last sample-specific columns
ncols <- length(colnames(peakList))
sample_cols <- length(rownames(xdata@phenoData)) # Number of samples

# Indices for the columns of interest
main_cols <- 1:(ncols - sample_cols - 3) # Main columns before sample columns
tail_cols <- (ncols - 2):ncols # The last 3 columns (adduct, isotope, pcgroup)

# Combine the selected columns from matgrp with the group names
variableMetadata <- data.frame(
  name = names_default,
  name_custom = names_custom,
  stringsAsFactors = FALSE
)

variableMetadata <- cbind(variableMetadata, peakList[, c(main_cols, tail_cols)])

if (!exists("RTinMinute")) RTinMinute <- FALSE

if (args$convertRTMinute && RTinMinute == FALSE) {
  RTinMinute <- TRUE
  variableMetadata <- RTSecondToMinute(variableMetadata = variableMetadata, convertRTMinute = args$convertRTMinute)
}

# Save the peak list to a TSV file
output_file_tsv <- "variableMetadata.tsv"
write.table(variableMetadata, file = output_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# Save the xsAnnotate object
output_file_RData <- "camera_fwhm.RData"
objects2save <- c("xa", "variableMetadata", "listOFlistArguments", "zipfile", "singlefile", "RTinMinute", "phenoData")
save(list = objects2save[objects2save %in% ls()], file = output_file_RData)

cat("Output files generated:", output_file_tsv, "and", output_file_RData, "\n")
