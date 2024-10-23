#!/usr/bin/env Rscript

# Load necessary libraries
library(CAMERA)
library(xcms)

# Retrieve command-line arguments
args_vec <- commandArgs(trailingOnly = TRUE)

print("Arguments retrieved from the command line:")
print(args_vec)

args <- list(
  image = args_vec[1], # the xsAnnotate object
  ppm = as.numeric(args_vec[2]), # ppm error for the search
  mzabs = as.numeric(args_vec[3]), # allowed variance for the search
  multiplier = as.numeric(args_vec[4]), # highest number(n) of allowed clusterion [nM+ion]
  polarity = args_vec[5], # Which polarity mode was used for measuring of the MS sample
  rules = args_vec[6], # custom ruleset or NULL for default ruleset
  max_peaks = as.numeric(args_vec[7]), # If run in parallel mode, this number defines how many peaks will be calculated in each thread
  psg_list = args_vec[8], # Vector of pseudospectra indices; correlation analysis will only be done for those groups
  intval = args_vec[9], # choose intensity values. Allowed values are "into", "maxo", "intb" (string)
  convertRTMinute = as.logical(args_vec[10]), # TRUE - FALSE
  numDigitsMZ = as.numeric(args_vec[11]), # Number of digits MZ
  numDigitsRT = as.numeric(args_vec[12]), # Number of digits RT
  singlefile_galaxyPath = args_vec[14], # @COMMAND_FILE_LOAD@
  singlefile_sampleName = args_vec[16] # @COMMAND_FILE_LOAD@
)

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

print("Converted arguments:")
print(args)

print("Argument types:")
print(sapply(args, class))

# Verify the arguments
if (!file.exists(args$image)) {
  stop("The provided RData file does not exist: ", args$image)
}

# Load the RData file
load(args$image)
args$image <- NULL

# Function to retrieve the raw file from the arguments
getRawfilePathFromArguments <- function(singlefile, zipfile, args) {
  if (!is.null(args$zipfile)) zipfile <- args$zipfile
  if (!is.null(args$zipfilePositive)) zipfile <- args$zipfilePositive
  if (!is.null(args$zipfileNegative)) zipfile <- args$zipfileNegative

  if (!is.null(args$singlefile_galaxyPath)) {
    singlefile_galaxyPaths <- args$singlefile_galaxyPath
    singlefile_sampleNames <- args$singlefile_sampleName
  }
  if (!is.null(args$singlefile_galaxyPathPositive)) {
    singlefile_galaxyPaths <- args$singlefile_galaxyPathPositive
    singlefile_sampleNames <- args$singlefile_sampleNamePositive
  }
  if (!is.null(args$singlefile_galaxyPathNegative)) {
    singlefile_galaxyPaths <- args$singlefile_galaxyPathNegative
    singlefile_sampleNames <- args$singlefile_sampleNameNegative
  }
  if (exists("singlefile_galaxyPaths")) {
    singlefile_galaxyPaths <- unlist(strsplit(singlefile_galaxyPaths, ","))
    singlefile_sampleNames <- unlist(strsplit(singlefile_sampleNames, ","))

    singlefile <- NULL
    for (singlefile_galaxyPath_i in seq_len(length(singlefile_galaxyPaths))) {
      singlefile_galaxyPath <- singlefile_galaxyPaths[singlefile_galaxyPath_i]
      singlefile_sampleName <- singlefile_sampleNames[singlefile_galaxyPath_i]
      singlefile[[singlefile_sampleName]] <- singlefile_galaxyPath
    }
  }
  for (argument in c(
    "zipfile", "zipfilePositive", "zipfileNegative",
    "singlefile_galaxyPath", "singlefile_sampleName",
    "singlefile_galaxyPathPositive", "singlefile_sampleNamePositive",
    "singlefile_galaxyPathNegative", "singlefile_sampleNameNegative"
  )) {
    args[[argument]] <- NULL
  }
  return(list(zipfile = zipfile, singlefile = singlefile, args = args))
}

# Save arguments for report generation
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[format(Sys.time(), "%y%m%d-%H:%M:%S_findAdducts")]] <- args

# Retrieve raw files
if (!exists("zipfile")) zipfile <- NULL
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
args <- rawFilePath$args

# Function to retrieve raw files in the working directory
retrieveRawfileInTheWorkingDir <- function(singlefile, zipfile) {
  if (!is.null(singlefile) && (length(singlefile) > 0)) {
    for (singlefile_sampleName in names(singlefile)) {
      singlefile_galaxyPath <- singlefile[[singlefile_sampleName]]
      if (!file.exists(singlefile_galaxyPath)) {
        error_message <- paste("Cannot access the sample:", singlefile_sampleName, "located at:", singlefile_galaxyPath)
        print(error_message)
        stop(error_message)
      }
      file.symlink(singlefile_galaxyPath, singlefile_sampleName)
    }
    directory <- "."
  }
  if (!is.null(zipfile) && (zipfile != "")) {
    if (!file.exists(zipfile)) {
      error_message <- paste("Cannot access the Zip file:", zipfile)
      print(error_message)
      stop(error_message)
    }
    suppressWarnings(unzip(zipfile, unzip = "unzip"))
    filesInZip <- unzip(zipfile, list = TRUE)
    directories <- unique(unlist(lapply(strsplit(filesInZip$Name, "/"), function(x) x[1])))
    directories <- directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
    directory <- "."
    if (length(directories) == 1) directory <- directories
    cat("files_root_directory\t", directory, "\n")
  }
  return(directory)
}

# Retrieve the files
directory <- retrieveRawfileInTheWorkingDir(singlefile, zipfile)

# @author G. Le Corguille
# This function convert if it is required the Retention Time in minutes
RTSecondToMinute <- function(variableMetadata, convertRTMinute) {
  if (convertRTMinute) {
    # converting the retention times (seconds) into minutes
    print("converting the retention times into minutes in the variableMetadata")
    variableMetadata[, "rt"] <- variableMetadata[, "rt"] / 60
    variableMetadata[, "rtmin"] <- variableMetadata[, "rtmin"] / 60
    variableMetadata[, "rtmax"] <- variableMetadata[, "rtmax"] / 60
  }
  return(variableMetadata)
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
