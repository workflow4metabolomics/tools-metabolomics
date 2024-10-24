#!/usr/bin/env Rscript

# Load necessary libraries
library(CAMERA)
library(xcms)

# Retrieve command line arguments
args_vec <- commandArgs(trailingOnly = TRUE)

print("Arguments retrieved from the command line:")
print(args_vec)

# Convert arguments to appropriate types
args <- list(
  image = args_vec[1], # RData input file
  sigma = as.numeric(args_vec[2]), # Convert sigma to numeric
  perfwhm = as.numeric(args_vec[3]), # Convert perfwhm to numeric
  intval = args_vec[4], # intval parameter ("into", "maxo", "intb")
  convertRTMinute = as.logical(args_vec[5]), # TRUE - FALSE
  numDigitsMZ = as.numeric(args_vec[6]), # Number of digits MZ
  numDigitsRT = as.numeric(args_vec[7]), # Number of digits RT
  singlefile_galaxyPath = args_vec[9], # @COMMAND_FILE_LOAD@
  singlefile_sampleName = args_vec[11] # @COMMAND_FILE_LOAD@
)

print("Converted arguments:")
print(args)

# Check if the image file exists
if (!file.exists(args$image)) {
  stop("The RData file does not exist: ", args$image)
}

# Load the RData file (it should contain the xset object, typically an xcmsSet or XCMSnExp)
load(args$image)
args$image <- NULL

# Function to get an xcmsSet or XCMSnExp object
getxcmsSetObject <- function(xobject) {
  # XCMS 1.x
  if (class(xobject) == "xcmsSet") {
    return(xobject)
  }
  # XCMS 3.x
  if (class(xobject) == "XCMSnExp") {
    suppressWarnings(xset <- as(xobject, "xcmsSet"))
    if (is.null(xset@phenoData$sample_group)) {
      sampclass(xset) <- "."
    } else {
      sampclass(xset) <- xset@phenoData$sample_group
    }
    if (!is.null(xset@phenoData$sample_name)) {
      rownames(xset@phenoData) <- xset@phenoData$sample_name
    }
    return(xset)
  }
}

# If the xdata object exists, convert it to xcmsSet
if (exists("xdata")) {
  xset <- getxcmsSetObject(xdata)
}

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

if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS 2.*")

# Verification of a group step before doing the fillpeaks job.
if (!hasFeatures(xdata)) stop("You must always do a group step after a retcor. Otherwise it won't work for the groupFWHM step")

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[format(Sys.time(), "%y%m%d-%H:%M:%S_groupFWHM")]] <- args

# Retrieve raw files
if (!exists("zipfile")) zipfile <- NULL
if (!exists("singlefile")) singlefile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
args <- rawFilePath$args

# Function to retrieve the raw file in the working directory
retrieveRawfileInTheWorkingDir <- function(singlefile, zipfile) {
  if (!is.null(singlefile) && (length(singlefile) > 0)) {
    for (singlefile_sampleName in names(singlefile)) {
      singlefile_galaxyPath <- singlefile[[singlefile_sampleName]]
      if (!file.exists(singlefile_galaxyPath)) {
        error_message <- paste("Cannot access the sample:", singlefile_sampleName, "located:", singlefile_galaxyPath)
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

# Retrieve files
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

# Convert the xset object to xsAnnotate using CAMERA
cat("Converting xset object to xsAnnotate...\n")
xsa <- xsAnnotate(xset)

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

# Save the peak list to a TSV file
output_file_tsv <- "variableMetadata.tsv"
write.table(variableMetadata, file = output_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# Save the xsAnnotate object
output_file_RData <- "camera_fwhm.RData"
objects2save <- c("xa", "variableMetadata", "listOFlistArguments", "zipfile", "singlefile", "RTinMinute")
save(list = objects2save[objects2save %in% ls()], file = output_file_RData)

cat("Output files generated:", output_file_tsv, "and", output_file_RData, "\n")
