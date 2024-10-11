#!/usr/bin/env Rscript

# Load necessary libraries
library(CAMERA)
library(xcms)

# Retrieve command line arguments
args_vec <- commandArgs(trailingOnly = TRUE)

print("Arguments retrieved from command line:")
print(args_vec)

args <- list(
	     image           = args_vec[1],                  # the xsAnnotate object
	     maxcharge       = as.numeric(args_vec[2]),      # max. number of the isotope charge
	     maxiso          = as.numeric(args_vec[3]),      # max. number of the isotope peaks
	     ppm             = as.numeric(args_vec[4]),      # ppm error for the search
	     mzabs           = as.numeric(args_vec[5]),      # allowed variance for the search
	     intval          = args_vec[6],                  # Intensity parameter ("into", "maxo", "intb")
	     minfrac         = as.numeric(args_vec[7]),      # In case of multiple samples, percentage value of samples, which must contain the correct C12/C13 ratio and are not NA
	     filter          = as.logical(args_vec[8]),      # Should C12/C13 filter be applied?
	     convertRTMinute = as.logical(args_vec[9]),      # TRUE - FALSE
             numDigitsMZ     = as.numeric(args_vec[10]),     # Number of digits MZ
             numDigitsRT     = as.numeric(args_vec[11])      # Number of digits RT
)

print("Converted arguments:")
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

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[format(Sys.time(), "%y%m%d-%H:%M:%S_findIsotopes")]] <- args

# Retrieve raw files
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

# Verify if the xa object is loaded
if (!exists("xa")) {
	stop("The xa object was not found in the RData file.")
}

print("xa object loaded:")
print(xa)

print("calcIsotopeMatrix")
calcIsotopeMatrix <- function(maxiso = 4){

	if(!is.numeric(maxiso)){
		stop("Parameter maxiso is not numeric!\n")
	} else if(maxiso < 1 | maxiso > 8){
		stop(paste("Parameter maxiso must be between 1 and 8. ",
			   "Otherwise, use your own IsotopeMatrix.\n"),sep="")
	}

	isotopeMatrix <- matrix(NA, 8, 4);
	colnames(isotopeMatrix) <- c("mzmin", "mzmax", "intmin", "intmax")

	isotopeMatrix[1, ] <- c(1.000, 1.0040, 1.0, 150)
	isotopeMatrix[2, ] <- c(0.997, 1.0040, 0.01, 200)
	isotopeMatrix[3, ] <- c(1.000, 1.0040, 0.001, 200)
	isotopeMatrix[4, ] <- c(1.000, 1.0040, 0.0001, 200)
	isotopeMatrix[5, ] <- c(1.000, 1.0040, 0.00001, 200)
	isotopeMatrix[6, ] <- c(1.000, 1.0040, 0.000001, 200)
	isotopeMatrix[7, ] <- c(1.000, 1.0040, 0.0000001, 200)
	isotopeMatrix[8, ] <- c(1.000, 1.0040, 0.00000001, 200)

	return(isotopeMatrix[1:maxiso, , drop=FALSE])
}

isotopeMatrix <- calcIsotopeMatrix(args$maxiso)

# Apply the findIsotopes function on the xsAnnotate object
print("Calling findIsotopes function:")
xa <- findIsotopes(xa, maxcharge=args$maxcharge, maxiso=args$maxiso, ppm=args$ppm, mzabs=args$mzabs, intval=args$intval, minfrac=args$minfrac, isotopeMatrix=isotopeMatrix, filter=args$filter)

print("Result of the findIsotopes function:")
print(xa)

# Extract the list of annotated peaks
peakList <- getPeaklist(xa, intval = args$intval)

# Generate group names with and without decimals for mz and rt
names_default <- groupnames(xa@xcmsSet, mzdec = 0, rtdec = 0)                                  # Names without decimals
names_custom <- groupnames(xa@xcmsSet, mzdec = args$numDigitsMZ , rtdec = args$numDigitsRT)    # Names with "x" decimals

# Calculate indices of the columns to include from peakList
# Select all columns except the last sample-specific columns
ncols <- length(colnames(peakList))
sample_cols <- length(xa@sample)   # Number of samples

# Indices for the columns of interest
main_cols <- 1:(ncols - (3 + sample_cols))  # Main columns before sample columns
tail_cols <- (ncols - 2):ncols              # The last 3 columns (adduct, isotope, pcgroup)

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
write.table(variableMetadata, file=output_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# Save the updated xsAnnotate object
output_file_RData <- "camera_findIsotopes.RData"
objects2save <- c("xa", "variableMetadata", "listOFlistArguments", "zipfile", "singlefile", "RTinMinute")
save(list = objects2save[objects2save %in% ls()], file = output_file_RData)

cat("Output files generated:", output_file_tsv, "and", output_file_RData, "\n")

