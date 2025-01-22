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

# Function to convert "NULL" strings to actual NULL values
convertNullString <- function(x) {
    if (x == "NULL") {
        return(NULL)
    }
    return(x)
}

# Function to convert string to numeric lists
convert_psg_list <- function(x) {
    # Check if x is NULL or has zero length before further processing
    if (is.null(x) || length(x) == 0) {
        return(NULL)
    }

    # Force conversion to character
    x <- as.character(x)

    if (grepl("^[0-9]+$", x)) {
        # If the string represents a single numeric value
        return(as.numeric(x))
    } else {
        # Convert string representation of a list to a numeric vector
        # Use a regular expression to split by common separators
        return(as.numeric(unlist(strsplit(x, "[,;\\s]+"))))
    }
}

for (arg in names(args)) {
    args[[arg]] <- convertNullString(args[[arg]])
}

# Convert only the 'psg_list' element in args
args$psg_list <- convert_psg_list(args$psg_list)

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

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[format(Sys.time(), "%y%m%d-%H:%M:%S_groupCorr")]] <- args

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

# Ensure the xa object is loaded
if (!exists("xa")) {
    stop("The xa object was not found in the RData file.")
}

print("xa object loaded:")
print(xa)

# Apply the groupCorr function to the xsAnnotate object
print("Calling groupCorr function:")
xa <- groupCorr(xa, cor_eic_th = args$cor_eic_th, pval = args$pval, graphMethod = args$graphMethod, calcIso = args$calcIso, calcCiS = args$calcCiS, calcCaS = args$calcCaS, psg_list = args$psg_list, cor_exp_th = args$cor_exp_th, intval = args$intval)

print("Result of groupCorr function:")
print(xa)

# Extract the list of annotated peaks
peakList <- getPeaklist(xa, intval = args$intval)

if (length(phenoData@data$sample_name) == 1) {
    peakList$name <- make.unique(paste0("M", round(peakList[, "mz"], 0), "T", round(peakList[, "rt"], 0)), "_")
    variableMetadata <- peakList[, c("name", setdiff(names(peakList), "name"))]
    variableMetadata <- formatIonIdentifiers(variableMetadata, numDigitsRT = args$numDigitsRT, numDigitsMZ = args$numDigitsMZ)
} else {
    names_default <- groupnames(xa@xcmsSet, mzdec = 0, rtdec = 0) # Names without decimals
    names_custom <- groupnames(xa@xcmsSet, mzdec = args$numDigitsMZ, rtdec = args$numDigitsRT) # Names with "x" decimals

    variableMetadata <- data.frame(
        name = names_default,
        name_custom = names_custom,
        stringsAsFactors = FALSE
    )
    variableMetadata <- cbind(variableMetadata, peakList[, !(make.names(colnames(peakList)) %in% c(make.names(sampnames(xa@xcmsSet))))])
}

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
objects2save <- c("xa", "variableMetadata", "listOFlistArguments", "zipfile", "singlefile", "RTinMinute", "phenoData")
save(list = objects2save[objects2save %in% ls()], file = output_file_RData)

cat("Output files generated:", output_file_tsv, "and", output_file_RData, "\n")
