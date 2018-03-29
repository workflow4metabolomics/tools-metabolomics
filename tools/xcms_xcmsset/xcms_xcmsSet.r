#!/usr/bin/env Rscript

# ----- LOG FILE -----
log_file <- file("log.txt", open="wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname){ argv <- commandArgs(trailingOnly=FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }
source_local("lib.r")

pkgs <- c("xcms","batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names=F, quote=F, sep='\t')

cat("\n\n")


# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

#saving the commun parameters
BPPARAM <- MulticoreParam(1)
if (!is.null(args$BPPARAM)){
    BPPARAM <- MulticoreParam(args$BPPARAM); args$BPPARAM <- NULL
}
register(BPPARAM)

#saving the specific parameters
if (!is.null(args$filterAcquisitionNum)){
    filterAcquisitionNumParam <- args$filterAcquisitionNum; args$filterAcquisitionNum <- NULL
}
if (!is.null(args$filterRt)){
    filterRtParam <- args$filterRt; args$filterRt <- NULL
}
if (!is.null(args$filterMz)){
    filterMzParam <- args$filterMz; args$filterMz <- NULL
}

method <- args$method; args$method <- NULL

if (!is.null(args$roiList)){
    cat("\t\troiList provided\n")
    args$roiList <- list(getDataFrameFromFile(args$roiList))
    print(args$roiList)
}

cat("\n\n")

# ----- INFILE PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
load(args$image); args$image=NULL
if (!exists("raw_data")) stop("\n\nERROR: The RData doesn't contain any object called 'raw_data' which is provided by the tool: MSnbase readMSData")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
args <- rawFilePath$args
directory <- retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)


cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tCOMPUTE\n")

## Get the full path to the files
files <- getMSFiles(directory)

cat("\t\t\tApply filter[s] (if asked)\n")
if (exists("filterAcquisitionNumParam")) {
    raw_data <- filterAcquisitionNum(raw_data, filterAcquisitionNumParam[1]:filterAcquisitionNumParam[2])
}
if (exists("filterRtParam")) {
    raw_data <- filterRt(raw_data, filterRtParam)
}
if (exists("filterMzParam")) {
    raw_data <- filterMz(raw_data, filterMzParam)
}

cat("\t\t\tChromatographic peak detection\n")
findChromPeaksParam <- do.call(paste0(method,"Param"), args)
print(findChromPeaksParam)
xdata <- findChromPeaks(raw_data, param=findChromPeaksParam)

# Check if there are no peaks
if (nrow(chromPeaks(xdata)) == 0) stop("No peaks were detected. You should review your settings")

# Transform the files absolute pathways into relative pathways
xdata@processingData@files <- sub(paste(getwd(), "/", sep="") , "", xdata@processingData@files)

# Create a sampleMetada file
sampleNamesList <- getSampleMetadata(xdata=xdata, sampleMetadataOutput="sampleMetadata.tsv")

cat("\n\n")

# ----- EXPORT -----

cat("\tXCMSnExp OBJECT INFO\n")
print(xdata)
cat("\n\n")

cat("\txcmsSet OBJECT INFO\n")
# Get the legacy xcmsSet object
xset <- getxcmsSetObject(xdata)
print(xset)
cat("\n\n")

#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xdata", "zipfile", "singlefile", "md5sumList", "sampleNamesList")
save(list=objects2save[objects2save %in% ls()], file="xcmsSet.RData")


cat("\tDONE\n")
