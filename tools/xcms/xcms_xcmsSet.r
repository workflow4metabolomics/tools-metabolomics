#!/usr/bin/env Rscript
# Authors:
#   - ABiMS Team
#   - LABERCA - PARC project founding

# ----- LOG FILE -----
log_file <- file("log.txt", open = "wt")
sink(log_file)
sink(log_file, type = "output")

# ----- PACKAGE -----
cat("\tSESSION INFO\n")

# Import the different functions
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("xcms", "batch", "Spectra", "MsExperiment")
loadAndDisplayPackages(pkgs)
cat("\n\n")

# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) # interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = FALSE, quote = FALSE, sep = "\t")

cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

# saving the commun parameters
BPPARAM <- MulticoreParam(1)
if (!is.null(args$BPPARAM)) {
    BPPARAM <- MulticoreParam(args$BPPARAM)
}
register(BPPARAM)

# saving the specific parameters
if (!is.null(args$filterAcquisitionNum)) filterAcquisitionNumParam <- args$filterAcquisitionNum
if (!is.null(args$filterRt)) filterRtParam <- args$filterRt
if (!is.null(args$filterMz)) filterMzParam <- args$filterMz
if (!is.null(args$peaklist)) peaklistParam <- args$peaklist

method <- args$method

if (!is.null(args$roiList)) {
    cat("\t\troiList provided\n")
    args$roiList <- list(getDataFrameFromFile(args$roiList))
    print(args$roiList)
}

cat("\n\n")

# ----- INFILE PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

# image is an .RData file necessary to use xset variable given by previous tools
load(args$image)
if (!exists("raw_data")) stop("\n\nERROR: The RData doesn't contain any object called 'raw_data' which is provided by the tool: readMsExperiment")

# Handle infiles
# rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
# zipfile <- rawFilePath$zipfile
# singlefile <- rawFilePath$singlefile

cat("\n\n")

# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

cat("\t\tCOMPUTE\n")

cat("\t\t\tApply filter[s] (if asked)\n")
if (exists("filterAcquisitionNumParam")) {
    raw_data <- filterSpectra(raw_data, filterAcquisitionNum, n = filterAcquisitionNumParam[1]:filterAcquisitionNumParam[2])
}
if (exists("filterRtParam")) {
    raw_data <- filterRt(raw_data, filterRtParam)
}
if (exists("filterMzParam")) {
    raw_data <- filterMz(raw_data, filterMzParam)
}
# Apply this filter only if file contain MS and MSn
if (length(unique(msLevel(spectra(raw_data)))) != 1) {
    raw_data <- filterSpectra(raw_data, filterMsLevel, msLevel. = 1)
}

cat("\t\t\tChromatographic peak detection\n")
# Clear the arguement list to remove unexpected key/value
args <- args[names(args) %in% slotNames(do.call(paste0(method, "Param"), list()))]

findChromPeaksParam <- do.call(paste0(method, "Param"), args)
print(findChromPeaksParam)
xdata <- findChromPeaks(raw_data, param = findChromPeaksParam)

# Check if there are no peaks
if (nrow(chromPeaks(xdata)) == 0) stop("No peaks were detected. You should review your settings")

# Create a sampleMetada file
sampleNamesList <- getSampleMetadata(xdata = xdata, sampleMetadataOutput = "sampleMetadata.tsv")

# Create a chromPeaks table if required
if (exists("peaklistParam")) {
    if (peaklistParam) {
        cat("\nCreating the chromatographic peaks' table...\n")
        write.table(chromPeaks(xdata), file = "chromPeak_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

cat("\n\n")

# ----- EXPORT -----

cat("\tMsExperiment OBJECT INFO\n")
print(xdata)
cat("\n\n")
# saving in .Rdata file
objects2save <- c("xdata", "md5sumList", "sampleNamesList")
save(list = objects2save[objects2save %in% ls()], file = "XcmsExpDatasets.RData")

cat("\tDONE\n")
