#!/usr/bin/env Rscript


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("xcms", "batch", "RColorBrewer")
loadAndDisplayPackages(pkgs)
cat("\n\n");


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = F, quote = F, sep = "\t")

cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

mergeXDataReturn <- mergeXData(args)
xdata <- mergeXDataReturn$xdata
singlefile <- mergeXDataReturn$singlefile
md5sumList <- mergeXDataReturn$md5sumList
sampleNamesList <- mergeXDataReturn$sampleNamesList
chromTIC <- mergeXDataReturn$chromTIC
chromBPI <- mergeXDataReturn$chromBPI
chromTIC_adjusted <- mergeXDataReturn$chromTIC_adjusted
chromBPI_adjusted <- mergeXDataReturn$chromBPI_adjusted

cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tDRAW GRAPHICS\n")

register(SerialParam())
if (!exists("chromTIC") || is.null(chromTIC)) {
  cat("\t\t\tCompute TIC\n"); chromTIC <- chromatogram(xdata, aggregationFun = "sum")
}
if (!exists("chromBPI") || is.null(chromBPI)) {
  cat("\t\t\tCompute BPI\n"); chromBPI <- chromatogram(xdata, aggregationFun = "max")
}

if (!is.null(chromTIC_adjusted)) chromTIC <- chromTIC_adjusted
if (!is.null(chromBPI_adjusted)) chromBPI <- chromBPI_adjusted

getPlotChromatogram(chromTIC, xdata, pdfname = "TICs.pdf", aggregationFun = "sum")
getPlotChromatogram(chromBPI, xdata, pdfname = "BPIs.pdf", aggregationFun = "max")

cat("\n\n")

# ----- EXPORT -----

cat("\tXCMSnExp OBJECT INFO\n")
print(xdata)
cat("\n\n")


cat("\tDONE\n")
