#!/usr/bin/env Rscript

library(batch) ## necessary for parseCommandArgs function

## ------------------------------
## test help option
## ------------------------------

# Prog. constants
argv.help <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=", "", argv.help[grep("--file=", argv.help)])
prog.name <- basename(script.path)

# Test Help
if (length(grep("-h", argv.help)) > 0) {
    cat(
        "Usage: Rscript ",
        prog.name,
        "{args} \n",
        "parameters: \n",
        "\tdataMatrix {file}: set the input data matrix file (mandatory) \n",
        "\tsampleMetadata {file}: set the input sample metadata file (mandatory) \n",
        "\tvariableMetadata {file}: set the input variable metadata file (mandatory) \n",
        "\tmethod {opt}: set the method; can set to \"all_loess_pool\" or \"all_loess_sample\" (mandatory) \n",
        "\tspan {condition}: set the span condition; (mandatory) \n",
        "\tdataMatrix_out {file}: set the output data matrix file (mandatory) \n",
        "\tvariableMetadata_out {file}: set the output variable metadata file (mandatory) \n",
        "\tgraph_output {file}: set the output graph file (mandatory) \n",
        "\trdata_output {file}: set the output Rdata file (mandatory) \n",
        "\tbatch_col_name {val}: the column name for batch. Default value is \"batch\".\n",
        "\tinjection_order_col_name {val}: the column name for the injection order. Default value is \"injectionOrder\".\n",
        "\tsample_type_col_name {val}: the column name for the sample types. Default value is \"sampleType\".\n",
        "\tsample_type_tags {val}: the tags used inside the sample type column, defined as key/value pairs separated by commas (example: blank=blank,pool=pool,sample=sample).\n",
        "\n"
    )
    quit(status = 0)
}

## ------------------------------
## init. params
## ------------------------------

args <- parseCommandArgs(evaluate = FALSE) # interpretation of arguments given in command line as an R list of objects

# Set default col names
if (!"batch_col_name" %in% names(args)) {
    args[["batch_col_name"]] <- "batch"
}
if (!"injection_order_col_name" %in% names(args)) {
    args[["injection_order_col_name"]] <- "injectionOrder"
}
if (!"sample_type_col_name" %in% names(args)) {
    args[["sample_type_col_name"]] <- "sampleType"
}
if (!"sample_type_tags" %in% names(args)) {
    args[["sample_type_tags"]] <- "blank=blank,pool=pool,sample=sample"
}

# Parse sample type tags
sample.type.tags <- list()
for (kv in strsplit(strsplit(args$sample_type_tags, ",")[[1]], "=")) {
    sample.type.tags[[kv[[1]]]] <- kv[[2]]
}
if (!all(c("pool", "blank", "sample") %in% names(sample.type.tags))) {
    stop("All tags pool, blank and sample must be defined in option sampleTypeTags.")
}
args$sample_type_tags <- sample.type.tags

## ------------------------------
## init. functions
## ------------------------------

source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}

## Import the different functions
source_local("batch_correction_all_loess_script.R")

argVc <- unlist(args)

##  argVc["method"] is either 'all_loess_pool' or 'all_loess_sample'
##  alternative version developped by CEA
##  all variables are treated with loess
##  the reference observations for loess are either 'pool'
## ('all_loess_pool') or 'sample' ('all_loess_sample')


## ------------------------------
## Initializing
## ------------------------------

## options
## --------

strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## libraries
## ----------

suppressMessages(library(ropls))

if (packageVersion("ropls") < "1.4.0") {
    stop("Please use 'ropls' versions of 1.4.0 and above")
}

## constants
## ----------

modNamC <- "Batch correction" ## module name

## log file
## ---------

## sink(argVc["information"]) ## not implemented

cat("\nStart of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n",
    sep = ""
)

## loading
## --------

rawMN <- t(as.matrix(read.table(argVc["dataMatrix"],
    header = TRUE,
    comment.char = "",
    row.names = 1,
    sep = "\t"
)))

samDF <- read.table(argVc["sampleMetadata"],
    header = TRUE,
    comment.char = "",
    row.names = 1,
    sep = "\t"
)

varDF <- read.table(argVc["variableMetadata"],
    check.names = FALSE,
    header = TRUE,
    comment.char = "",
    row.names = 1,
    sep = "\t"
) ## not used; for compatibility only

refC <- tolower(gsub("all_loess_", "", argVc["method"]))

spnN <- as.numeric(argVc["span"])

## checking
## ---------

stopifnot(refC %in% c("pool", "sample"))
refC <- args$sample_type_tags[[refC]]

if (refC == args$sample_type_tags$pool &&
    !any(args$sample_type_tags$pool %in% samDF[, args$sample_type_col_name])) {
    stop("No 'pool' found in the 'sampleType' column; use the samples as normalization reference instead")
}

refMN <- rawMN[samDF[, args$sample_type_col_name] == refC, ]
refNasZerVl <- apply(
    refMN, 2,
    function(refVn) {
        all(sapply(
            refVn,
            function(refN) {
                is.na(refN) || refN == 0
            }
        ))
    }
)

if (sum(refNasZerVl)) {
    refNasZerVi <- which(refNasZerVl)
    cat("The following variables have 'NA' or 0 values in all reference samples; they will be removed from the data:\n", sep = "")
    rawMN <- rawMN[, !refNasZerVl, drop = FALSE]
    varDF <- varDF[!refNasZerVl, , drop = FALSE]
}

## ------------------------------
## Computation
## ------------------------------


## ordering (batch and injection order)
## -------------------------------------

samDF[, "ordIniVi"] <- 1:nrow(rawMN)
ordBatInjVi <- order(samDF[, args$batch_col_name], samDF[, args$injection_order_col_name])
rawMN <- rawMN[ordBatInjVi, ]
samDF <- samDF[ordBatInjVi, ]

## signal drift and batch-effect correction
## -----------------------------------------

nrmMN <- shiftBatchCorrectF(
    rawMN,
    samDF,
    refC,
    spnN
)

## figure
## -------

cat("\nPlotting\n")

pdf(argVc["graph_output"], onefile = TRUE, width = 11, height = 7)
plotBatchF(rawMN, samDF, spnN)
plotBatchF(nrmMN, samDF, spnN)
dev.off()

## returning to initial order
## ---------------------------

ordIniVi <- order(samDF[, "ordIniVi"])
nrmMN <- nrmMN[ordIniVi, ]
samDF <- samDF[ordIniVi, ]
samDF <- samDF[, colnames(samDF) != "ordIniVi", drop = FALSE]


## ------------------------------
## Ending
## ------------------------------


## saving
## -------

datMN <- nrmMN

datDF <- cbind.data.frame(
    dataMatrix = colnames(datMN),
    as.data.frame(t(datMN))
)
write.table(datDF,
    file = argVc["dataMatrix_out"],
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
)

varDF <- cbind.data.frame(
    variableMetadata = rownames(varDF),
    varDF
) ## not modified; for compatibility only
write.table(varDF,
    file = argVc["variableMetadata_out"],
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
)


res <- list(
    dataMatrix_raw = rawMN,
    dataMatrix_normalized = nrmMN,
    sampleMetadata = samDF
)
save(res,
    file = argVc["rdata_output"]
) ## for compatibility

## closing
## --------

cat("\nEnd of '", modNamC, "' Galaxy module call: ",
    as.character(Sys.time()), "\n",
    sep = ""
)

## sink()

options(stringsAsFactors = strAsFacL)


rm(argVc)
