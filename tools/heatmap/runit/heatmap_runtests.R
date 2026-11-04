#!/usr/bin/env Rscript

## Package
## --------

library(RUnit)

## Constants
## ----------

testOutDirC <- "output"
argVc <- commandArgs(trailingOnly = FALSE)
scriptPathC <- sub("--file=", "", argVc[grep("--file=", argVc)])


## Functions
## -----------

## Reading tables (matrix or data frame)
readTableF <- function(fileC, typeC = c("matrix", "dataframe")[1]) {
    file.exists(fileC) || stop(paste0("No output file \"", fileC, "\"."))

    switch(typeC,
        matrix = return(t(as.matrix(read.table(
            file = fileC,
            header = TRUE,
            row.names = 1,
            sep = "\t",
            stringsAsFactors = FALSE
        )))),
        dataframe = return(read.table(
            file = fileC,
            header = TRUE,
            row.names = 1,
            sep = "\t",
            stringsAsFactors = FALSE
        ))
    )
}

## Call wrapper
wrapperCallF <- function(paramLs) {
    ## Set program path
    wrapperPathC <- file.path(dirname(scriptPathC), "..", "heatmap_wrapper.R")

    ## Set arguments
    argLs <- NULL
    for (parC in names(paramLs)) {
        argLs <- c(argLs, parC, paramLs[[parC]])
    }

    ## Call
    wrapperCallC <- paste(c(wrapperPathC, argLs), collapse = " ")

    if (.Platform$OS.type == "windows") {
        wrapperCallC <- paste("Rscript", wrapperCallC)
    }

    wrapperCodeN <- system(wrapperCallC)

    if (wrapperCodeN != 0) {
        stop("Error when running heatmap_wrapper.R.")
    }

    ## Get output
    outLs <- list()
    if ("data_matrix_out" %in% names(paramLs)) {
        outLs[["datMN"]] <- readTableF(paramLs[["data_matrix_out"]], "matrix")
    }
    if ("sample_metadata_out" %in% names(paramLs)) {
        outLs[["samDF"]] <- readTableF(paramLs[["sample_metadata_out"]], "dataframe")
    }
    if ("variable_metadata_out" %in% names(paramLs)) {
        outLs[["varDF"]] <- readTableF(paramLs[["variable_metadata_out"]], "dataframe")
    }
    if ("information" %in% names(paramLs)) {
        outLs[["infVc"]] <- readLines(paramLs[["information"]])
    }

    return(outLs)
}

## Setting default parameters
defaultArgF <- function(testInDirC) {
    defaultArgLs <- list()
    if (file.exists(file.path(dirname(scriptPathC), testInDirC, "dataMatrix.tsv"))) {
        defaultArgLs[["data_matrix_in"]] <- file.path(dirname(scriptPathC), testInDirC, "dataMatrix.tsv")
    }
    if (file.exists(file.path(dirname(scriptPathC), testInDirC, "sampleMetadata.tsv"))) {
        defaultArgLs[["sample_metadata_in"]] <- file.path(dirname(scriptPathC), testInDirC, "sampleMetadata.tsv")
    }
    if (file.exists(file.path(dirname(scriptPathC), testInDirC, "variableMetadata.tsv"))) {
        defaultArgLs[["variable_metadata_in"]] <- file.path(dirname(scriptPathC), testInDirC, "variableMetadata.tsv")
    }

    defaultArgLs[["data_matrix_out"]] <- file.path(dirname(scriptPathC), testOutDirC, "dataMatrix.tsv")
    defaultArgLs[["sample_metadata_out"]] <- file.path(dirname(scriptPathC), testOutDirC, "sampleMetadata.tsv")
    defaultArgLs[["variable_metadata_out"]] <- file.path(dirname(scriptPathC), testOutDirC, "variableMetadata.tsv")
    defaultArgLs[["figure"]] <- file.path(dirname(scriptPathC), testOutDirC, "figure.pdf")
    defaultArgLs[["information"]] <- file.path(dirname(scriptPathC), testOutDirC, "information.txt")

    defaultArgLs
}

## Main
## -----

## Create output folder
file.exists(testOutDirC) || dir.create(testOutDirC)

## Run tests
test.suite <- defineTestSuite("tests", dirname(scriptPathC), testFileRegexp = paste0("^.*_tests\\.R$"), testFuncRegexp = "^.*$")
isValidTestSuite(test.suite)
test.results <- runTestSuite(test.suite)
print(test.results)
