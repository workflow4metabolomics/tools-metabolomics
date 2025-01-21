#!/usr/bin/env Rscript

## 070115_NmrBucketing2galaxy_v1.R
## Marie Tremblay-Franco
## MetaboHUB: The French Infrastructure for Metabolomics and Fluxomics
## www.metabohub.fr/en
## marie.tremblay-franco@toulouse.inra.fr

runExampleL <- FALSE


## ------------------------------
## Options
## ------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)
options(scipen = 999)


## ------------------------------
## Libraries laoding
## ------------------------------
# For parseCommandArgs function
library(batch)

# Constants
argv <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=", "", argv[grep("--file=", argv)])
prog.name <- basename(script.path)

# Print help
if (length(grep("-h", argv)) > 0) {
    cat(
        "Usage:", prog.name,
        "dataMatrix myDataMatrix.tsv",
        "scalingMethod PQN|QuantitativeVariable",
        "graphType None|Overlay|One_per_individual",
        "logOut myLog.txt",
        "dataMatrixOut myDataMatrixOutput.tsv",
        "graphOut myGraph.pdf",
        "\n"
    )
    quit(status = 0)
}

# R script call
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
# Import the different functions
source_local("NmrNormalization_script.R")
source_local("DrawSpec.R")


## ------------------------------
## Errors ?????????????????????
## ------------------------------


## ------------------------------
## Constants
## ------------------------------
topEnvC <- environment()
flagC <- "\n"


## ------------------------------
## Script
## ------------------------------
if (!runExampleL) {
    argLs <- parseCommandArgs(evaluate = FALSE)
}


## Parameters Loading
## -------------------
# Inputs
data <- read.table(argLs[["dataMatrix"]], check.names = FALSE, header = TRUE, sep = "\t", row.names = 1)
names <- rownames(data)
## Add a test to check if all values are numercical
if (!all(vapply(data, is.numeric, FUN.VALUE = FALSE))) {
    stop("Data are not numeric")
}
## Integer conversion to avoid stack overflow when computin the sum
data <- as.data.frame(lapply(data, as.numeric))
rownames(data) <- names

scaling <- argLs[["scalingMethod"]]
graphique <- argLs[["graphType"]]

if (scaling == "PQN") {
    metadataSample <- read.table(argLs[["sampleMetadata"]], check.names = FALSE, header = TRUE, sep = "\t")
    factor <- argLs[["factor"]]
    ControlGroup <- argLs[["controlGroup"]]
}
if (scaling == "QuantitativeVariable") {
    metadataSample <- read.table(argLs[["sampleMetadata"]], check.names = FALSE, header = TRUE, sep = "\t")
    factor <- argLs[["factor"]]
}

# Outputs
nomGraphe <- argLs[["graphOut"]]
dataMatrixOut <- argLs[["dataMatrixOut"]]
log <- argLs[["logOut"]]


## Checking R packages
## --------------------
sink(log)
cat("\tPACKAGE INFO\n")
pkgs <- c("batch")
for (pkg in pkgs) {
    suppressPackageStartupMessages(stopifnot(library(pkg, quietly = TRUE, logical.return = TRUE, character.only = TRUE)))
    cat(pkg, "\t", as.character(packageVersion(pkg)), "\n", sep = "")
}
cat("\n")


## Checking arguments
## -------------------
error.stock <- "\n"
if (length(error.stock) > 1) {
    stop(error.stock)
}


## Computation
## ------------
NormalizationResults <- NmrNormalization(
    dataMatrix = data, scalingMethod = scaling, sampleMetadata = metadataSample,
    bioFactor = factor, ControlGroup = ControlGroup,
    graph = graphique, nomFichier = nomGraphe, savLog.txtC = log
)

data_normalized <- NormalizationResults[[1]]


## Graphical outputs
## ------------------
if (graphique != "None") {
    # Graphic Device opening
    pdf(nomGraphe, onefile = TRUE)

    if (graphique == "Overlay") {
        # Global spectral window
        spectra <- data.frame(t(data_normalized))
        drawSpec(spectra, xlab = "", ylab = "Intensity", main = "")
    } else {
        for (i in 1:ncol(data_normalized))
        {
            spectra <- t(data_normalized[, i])
            drawSpec(spectra, xlab = "", ylab = "Intensity", main = colnames(data_normalized)[i])
        }
    }
    dev.off()
}


## Saving
## -------
# Data
data_normalized <- cbind(rownames(data_normalized), data_normalized)
colnames(data_normalized) <- c("Variable", colnames(data_normalized)[-1])
write.table(data_normalized, file = argLs$dataMatrixOut, quote = FALSE, row.names = FALSE, sep = "\t")


## Ending
## ---------------------
cat("\nEnd of 'Normalization' Galaxy module call: ", as.character(Sys.time()), sep = "")
sink()
options(stringsAsFactors = strAsFacL)
rm(list = ls())
