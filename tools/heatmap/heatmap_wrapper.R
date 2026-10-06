#!/usr/bin/env Rscript

library(batch) ## parseCommandArgs

source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}

source_local("heatmap_script.R")

argVc <- unlist(parseCommandArgs(evaluate = FALSE))


## ------------------------------
## Initializing
## ------------------------------

## options
## --------

strAsFacL <- options()[["stringsAsFactors"]]
options(stringsAsFactors = FALSE)

## constants
## ----------

modNamC <- "Heatmap" ## module name

## log file
## ---------

sink(argVc["information"])

cat("\nStart of the '", modNamC, "' module: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n",
    sep = ""
)

## loading
## --------

proMN <- t(as.matrix(read.table(argVc["data_matrix_in"],
    check.names = FALSE,
    header = TRUE,
    row.names = 1,
    sep = "\t"
)))

obsDF <- read.table(argVc["sample_metadata_in"],
    check.names = FALSE,
    header = TRUE,
    row.names = 1,
    sep = "\t"
)

feaDF <- read.table(argVc["variable_metadata_in"],
    check.names = FALSE,
    header = TRUE,
    row.names = 1,
    sep = "\t"
)

## adding default parameter values
## --------------------------------


if (!("cor_met_c" %in% names(argVc))) {
    argVc["cor_met_c"] <- "pearson"
}
if (!("agg_met_c" %in% names(argVc))) {
    argVc["agg_met_c"] <- "ward"
}
if (!("col_c" %in% names(argVc))) {
    argVc["col_c"] <- "blueOrangeRed"
}
if (!("sca_l" %in% names(argVc))) {
    argVc["sca_l"] <- "TRUE"
}
if (!("cex_n" %in% names(argVc))) {
    argVc["cex_n"] <- "0.8"
}

## checking
## ---------

if (as.numeric(argVc["cut_sam_n"]) > nrow(proMN)) {
    stop("Number of sample clusters must be inferior to the number of samples")
}
if (as.numeric(argVc["cut_var_n"]) > ncol(proMN)) {
    stop("Number of variable clusters must be inferior to the number of variables")
}

## printing arguments
## -------------------

cat("\nArguments used:\n\n")
argMC <- as.matrix(argVc)
colnames(argMC) <- "value"
argDatVl <- grepl("\\.dat$", argVc) ## discarding dataset file names
if (sum(argDatVl)) {
    argMC <- argMC[!argDatVl, , drop = FALSE]
}
print(argMC)


## ------------------------------
## Computation
## ------------------------------


heaLs <- heatmapF(
    proMN = proMN,
    obsDF = obsDF,
    feaDF = feaDF,
    dis_c = argVc["dis_c"],
    cut_sam_n = as.numeric(argVc["cut_sam_n"]),
    cut_var_n = as.numeric(argVc["cut_var_n"]),
    fig.pdfC = argVc["figure"],
    cor_met_c = argVc["cor_met_c"],
    agg_met_c = argVc["agg_met_c"],
    col_c = argVc["col_c"],
    sca_l = as.logical(argVc["sca_l"]),
    cex_n = as.numeric(argVc["cex_n"])
)


## ------------------------------
## Ending
## ------------------------------


## saving
## -------

proDF <- cbind.data.frame(
    dataMatrix = colnames(heaLs[["proMN"]]),
    as.data.frame(t(heaLs[["proMN"]]))
)
write.table(proDF,
    file = argVc["data_matrix_out"],
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
)

obsDF <- cbind.data.frame(
    sampleMetadata = rownames(heaLs[["obsDF"]]),
    heaLs[["obsDF"]]
)
write.table(obsDF,
    file = argVc["sample_metadata_out"],
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
)

feaDF <- cbind.data.frame(
    variableMetadata = rownames(heaLs[["feaDF"]]),
    heaLs[["feaDF"]]
)
write.table(feaDF,
    file = argVc["variable_metadata_out"],
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
)

## Ending
## -------

cat("\nEnd of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n",
    sep = ""
)

cat("\n\n\n============================================================================")
cat("\nAdditional information about the call:\n")
cat("\n1) Parameters:\n")
print(cbind(value = argVc))

cat("\n2) Session Info:\n")
sessioninfo <- sessionInfo()
cat(sessioninfo$R.version$version.string, "\n")
cat("Main packages:\n")
for (pkg in names(sessioninfo$otherPkgs)) {
    cat(paste(pkg, packageVersion(pkg)), "\t")
}
cat("\n")
cat("Other loaded packages:\n")
for (pkg in names(sessioninfo$loadedOnly)) {
    cat(paste(pkg, packageVersion(pkg)), "\t")
}
cat("\n")

cat("============================================================================\n")

sink()

options(stringsAsFactors = strAsFacL)

rm(list = ls())
