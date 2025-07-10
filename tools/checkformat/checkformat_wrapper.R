#!/usr/bin/env Rscript

library(batch) ## parseCommandArgs

source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep = "/"))
}

source_local("checkformat_script.R")

argVc <- unlist(parseCommandArgs(evaluate = FALSE))


## ------------------------------
## Initializing
## ------------------------------

## options
## --------

strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## constants
## ----------

modNamC <- "Check Format" ## module name

## log file
## ---------

sink(argVc["information"])

cat("\nStart of the '", modNamC, "' Galaxy module call: ",
  format(Sys.time(), "%a %d %b %Y %X"), "\n",
  sep = ""
)


## ------------------------------
## Computation
## ------------------------------


resLs <- readAndCheckF(
  argVc["dataMatrix_in"],
  argVc["sampleMetadata_in"],
  argVc["variableMetadata_in"],
  as.logical(argVc["makeNameL"])
)


## ------------------------------
## Ending
## ------------------------------


## dataMatrix

datMN <- resLs[["datMN"]]
datDF <- cbind.data.frame(
  dataMatrix = colnames(datMN),
  as.data.frame(t(datMN))
)
write.table(datDF,
  file = argVc[["dataMatrix_out"]],
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

## sampleMetadata

samDF <- resLs[["samDF"]]
samDF <- cbind.data.frame(
  sampleMetadata = rownames(samDF),
  samDF
)
write.table(samDF,
  file = argVc["sampleMetadata_out"],
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

## variableMetadata

varDF <- resLs[["varDF"]]
varDF <- cbind.data.frame(
  variableMetadata = rownames(varDF),
  varDF
)
write.table(varDF,
  file = argVc["variableMetadata_out"],
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

if (resLs[["chkL"]]) {
  if (resLs[["newL"]]) {
    cat("\nWarning: The sample and/or variable names or orders from the input tables have been modified\n(see the information file for details); please use the new output tables for your analyses.\n")
  } else {
    cat("\nThe input tables have a correct format and can be used for your analyses.\n")
  }
}

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

if (!resLs[["chkL"]]) {
  stop("Please check the generated 'information' file")
}


## closing
## --------

options(stringsAsFactors = strAsFacL)

rm(list = ls())
