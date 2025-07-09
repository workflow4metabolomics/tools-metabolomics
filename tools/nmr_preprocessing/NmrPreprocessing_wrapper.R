#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 170116_NmrPreprocessing.R
## Manon Martin and Marie Tremblay-Franco

## ======================================================
## ======================================================
# Preamble
## ======================================================
## ======================================================

runExampleL <- FALSE


## ------------------------------
## Options
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## Libraries laoding
## ------------------------------
# library(batch)
suppressPackageStartupMessages(library(ptw))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(reshape2))

# In-house function for argument parsing
parse_args <- function() {
  args <- commandArgs()
  start <- which(args == "--args")[1] + 1
  if (is.na(start)) {
    return(list())
  }
  seq_by2 <- seq(start, length(args), by = 2)
  result <- as.list(args[seq_by2 + 1])
  names(result) <- args[seq_by2]
  return(result)
}

# R script call
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep = "/"))
}
# Import the different functions
source_local("NmrPreprocessing_script.R")
source_local("DrawFunctions.R")

## ------------------------------
## Script
## ------------------------------
runExampleL <- FALSE

if (!runExampleL) {
  #  argLs <- parseCommandArgs(evaluate=FALSE)
  argLs <- unlist(parse_args())
}

print(argLs[["logOut"]])

# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")
print(argLs)

## ------------------------------
## Errors ?????????????????????
## ------------------------------


## ------------------------------
## Constants
## ------------------------------
topEnvC <- environment()
flagC <- "\n"

## Starting
cat("\nStart of 'Preprocessing' Galaxy module call: ", as.character(Sys.time()), "\n", sep = "")


## ======================================================
## ======================================================
## Parameters Loading
## ======================================================
## ======================================================

# graphical inputs
FirstOPCGraph <- argLs[["FirstOPCGraph"]]
SSGraph <- argLs[["SSGraph"]]
ApodGraph <- argLs[["ApodGraph"]]
FTGraph <- argLs[["FTGraph"]]
SRGraph <- argLs[["SRGraph"]]
ZeroOPCGraph <- argLs[["ZeroOPCGraph"]]
BCGraph <- argLs[["BCGraph"]]
FinalGraph <- argLs[["FinalGraph"]]

# 1rst order phase correction ------------------------
# Inputs
## Data matrix
Fid_data0 <- read.table(argLs[["dataMatrixFid"]], header = TRUE, check.names = FALSE, sep = "\t")
# Fid_data0 <- Fid_data0[,-1]
Fid_data0 <- as.matrix(Fid_data0)

## Samplemetadata
samplemetadataFid <- read.table(argLs[["sampleMetadataFid"]], check.names = FALSE, header = TRUE, sep = "\t")
samplemetadataFid <- as.matrix(samplemetadataFid)

# water and solvent(s) correction ------------------------
# Inputs
lambda <- argLs[["lambda"]]

# apodization -----------------------------------------
# Inputs
phase <- 0
rectRatio <- 1 / 2
gaussLB <- 1
expLB <- 1
apodization <- argLs[["apodizationMethod"]]

if (apodization == "exp") {
  expLB <- argLs[["expLB"]]
} else if (apodization == "cos2") {
  phase <- argLs[["phase"]]
} else if (apodization == "hanning") {
  phase <- argLs[["phase"]]
} else if (apodization == "hamming") {
  phase <- argLs[["phase"]]
} else if (apodization == "blockexp") {
  rectRatio <- argLs[["rectRatio"]]
  expLB <- argLs[["expLB"]]
} else if (apodization == "blockcos2") {
  rectRatio <- argLs[["rectRatio"]]
} else if (apodization == "gauss") {
  rectRatio <- argLs[["rectRatio"]]
  gaussLB <- argLs[["gaussLB"]]
}

# Fourier transform ----------------------------------
# Inputs

# Zero Order Phase Correction -------------------------------
# Inputs

angle <- NULL
excludeZOPC <- NULL


zeroOrderPhaseMethod <- argLs[["zeroOrderPhaseMethod"]]

if (zeroOrderPhaseMethod == "manual") {
  angle <- argLs[["angle"]]
}

excludeZoneZeroPhase <- argLs[["excludeZoneZeroPhase.choice"]]
if (excludeZoneZeroPhase == "YES") {
  excludeZoneZeroPhaseList <- list()
  for (i in which(names(argLs) == "excludeZoneZeroPhase_left")) {
    excludeZoneZeroPhaseLeft <- argLs[[i]]
    excludeZoneZeroPhaseRight <- argLs[[i + 1]]
    excludeZoneZeroPhaseList <- c(excludeZoneZeroPhaseList, list(c(excludeZoneZeroPhaseLeft, excludeZoneZeroPhaseRight)))
  }
  excludeZOPC <- excludeZoneZeroPhaseList
}


# Internal referencering ----------------------------------
# Inputs
shiftTreshold <- 2
ppm <- TRUE
shiftReferencingRangeList <- NULL # fromto.RC
pctNearValue <- 0.02 # pc
rowindex_graph <- NULL
ppm_ref <- 0 # ppm.ref


# shiftReferencing <- argLs[["shiftReferencing"]]
# print(shiftReferencing)
#
# if (shiftReferencing=="YES") {
#
# shiftReferencingMethod <- argLs[["shiftReferencingMethod"]]
#
# if (shiftReferencingMethod == "thres")	{
# 	shiftTreshold <- argLs[["shiftTreshold"]]
# }

shiftReferencingRange <- argLs[["shiftReferencingRange"]]

if (shiftReferencingRange == "near0") {
  pctNearValue <- argLs[["pctNearValue"]]
}

if (shiftReferencingRange == "window") {
  shiftReferencingRangeList <- list()
  for (i in which(names(argLs) == "shiftReferencingRangeLeft"))
  {
    shiftReferencingRangeLeft <- argLs[[i]]
    shiftReferencingRangeRight <- argLs[[i + 1]]
    shiftReferencingRangeList <- c(shiftReferencingRangeList, list(c(shiftReferencingRangeLeft, shiftReferencingRangeRight)))
  }
}

shiftHandling <- argLs[["shiftHandling"]]

ppmvalue <- argLs[["ppmvalue"]]
# }

# Baseline Correction -------------------------------
# Inputs
lambdaBc <- argLs[["lambdaBc"]]
pBc <- argLs[["pBc"]]
epsilon <- argLs[["epsilon"]]

excludeBC <- NULL

excludeZoneBC <- argLs[["excludeZoneBC.choice"]]
if (excludeZoneBC == "YES") {
  excludeZoneBCList <- list()
  for (i in which(names(argLs) == "excludeZoneBC_left")) {
    excludeZoneBCLeft <- argLs[[i]]
    excludeZoneBCRight <- argLs[[i + 1]]
    excludeZoneBCList <- c(excludeZoneBCList, list(c(excludeZoneBCLeft, excludeZoneBCRight)))
  }
  excludeBC <- excludeZoneBCList
}

# transformation of negative values -------------------------------
# Inputs
NegativetoZero <- argLs[["NegativetoZero"]]

# Outputs
nomGraphe <- argLs[["graphOut"]]
# dataMatrixOut <- argLs[["dataMatrixOut"]]
log <- argLs[["logOut"]]

## Checking arguments
## -------------------
error.stock <- "\n"

if (length(error.stock) > 1) {
  stop(error.stock)
}

## ======================================================
## Computation
## ======================================================
pdf(nomGraphe, onefile = TRUE, width = 13, height = 13)

# FirstOrderPhaseCorrection ---------------------------------
# Fid_data <- GroupDelayCorrection(Fid_data0, Fid_info = samplemetadataFid, group_delay = NULL)

if (FirstOPCGraph == "YES") {
  title <- "FIDs after Group Delay Correction"
  DrawSignal(Fid_data,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}

# SolventSuppression ---------------------------------
Fid_data <- SolventSuppression(Fid_data, lambda.ss = lambda, ptw.ss = TRUE, plotSolvent = F, returnSolvent = F)

if (SSGraph == "YES") {
  title <- "FIDs after Solvent Suppression "
  DrawSignal(Fid_data,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}


# Apodization ---------------------------------
Fid_data <- Apodization(Fid_data,
  Fid_info = samplemetadataFid, DT = NULL,
  type.apod = apodization, phase = phase, rectRatio = rectRatio, gaussLB = gaussLB, expLB = expLB, plotWindow = F, returnFactor = F
)

if (ApodGraph == "YES") {
  title <- "FIDs after Apodization"
  DrawSignal(Fid_data,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}


# FourierTransform ---------------------------------
Spectrum_data <- FourierTransform(Fid_data, Fid_info = samplemetadataFid, reverse.axis = TRUE)


if (FTGraph == "YES") {
  title <- "Fourier transformed spectra"
  DrawSignal(Spectrum_data,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}


# if (FTGraph == "YES") {
#   title = "Fourier transformed spectra"
# DrawSignal(Spectrum_data, subtype = "stacked",
#            ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
#            xlab = "Frequency", num.stacked = 4,
#            main = title, createWindow=FALSE)
# }

# ZeroOrderPhaseCorrection ---------------------------------
Spectrum_data <- ZeroOrderPhaseCorrection(Spectrum_data,
  type.zopc = zeroOrderPhaseMethod,
  plot_rms = NULL, returnAngle = FALSE,
  createWindow = TRUE, angle = angle,
  plot_spectra = FALSE,
  ppm.zopc = TRUE, exclude.zopc = excludeZOPC
)


# InternalReferencing ---------------------------------
# if (shiftReferencing=="YES") {
Spectrum_data <- InternalReferencing(Spectrum_data, samplemetadataFid,
  method = "max", range = shiftReferencingRange,
  ppm.value = ppmvalue, shiftHandling = shiftHandling, ppm.ir = TRUE,
  fromto.RC = shiftReferencingRangeList, pc = pctNearValue
)

if (SRGraph == "YES") {
  title <- "Spectra after Shift Referencing"
  DrawSignal(Spectrum_data,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}

# }

if (ZeroOPCGraph == "YES") {
  title <- "Spectra after Zero Order Phase Correction"
  DrawSignal(Spectrum_data,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}


# BaselineCorrection ---------------------------------
Spectrum_data <- BaselineCorrection(Spectrum_data,
  ptw.bc = TRUE, lambda.bc = lambdaBc,
  p.bc = pBc, eps = epsilon, ppm.bc = TRUE,
  exclude.bc = excludeBC,
  returnBaseline = F
)



if (BCGraph == "YES") {
  title <- "Spectra after Baseline Correction"
  DrawSignal(Spectrum_data,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}

# if (BCGraph == "YES") {
# title = "Spectra after Baseline Correction"
# DrawSignal(Spectrum_data, subtype = "stacked",
#          ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
#          xlab = "Frequency", num.stacked = 4,
#          main = title, createWindow=FALSE)
# }

# NegativeValuesZeroing ---------------------------------
if (NegativetoZero == "YES") {
  Spectrum_data <- NegativeValuesZeroing(Spectrum_data)
}

if (FinalGraph == "YES") {
  title <- "Final preprocessed spectra"
  DrawSignal(Spectrum_data,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}

# invisible(dev.off())

# data_variable <- matrix(NA, nrow = 1, ncol = dim(Spectrum_data)[2], dimnames = list("ID", NULL))
# colnames(data_variable) <- colnames(Spectrum_data)
# data_variable[1,] <- colnames(data_variable)

data_variable <- matrix(NA, nrow = 1, ncol = dim(Spectrum_data)[2], dimnames = list("ID", NULL))
colnames(data_variable) <- colnames(Spectrum_data)
data_variable[1, ] <- colnames(data_variable)


## ======================================================
## ======================================================
## Saving
## ======================================================
## ======================================================

# Data Matrix
write.table(round(t(Re(Spectrum_data)), 6), file = argLs$dataMatrix, quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# Variable metadata
write.table(data_variable, file = argLs$variableMetadata, quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# log file
# write.table(t(data.frame(argLs)), file = argLs$logOut, col.names = FALSE, quote=FALSE)

# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")

argLs


## Ending
cat("\nVersion of R librairies")
sessionInfo()
cat("\nEnd of 'Preprocessing' Galaxy module call: ", as.character(Sys.time()), sep = "")

sink()

options(stringsAsFactors = strAsFacL)

rm(list = ls())
