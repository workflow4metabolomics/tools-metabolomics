#!/usr/bin/env Rscript

# R version 4.3.3
# Conversion tool mzXMLShaper.R (former cdf2mzml)
# Author: Quentin RUIN
# Creation: 25/04/2023
# Last updated: 09/09/2024

cat("\nJob starting time:\n", format(Sys.time(), "%a %d %b %Y %X"), "\n\n")

library(mzR)
library(msdata)
library("W4MRUtils")
library("tools")

args <- W4MRUtils::parse_args(args = commandArgs())

cat("\n\n--------------------------------------------------------------------",
    "\nParameters used by the 'mz(X)MLShaper' tool:\n\n")
cat("--------------------------------------------------------------------\n\n")
print(args)

inputfilename <- args[[1]]
outputfileformat <- args[[2]]
outputfilename <- args[[3]]
spectrum <- mzR::openMSfile(inputfilename)

## Get the spectra
pks <- mzR::spectra(spectrum)

## Get the header
hdr <- mzR::header(spectrum)

  if (outputfileformat == 'mzml') 
  {  
    writeMSData(pks, file = outputfilename, outformat = 'mzml', header = hdr)
  }

  if (outputfileformat == 'mzXml') 
  { 
    writeMSData(pks, file = outputfilename, outformat = 'mzxml', header = hdr)
  }


cat("\n--------------------------------------------------------------------",
    "\nInformation about R (version, Operating System, attached or loaded packages):\n\n")
sessionInfo()
cat("--------------------------------------------------------------------\n",
    "\nJob ending time:\n", format(Sys.time(), "%a %d %b %Y %X"))
