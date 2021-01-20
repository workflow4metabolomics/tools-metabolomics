#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 201919016 2DNmrAnnotation_1.0.0.R
## Marie Tremblay-Franco
## MetaboHUB: The French Infrastructure for Metabolomics and Fluxomics
## www.metabohub.fr/en
## marie.tremblay-franco@toulouse.inra.fr

runExampleL <- FALSE

if(runExampleL) {
##------------------------------
## Example of arguments
##------------------------------
}


##------------------------------
## Options
##------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

##------------------------------
## Constants
##------------------------------
topEnvC <- environment()
flagC <- "\n"


##-------------------------
## Input parameters reading
##-------------------------

##------------------------------
## R libraries laoding
##------------------------------
library(batch)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(tidyr)
library(curl)
library(jsonlite)
## library(tidyverse)

if(!runExampleL)
    argLs <- parseCommandArgs(evaluate=FALSE)
logFile <- argLs[["logOut"]]
sink(logFile)

cat("\tPACKAGE INFO\n")
sessionInfo()

##------------------------------
## Functions
##------------------------------
source_local <- function(fname)
{
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	source(paste(base_dir, fname, sep="/"))
}
#Import the different functions
source_local("annotationRmn2D.R")
source_local("annotationRmn2DGlobale.R")
source_local("viridis.R")

## Input parameter values
fileToAnnotate <- argLs[[1]]
  # Constraints values
ph <- argLs$pH
field <- argLs$magneticField

  # Chosen sequence(s)
cosy <- 0
hmbc <- 0
hsqc <- 0
jres <- 0
tocsy <- 0

if (argLs$cosy_2dsequences=='yes')
{
  cosy <- 1
  peakforestSpectra <- readLines(curl("https://metabohub.peakforest.org/rest/spectra/nmr2d/search?query=cosy&token=9131jq9l8gsjn1j14t351h716u&max=500"))
  peakforestSpectra <- fromJSON(peakforestSpectra, simplifyDataFrame = TRUE)
  if (ph != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$sampleNMRTubeConditionsMetadata$potentiaHydrogenii==ph, ]
  if (field != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$analyzerNMRSpectrometerDevice$magneticFieldStrenght==field, ]
  
  if (nrow(peakforestSpectra) != 0)
  {
    BdDReference_COSY <- peakforestSpectra$peaks
    names(BdDReference_COSY) <- str_split(peakforestSpectra[, 1], simplify=TRUE, pattern=";")[, 1]
    for (k in 1:length(BdDReference_COSY))
    {
      peakforestSpectra_df <- data.frame(ppm.dim1=BdDReference_COSY[[k]][, 2], ppm.dim2=BdDReference_COSY[[k]][, 1], 
                                         BdDReference_COSY[[k]][, 3:ncol(BdDReference_COSY[[k]])])
      BdDReference_COSY[[k]] <- peakforestSpectra_df
    }
  }
  else{
    stop("No COSY spectra correspond to requested pH and/or magnetic field", call.=FALSE)
  }
  rm(peakforestSpectra)
  rm(peakforestSpectra_df)
}

if (argLs$hmbc_2dsequences=='yes')
{
  hmbc <- 1
  peakforestSpectra <- readLines(curl("https://metabohub.peakforest.org/rest/spectra/nmr2d/search?query=hmbc&token=9131jq9l8gsjn1j14t351h716u&max=500"))
  peakforestSpectra <- fromJSON(peakforestSpectra, simplifyDataFrame = TRUE)
  if (ph != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$sampleNMRTubeConditionsMetadata$potentiaHydrogenii==ph, ]
  if (field != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$analyzerNMRSpectrometerDevice$magneticFieldStrenght==field, ]
  
  if (nrow(peakforestSpectra) != 0)
  {
    BdDReference_HMBC <- peakforestSpectra$peaks
    names(BdDReference_HMBC) <- str_split(peakforestSpectra[, 1], simplify=TRUE, pattern=";")[, 1]
    peakforestSpectra_df <- data.frame()
    for (k in 1:length(BdDReference_HMBC))
    {  
      peakforestSpectra_df <- data.frame(ppm.dim1=BdDReference_HMBC[[k]][, 2], ppm.dim2=BdDReference_HMBC[[k]][, 1], 
                                         BdDReference_HMBC[[k]][, 3:ncol(BdDReference_HMBC[[k]])])
      BdDReference_HMBC[[k]] <- peakforestSpectra_df
    }
  }
  else{
    stop("No HMBC spectra correspond to requested pH and/or magnetic field", call.=FALSE)
  }
  rm(peakforestSpectra)
  rm(peakforestSpectra_df)
}

if (argLs$hsqc_2dsequences=='yes')
{
  hsqc <- 1
  peakforestSpectra <- readLines(curl("https://metabohub.peakforest.org/rest/spectra/nmr2d/search?query=hsqc&token=9131jq9l8gsjn1j14t351h716u&max=500"))
  peakforestSpectra <- fromJSON(peakforestSpectra, simplifyDataFrame = TRUE)
  if (ph != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$sampleNMRTubeConditionsMetadata$potentiaHydrogenii==ph, ]
  if (field != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$analyzerNMRSpectrometerDevice$magneticFieldStrenght==field, ]
  
  if (nrow(peakforestSpectra) != 0)
  {
    BdDReference_HSQC <- peakforestSpectra$peaks
    names(BdDReference_HSQC) <- str_split(peakforestSpectra[, 1], simplify=TRUE, pattern=";")[, 1]
    for (k in 1:length(BdDReference_HSQC))
    {
      peakforestSpectra_df <- data.frame(ppm.dim1=BdDReference_HSQC[[k]][, 2], ppm.dim2=BdDReference_HSQC[[k]][, 1], 
                                         BdDReference_HSQC[[k]][, 3:ncol(BdDReference_HSQC[[k]])])
      BdDReference_HSQC[[k]] <- peakforestSpectra_df
    }
  }
  else{
    stop("No HSQC spectra correspond to requested pH and/or magnetic field", call.=FALSE)
  }  
  rm(peakforestSpectra)
  rm(peakforestSpectra_df)
}

if (argLs$jres_2dsequences=='yes')
{
  jres <- 1
  peakforestSpectra <- readLines(curl("https://metabohub.peakforest.org/rest/spectra/nmr2d/search?query=jres&token=9131jq9l8gsjn1j14t351h716u&max=500"))
  peakforestSpectra <- fromJSON(peakforestSpectra, simplifyDataFrame = TRUE)
  if (ph != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$sampleNMRTubeConditionsMetadata$potentiaHydrogenii==ph, ]
  if (field != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$analyzerNMRSpectrometerDevice$magneticFieldStrenght==field, ]
  
  if (nrow(peakforestSpectra) != 0)
  {
    BdDReference_JRES <- peakforestSpectra$peaks
    names(BdDReference_JRES) <- str_split(peakforestSpectra[, 1], simplify=TRUE, pattern=";")[, 1]
    for (k in 1:length(BdDReference_JRES))
    {
      peakforestSpectra_df <- data.frame(ppm.dim1=BdDReference_JRES[[k]][, 2], ppm.dim2=BdDReference_JRES[[k]][, 1], 
                                         BdDReference_JRES[[k]][, 3:ncol(BdDReference_JRES[[k]])])
      BdDReference_JRES[[k]] <- peakforestSpectra_df
    }
  }
  else{
    stop("No JRES spectra correspond to requested pH and/or magnetic field", call.=FALSE)
  }  
  rm(peakforestSpectra)
  rm(peakforestSpectra_df)
}

if (argLs$tocsy_2dsequences=='yes')
{
  tocsy <- 1
  peakforestSpectra <- readLines(curl("https://metabohub.peakforest.org/rest/spectra/nmr2d/search?query=tocsy&token=9131jq9l8gsjn1j14t351h716u&max=500"))
  peakforestSpectra <- fromJSON(peakforestSpectra, simplifyDataFrame = TRUE)
  if (ph != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$sampleNMRTubeConditionsMetadata$potentiaHydrogenii==ph, ]
  if (field != 0)
    peakforestSpectra <- peakforestSpectra[peakforestSpectra$analyzerNMRSpectrometerDevice$magneticFieldStrenght==field, ]
  
  if (nrow(peakforestSpectra) != 0)
  {
    BdDReference_TOCSY <- peakforestSpectra$peaks
    names(BdDReference_TOCSY) <- str_split(peakforestSpectra[, 1], simplify=TRUE, pattern=";")[, 1]
    for (k in 1:length(BdDReference_TOCSY))
    {
      peakforestSpectra_df <- data.frame(ppm.dim1=BdDReference_TOCSY[[k]][, 2], ppm.dim2=BdDReference_TOCSY[[k]][, 1], 
                                         BdDReference_TOCSY[[k]][, 3:ncol(BdDReference_TOCSY[[k]])])
      BdDReference_TOCSY[[k]] <- peakforestSpectra_df
    }
  }
  else{
    stop("No TOCSY spectra correspond to requested pH and/or magnetic field", call.=FALSE)
  }
  rm(peakforestSpectra)
  rm(peakforestSpectra_df)
}

if (argLs$cosy_2dsequences=='no' & argLs$hmbc_2dsequences=='no' & argLs$hsqc_2dsequences=='no' & argLs$jres_2dsequences=='no' & 
    argLs$tocsy_2dsequences=='no')
  stop("No chosen sequence. You have to choose at least 1 sequence", call.=FALSE)


  # User database


  # Allowed chemical shifts
tolPpm1 <- argLs$tolppm1
tolPpm2HJRes <- argLs$tolppmJRES
tolPpm2C <- argLs$tolppm2
  # Threshold to remove metabolites (probability score < threshold)
seuil <- argLs$threshold
# Remove metabolites when multiple assignations?
unicite <- str_to_upper(argLs$unicity)

## Output paramater values
AnnotationGraph <- argLs[["AnnotationGraph"]]

print(argLs)

## ANNOTATION
st0=Sys.time()
pdf(AnnotationGraph,onefile=TRUE)
annotationMelange <- annotationRmn2DGlobale(fileToAnnotate, tolPpm1=tolPpm1, tolPpm2HJRes=tolPpm2HJRes, 
                                             tolPpm2C=tolPpm2C, cosy=cosy, hmbc=hmbc, hsqc=hsqc, 
                                             jres=jres, tocsy=tocsy, seuil=seuil, unicite=unicite)
dev.off()

if (cosy==1)
{
  write.table(annotationMelange$COSY$liste_resultat, file=argLs[["annotationCOSY"]], quote=FALSE, 
              row.names=FALSE,sep="\t")
  if (nrow(annotationMelange$COSY$listing_ppm_commun) != 0)
      write.table(annotationMelange$COSY$listing_ppm_commun, file=argLs[["ppmCommunCOSY"]], quote=FALSE, 
                  row.names=FALSE,sep="\t")
}

if (hmbc==1)
{
  write.table(annotationMelange$HMBC$liste_resultat, file=argLs[["annotationHMBC"]], quote=FALSE, 
              row.names=FALSE,sep="\t")
  if (nrow(annotationMelange$HMBC$listing_ppm_commun) != 0)
    write.table(annotationMelange$HMBC$listing_ppm_commun, file=argLs[["ppmCommunHMBC"]], quote=FALSE, 
                row.names=FALSE,sep="\t")
}

if (hsqc==1)
{
  write.table(annotationMelange$HSQC$liste_resultat, file=argLs[["annotationHSQC"]], quote=FALSE, 
              row.names=FALSE,sep="\t")
  if (nrow(annotationMelange$HSQC$listing_ppm_commun) != 0)
    write.table(annotationMelange$HSQC$listing_ppm_commun, file=argLs[["ppmCommunHSQC"]], quote=FALSE, 
                row.names=FALSE,sep="\t")
}

if (jres==1)
{
  write.table(annotationMelange$JRES$liste_resultat, file=argLs[["annotationJRES"]], quote=FALSE, 
              row.names=FALSE,sep="\t")
  if (nrow(annotationMelange$JRES$listing_ppm_commun) != 0)
    write.table(annotationMelange$JRES$listing_ppm_commun, file=argLs[["ppmCommunJRES"]], quote=FALSE, 
                row.names=FALSE,sep="\t")
}

if (tocsy==1)
{
  write.table(annotationMelange$TOCSY$liste_resultat, file=argLs[["annotationTOCSY"]], quote=FALSE, 
              row.names=FALSE,sep="\t")
  if (nrow(annotationMelange$TOCSY$listing_ppm_commun) != 0)
    write.table(annotationMelange$TOCSY$listing_ppm_commun, file=argLs[["ppmCommunTOCSY"]], quote=FALSE, 
                row.names=FALSE,sep="\t")
}

## Combinaison de sequences
if (cosy + jres + hmbc + hsqc + tocsy > 1)
{
  write.table(annotationMelange$combination, file=argLs[["annotationCombination"]], quote=FALSE, 
              row.names=FALSE,sep="\t")
}
st1=Sys.time()
print(st1-st0)

## Ending
##--------
cat("\nEnd of '2D NMR annotation' Galaxy module call: ", as.character(Sys.time()), sep = "")
sink()
options(stringsAsFactors = strAsFacL)
rm(list = ls())
