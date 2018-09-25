#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 070115_NmrBucketing2galaxy_v1.R
## Marie Tremblay-Franco
## MetaboHUB: The French Infrastructure for Metabolomics and Fluxomics
## www.metabohub.fr/en
## marie.tremblay-franco@toulouse.inra.fr

runExampleL <- FALSE

if(runExampleL) {
##------------------------------
## Example of arguments
##------------------------------
argLs <- list(StudyDir = "Tlse_BPASourisCerveau",
              upper = "10.0",
              lower = "0.50",
              bucket.width = "0.01",
              exclusion = "TRUE",
              exclusion.zone = list(c(6.5,4.5)),
              graph="Overlay")

argLs <- c(argLs,
           list(dataMatrixOut = paste(directory,"_NmrBucketing_dataMatrix.tsv",sep=""),
                sampleMetadataOut = paste(directory,"_NmrBucketing_sampleMetadata.tsv",sep=""),
                variableMetadataOut = paste(directory,"_NmrBucketing_variableMetadata.tsv",sep=""),
                graphOut = paste(directory,"_NmrBucketing_graph.pdf",sep=""),
                logOut = paste(directory,"_NmrBucketing_log.txt",sep="")))
}

##------------------------------
## Options
##------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)


##------------------------------
## Libraries laoding
##------------------------------
# For parseCommandArgs function
library(batch)
# For cumtrapz function
library(pracma)

# R script call
source_local <- function(fname)
{
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	source(paste(base_dir, fname, sep="/"))
}
#Import the different functions
source_local("NmrBucketing_script.R")
source_local("DrawSpec.R")

##------------------------------
## Errors ?????????????????????
##------------------------------


##------------------------------
## Constants
##------------------------------
topEnvC <- environment()
flagC <- "\n"


##------------------------------
## Script
##------------------------------
if(!runExampleL)
    argLs <- parseCommandArgs(evaluate=FALSE)


## Parameters Loading
##-------------------
  # Inputs
if (!is.null(argLs[["zipfile"]])){
	fileType="zip"
	zipfile= argLs[["zipfile"]]
	directory=unzip(zipfile, list=F)
	directory=paste(getwd(),strsplit(directory[1],"/")[[1]][2],sep="/")
} else if (!is.null(argLs[["tsvfile"]])){
	fileType="tsv"
	directory <- read.table(argLs[["tsvfile"]],check.names=FALSE,header=TRUE,sep="\t")
}

leftBorder <- argLs[["left_border"]]
rightBorder <- argLs[["right_border"]]
bucketSize <- argLs[["bucket_width"]]
exclusionZones <- argLs[["zone_exclusion_choices.choice"]]

exclusionZonesBorders <- NULL
if (!is.null(argLs$zone_exclusion_left))
{
   for(i in which(names(argLs)=="zone_exclusion_left"))
   {
     exclusionZonesBorders <- c(exclusionZonesBorders,list(c(argLs[[i]],argLs[[i+1]])))
   }
}

graphique <- argLs[["graphType"]]

  # Outputs
nomGraphe <- argLs[["graphOut"]]
dataMatrixOut <- argLs[["dataMatrixOut"]]
logFile <- argLs[["logOut"]]
if (fileType=="zip")
{
  sampleMetadataOut <- argLs[["sampleOut"]]
  variableMetadataOut <- argLs[["variableOut"]]
}

## Checking arguments
##-------------------
error.stock <- "\n"

if(length(error.stock) > 1)
  stop(error.stock)


## Computation
##------------
outputs <- NmrBucketing(fileType=fileType, fileName=directory, leftBorder=leftBorder, rightBorder=rightBorder, bucketSize=bucketSize,
						exclusionZones=exclusionZones, exclusionZonesBorders=exclusionZonesBorders, graph=graphique, nomFichier=nomGraphe,
						savLog.txtC=logFile)
data_bucket1 <- outputs[[1]]
somme <- apply(data_bucket1, 1, sum)
data_bucket <- data_bucket1[somme !=0, ]

data_sample <- outputs[[2]]
data_variable <- outputs[[3]]
ppm <- outputs[[4]]
ppm <- round(ppm,2)

## Graphical outputs
##------------------
if (graphique != "None")
{
  excludedZone <- NULL
  for (c in 1:length(exclusionZonesBorders))
  {
    excludedZone <- c(excludedZone,exclusionZonesBorders[[c]])
    excludedZone <- sort(excludedZone)
  }
  nbZones <- length(excludedZone)/2
  n <- length(excludedZone)
  
  # Graphic Device opening
  pdf(nomGraphe,onefile=TRUE)
  
  if (graphique == "Overlay")
  {
    # Global spectral window
    spectra <- data.frame(t(data_bucket))
    drawSpec(spectra,xlab="", ylab="Intensity", main="")
    
    ## Zoomed spectral window depending on exclusion zone(s)
    if (nbZones != 0)
    {
      BInf <- excludedZone[n]
      if (round(BInf,1) == BInf)
      {
        BInf <- BInf+0.01
      }
      spectra <- data.frame(t(data_bucket[1:(which(ppm == BInf)[[1]]),]))
      drawSpec(spectra,xlab="", ylab="Intensity", main="")			
      n <- n - 1
      
      while (n >= nbZones & nbZones > 1)
      {
        BInf <- excludedZone[n-1]
        if (round(BInf,1) > BInf)
        {
          BInf <- BInf+0.01
        }
        spectra <- data.frame(t(data_bucket[(which(ppm == excludedZone[n])[[1]]):(which(ppm == BInf)[[1]]),]))
        drawSpec(spectra,xlab="", ylab="Intensity", main="")
        n <- n - 2
      }
      
      BInf <- excludedZone[1]
      if (round(BInf,1) <= BInf)
      {
        BInf <- BInf+0.01
      }
      spectra <- data.frame(t(data_bucket[(which(ppm == BInf)[[1]]):nrow(data_bucket),]))
      drawSpec(spectra,xlab="", ylab="Intensity", main="")
    }
  }
  else
  {
    for (i in 1:ncol(data_bucket))
    {
      par(mfrow=c((nbZones+2),1))
      n <- length(excludedZone)
      spectra <- t(data_bucket[,i])
	  names(spectra) <- rownames(data_bucket)
      plot(1:length(spectra), spectra, type='l', xlab="", ylab="Intensity", main=colnames(data_bucket)[i], xaxt = "n")
	  xPos <- 1
	  nAxisPos <- 4
	  startP <- length(nAxisPos) 
	  endP <- nrow(data_bucket)
	  GraphRange <- c(startP:endP)
	  tempVal = trunc(length(GraphRange)/nAxisPos)
	  xPos = c(0:nAxisPos) * tempVal
	  axis(1, at = xPos, labels = rownames(data_bucket)[xPos + startP])
     
      ## Zoomed spectral window depending on exclusion zone(s)
      if (nbZones != 0)
      {
        BInf <- excludedZone[n]
        if (round(BInf,1) == BInf)
        {
          BInf <- BInf+0.01
        }
        spectra <- t(data_bucket[1:(which(ppm == BInf)[[1]]),i])
		names(spectra) <- rownames(data_bucket)[1:(which(ppm == BInf)[[1]])]
		plot(1:length(spectra), spectra, type='l',xlab="", ylab="Intensity", main="", xaxt = "n")			
		xPos <- 1
		nAxisPos <- 4
		startP <- length(nAxisPos) 
		endP <- length(spectra)
		GraphRange <- c(startP:endP)
		tempVal = trunc(length(GraphRange)/nAxisPos)
		xPos = c(0:nAxisPos) * tempVal
		axis(1, at = xPos, labels = rownames(data_bucket)[xPos + startP])
        n <- n - 1
        
        while (n >= nbZones & nbZones > 1)
        {
          BInf <- excludedZone[n-1]
          if (round(BInf,1) > BInf)
          {
            BInf <- BInf+0.01
          }
          spectra <- t(data_bucket[(which(ppm == excludedZone[n])[[1]]):(which(ppm == BInf)[[1]]),i])
		  names(spectra) <- rownames(data_bucket)[(which(ppm == excludedZone[n])[[1]]):(which(ppm == BInf)[[1]])]
          plot(1:length(spectra), spectra, type='l',xlab="", ylab="Intensity", main="", xaxt = "n")
		  xPos <- 1
		  nAxisPos <- 4
		  startP <- length(nAxisPos) 
		  endP <- length(spectra)
		  GraphRange <- c(startP:endP)
		  tempVal = trunc(length(GraphRange)/nAxisPos)
		  xPos = c(0:nAxisPos) * tempVal
		  axis(1, at = xPos, labels = rownames(data_bucket)[xPos + startP])
          n <- n - 2
        }
        
        BInf <- excludedZone[1]
        if (round(BInf,1) <= BInf)
        {
          BInf <- BInf+0.01
        }
        spectra <- t(data_bucket[(which(ppm == BInf)[[1]]):nrow(data_bucket),i])
		names(spectra) <- rownames(data_bucket)[(which(ppm == BInf)[[1]]):nrow(data_bucket)]
        plot(1:length(spectra), spectra, type='l',xlab="", ylab="Intensity", main="", xaxt = "n")
		xPos <- 1
		nAxisPos <- 4
		startP <- length(nAxisPos) 
		endP <- length(spectra)
		GraphRange <- c(startP:endP)
		tempVal = trunc(length(GraphRange)/nAxisPos)
		xPos = c(0:nAxisPos) * tempVal
		axis(1, at = xPos, labels = rownames(data_bucket)[xPos + startP])
      }
    }
  }
  dev.off()
}
## Saving
##-------
  # Data
data_bucket <- cbind(rownames(data_bucket),data_bucket)
colnames(data_bucket) <- c("Bucket",colnames(data_bucket)[-1])
write.table(data_bucket,file=argLs$dataMatrixOut,quote=FALSE,row.names=FALSE,sep="\t")
  # Sample
data_sample <- cbind(rownames(data_sample),data_sample)
colnames(data_sample) <- c("Sample",colnames(data_sample)[-1])
write.table(data_sample,file=argLs$sampleOut,quote=FALSE,row.names=FALSE,sep="\t")
  # Variable
data_variable <- cbind(rownames(data_variable),data_variable)
colnames(data_variable) <- c("Bucket",colnames(data_variable)[-1])
write.table(data_variable,file=argLs$variableOut,quote=FALSE,row.names=FALSE,sep="\t")


## Ending
##---------------------

cat("\nEnd of 'NMR bucketing' Galaxy module call: ", as.character(Sys.time()), sep = "")

## sink(NULL)

options(stringsAsFactors = strAsFacL)

rm(list = ls())
