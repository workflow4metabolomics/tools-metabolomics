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

## sink(argLs[["logOut"]])


## Parameters Loading
##-------------------
  # Inputs
if (!is.null(argLs[["zipfile"]])){
	fileType="zip"
	zipfile= argLs[["zipfile"]]
	directory=unzip(zipfile, list=F)
	directory=paste(getwd(),strsplit(directory[1],"/")[[1]][2],sep="/")
} else if (!is.null(argLs[["library"]])){
	fileType="zip"
	directory=argLs[["library"]]
    	if(!file.exists(directory)){
		error_message=paste("Cannot access the directory :",directory,".Please verify if the directory exists or not.")
		print(error_message)
		stop(error_message)
	}
}
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
sampleMetadataOut <- argLs[["sampleOut"]]
variableMetadataOut <- argLs[["variableOut"]]
logFile <- argLs[["logOut"]]

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
data_bucket <- outputs[[1]]
data_sample <- outputs[[2]]
data_variable <- outputs[[3]]

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
