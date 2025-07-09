#####################################################
# Mini tools for Galaxy scripting
# Coded by: M.Petera, 
# - -
# R functions to use in R scripts and wrappers
# to make things easier (lightening code, reducing verbose...)
# - -
# V0: script structure + first functions
# V1: addition of functions to handle special characters in identifiers
#####################################################


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Function to call packages without printing all the verbose
# (only getting the essentials, like warning messages for example)

shyLib <- function(...){
	for(i in 1:length(list(...))){
		suppressPackageStartupMessages(library(list(...)[[i]],character.only=TRUE))
	}
}

#example: shyLib("xcms","pcaMethods")



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fonction pour sourcer les scripts R requis
# /!\ ATTENTION : actuellement la fonction n'est pas chargee au lancement du script,
# il faut donc la copier-coller dans le wrapper R pour pouvoir l'utiliser. 

if(FALSE){
source_local <- function(...){
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	for(i in 1:length(list(...))){
		source(paste(base_dir, list(...)[[i]], sep="/"))
	}
}
}

#example: source_local("filter_script.R","RcheckLibrary.R")



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Functions to stock identifiers before applying make.names() and
# to reinject it into final matrices
# Note: it reproduces the original order of datasets' identifiers
# - - -
# stockID: stocks original identifiers and original order
# -> needs checked data regarding table match
# reproduceID: reinjects original identifiers and original order into final tables
# -> function to be used at the very end, when exporting tables

stockID <- function(dataMatrix, Metadata, Mtype){
  # dataMatrix = data.frame containing dataMatrix
  # Metadata = data.frame containing sampleMetadata or variableMetadata
  # Mtype = "sample" or "variable" depending on Metadata content
  cname <- colnames(dataMatrix)[1]
      # dataMatrix temporary-stock + transfo - - - -
  if(Mtype=="sample"){
    id.ori <- colnames(dataMatrix)[-1] 
    colnames(dataMatrix) <- make.names(colnames(dataMatrix)) 
  }
  if(Mtype=="variable"){
    id.ori <- dataMatrix[,1] 
    dataMatrix[,1] <- make.names(dataMatrix[,1]) 
  }
      # global stock - - - - - - - - - - - - - - - -
  id.new <- data.frame(order.ori=c(1:length(Metadata[,1])),Metadata[,1],
					   id.new=make.names(Metadata[,1]),id.ori,
					   id.new.DM=make.names(id.ori),stringsAsFactors=FALSE)
  colnames(id.new)[c(2,4)] <- c(colnames(Metadata)[1],cname)
      # Metadata transfo + returning data - - - - - 
  Metadata[,1] <- make.names(Metadata[,1]) 
  return(list(id.match=id.new, dataMatrix=dataMatrix, Metadata=Metadata))
}
#example: A<-stockID(myDM,mysM,"sample") ; myDM<-A$dataMatrix ; mysM<-A$Metadata ; A<-A$id.match

reproduceID <- function(dataMatrix, Metadata, Mtype, id.match){
  # dataMatrix = data.frame containing dataMatrix
  # Metadata = data.frame containing sampleMetadata or variableMetadata
  # Mtype = "sample" or "variable" depending on Metadata content
  # id.match = 'id.match' element produced by stockID
      #Metadada - - - - - - - - - - - - - - 
  temp.table <- id.match[,c(1,2,3)]
  ## Removing deleted rows
  for(i in 1:(dim(id.match)[1])){
	if(!(temp.table[i,3]%in%Metadata[,1])){temp.table[i,1] <- 0}
  }
  if(length(which(temp.table[,1]==0))!=0){
	temp.table <- temp.table[-c(which(temp.table[,1]==0)),]
  }
  ## Restoring original identifiers and order
  temp.table <- merge(x=temp.table,y=Metadata,by.x=3,by.y=1)
  temp.table <- temp.table[order(temp.table$order.ori),]
  Metadata <- temp.table[,-c(1,2)]
  rownames(Metadata) <- NULL
      #dataMatrix - - - - - - - - - - - - - 
  rownames(dataMatrix)<-dataMatrix[,1]
  if(Mtype=="sample"){
    dataMatrix <- t(dataMatrix[,-1])
  }
  temp.table <- id.match[,c(1,4,5)]
  ## Removing deleted rows
  for(i in 1:(dim(id.match)[1])){
	if(!(temp.table[i,3]%in%rownames(dataMatrix))){temp.table[i,1] <- 0}
  }
  if(length(which(temp.table[,1]==0))!=0){
	temp.table <- temp.table[-c(which(temp.table[,1]==0)),]
  }
  ## Restoring original identifiers and order
  temp.table <- merge(x=temp.table,y=dataMatrix,by.x=3,by.y=0)
  temp.table <- temp.table[order(temp.table$order.ori),]
  if(Mtype=="variable"){
	dataMatrix <- temp.table[,-c(1,2,4)]
	colnames(dataMatrix)[1] <- colnames(id.match)[4]
  } else {
	rownames(temp.table) <- temp.table[,3]
	temp.table <- t(temp.table[,-c(1,2,3)])
	dataMatrix <- data.frame(rownames(temp.table),temp.table)
	colnames(dataMatrix)[1] <- colnames(id.match)[4]
  }
  rownames(dataMatrix) <- NULL
      # return datasets - - - - - - - - - - - 
  return(list(dataMatrix=dataMatrix, Metadata=Metadata))
}
#example: B<-reproduceID(myDM,mysM,"sample",A) ; myDM<-B$dataMatrix ; mysM<-B$Metadata



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
