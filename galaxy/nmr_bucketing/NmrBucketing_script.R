################################################################################################
# SPECTRA BUCKETING AND INTEGRATION FROM RAW BRUKER FILES                                      #
# User : Galaxy                                                                                #
# Original data : --                                                                           #
# Starting date : 20-10-2014                                                                   #
# Version 1 : 18-12-2014                                                                       #
# Version 2 : 07-01-2015                                                                       #
#                                                                                              #
# Input files : files in included in user-defined directory                                    #
################################################################################################
NmrBucketing <- function(directory,leftBorder = 10.0,rightBorder = 0.5,bucketSize = 0.04,exclusionZones,exclusionZonesBorders=NULL,
                         graph=c("None","Overlay","One_per_individual"),nomFichier,savLog.txtC = NULL) 
{
  ## Option
  ##---------------
  strAsFacL <- options()$stringsAsFactors
  options(stingsAsFactors = FALSE)
  options(warn = -1)
  
  
  ## Constants
  ##---------------
  topEnvC <- environment()
  flgC <- "\n"
  
  ## Log file (in case of integration into Galaxy)
  ##----------------------------------------------
  if(!is.null(savLog.txtC))
    sink(savLog.txtC, append = TRUE)
  
  ## Functions definition
  ##---------------------  
    ## RAW BRUKER FILE READING FUNCTION
  NmRBrucker_read <- function(DataDir,SampleSpectrum)
  {
    
    bruker.get_param <- function (ACQ,paramStr)
    {
      regexpStr <- paste("^...",paramStr,"=",sep="")
      as.numeric(gsub("^[^=]+= ","" ,ACQ[which(simplify2array(regexpr(regexpStr,ACQ))>0)]))
    }
    
    ACQFILE <- "acqus"
    SPECFILE <- paste(DataDir,"/1r",sep="")
    PROCFILE <- paste(DataDir,"/procs",sep="")
    
    ACQ <- readLines(ACQFILE)
    TD      <- bruker.get_param(ACQ,"TD")
    SW      <- bruker.get_param(ACQ,"SW")
    SWH     <- bruker.get_param(ACQ,"SW_h")
    DTYPA   <- bruker.get_param(ACQ,"DTYPA")
    BYTORDA <- bruker.get_param(ACQ,"BYTORDA")
    #ENDIAN = ifelse( BYTORDA==0, "little", "big")
    ENDIAN <- "little"
    SIZE = ifelse( DTYPA==0, 4, 8)
    
    PROC <- readLines(PROCFILE)
    OFFSET <- bruker.get_param(PROC,"OFFSET")
    SI <- bruker.get_param(PROC,"SI")
    
    to.read = file(SPECFILE,"rb")
    maxTDSI = max(TD,SI)
    #  signal<-rev(readBin(to.read, what="int",size=SIZE, n=TD, signed = TRUE, endian = ENDIAN))
    signal<-rev(readBin(to.read, what="int",size=SIZE, n=maxTDSI, signed = TRUE, endian = ENDIAN))
    close(to.read)
    
    td <- length(signal)
    
    #  dppm <- SW/(TD-1)
    dppm <- SW/(td-1)
    pmax <- OFFSET
    pmin <- OFFSET - SW
    ppmseq <- seq(from=pmin, to=pmax, by=dppm)
    signal <- 100*signal/max(signal)
    
    SampleSpectrum <- cbind(ppmseq,signal)
    return(SampleSpectrum)
  }
  
    ## SPECTRUM BUCKETING
  NmrBrucker_bucket <- function(spectrum)
  {
    # Initialisations
    b <- 1
    j <- 1
    # Variable number
    J <- round((spectrum[1,1]-spectrum[dim(spectrum)[1],1])/bucketSize)
    f.bucket <- matrix(rep(0,J*2),ncol=2)
    colnames(f.bucket) <- c("Bucket",FileNames[i])
    
    
    # Data bucketing
    while (j < dim(spectrum)[1])
    {
      # chemical shift
      BUB <- spectrum[j,1]
      
      # In zone exclusion?
      exclusion.in <- FALSE
      if (!is.null(exclusionZonesBorders))
      {
        for (k in 1:nrow(exclusion.zone.m))
             if (BUB <= exclusion.zone.m[k,1] && exclusion.zone.m[k,2] < BUB)
                 exclusion.in <- TRUE
      }

      if (exclusion.in)
        j <- j + 1
      
      if (!exclusion.in)
        # Bucketing
      {
        BLB <- BUB - bucketSize
        bucket <- spectrum[j,]
        while (j < dim(spectrum)[1] && spectrum[j,1] > BLB)
        {
          j <- j + 1
          if (spectrum[j,1] > BLB)
            bucket <- rbind(bucket,spectrum[j,])
        }
        
        # Integration (trapezoid method)
        s <- cumtrapz(bucket[,1],bucket[,2])
        f.bucket[b,] <- c(round(mean(bucket[,1]),3),abs(s[length(s)][[1]]))
        
        # Next bucket boundary
        BUB <- spectrum[j,1]
        b <- b + 1
      }
    }
    return(f.bucket)
  }
  
    
  # File names
  FileNames <- list.files(directory)
  n <- length(FileNames)
  
  # Exclusion zones
##  if (exclusionZones == "yes")
  if (!is.null(exclusionZonesBorders))
  {
    exclusion.zone.m <- matrix(exclusionZonesBorders[[1]],nrow=1)
    if (length(exclusionZonesBorders) > 1)
      for (k in 2:length(exclusionZonesBorders))
        exclusion.zone.m <- rbind(exclusion.zone.m,exclusionZonesBorders[[k]])
  }

  # Reading and Bucketing
  directory <- paste(directory,"/",sep="")

  i <- 1
  while (i <= n)
  {
    # File reading
    SampleDir <- paste(directory,FileNames[i],"/1/",sep="")
    setwd(SampleDir)
    DataDir <- "pdata/1"

    rawSpectrum <- NmRBrucker_read(DataDir,rawSpectrum)

    orderedSpectrum <- rawSpectrum[order(rawSpectrum[,1],decreasing=T), ]
    
    # Removal of chemical shifts > leftBorder or < rightBorder boundaries
    truncatedSpectrum <- orderedSpectrum[orderedSpectrum[,1] < leftBorder & orderedSpectrum[,1] > rightBorder, ]
    truncatedSpectrum[,1] <- round(truncatedSpectrum[,1],3)
    
    # Bucketing
    spectrum.bucket <- NmrBrucker_bucket(truncatedSpectrum)
    
    # spectrum Concatenation
    if (i == 1)
      bucketedSpectra <- spectrum.bucket
    if (i > 1)
      bucketedSpectra <- cbind(bucketedSpectra,spectrum.bucket[,2])
    colnames(bucketedSpectra)[i+1] <- FileNames[i]
    
    # Next sample
    rm(spectrum.bucket)
    i <- i +1
  }
  identifiants <- gsub("([- , * { } | \\[ ])","_",colnames(bucketedSpectra)[-1])
  colnames(bucketedSpectra) <- c(colnames(bucketedSpectra)[1],identifiants)

  bucketedSpectra <- bucketedSpectra[bucketedSpectra[,1]!=0,]
  rownames(bucketedSpectra) <- paste("B",bucketedSpectra[,1],sep="")
  bucketedSpectra <- bucketedSpectra[,-1]
  
  # Metadata matrice outputs
  sampleMetadata <- data.frame(1:n)
  rownames(sampleMetadata) <- colnames(bucketedSpectra)
  colnames(sampleMetadata) <- "SampleOrder"
  
  variableMetadata <- data.frame(1:nrow(bucketedSpectra))
  rownames(variableMetadata) <- rownames(bucketedSpectra)
  colnames(variableMetadata) <- "VariableOrder"

  # Directory
  cd(directory)  
  
  # Bucketed spectra graph
  if (graph != "None")
  {
    # Graphic Device opening
    pdf(nomFichier,onefile=TRUE)

    if (graph == "Overlay")
    {
      x <- 1:length(BucketedData[,1])
      ymax <- max(bucketedSpectra)
      plot(x,BucketedData[,1],ylim=c(0,ymax),type='l',col=1,xlab="",xaxt="n",ylab="Intensity")
      # x-axis labels
      axis(1, at=seq(1,length(x),by=50),labels=gsub("B","",rownames(BucketedData)[seq(1,length(x),by=50)]), las=2)
      for (i in 2:ncol(bucketedSpectra))
      {
        spectre <- bucketedSpectra[,i]
        lines(spectre,col=i)
      }
#      legend(0,ymax,lty=c(1,1),legend=colnames(bucketedSpectra),col=1:ncol(bucketedSpectra))
    }
    else
    {
      for (i in 1:ncol(bucketedSpectra))
      {
        x <- 1:length(BucketedData[,1])
        plot(x,bucketedSpectra[,i],type='l',col=1,xlab="",xaxt="n",ylab="Intensity")
        axis(1, at=seq(1,length(x),by=50),labels=gsub("B","",rownames(BucketedData)[seq(1,length(x),by=50)]), las=2)
      }
    }
    dev.off()
  }
  return(list(bucketedSpectra,sampleMetadata,variableMetadata)) # ,truncatedSpectrum_matrice
}


#################################################################################################################
## Typical function call
#################################################################################################################
## StudyDir <- "K:/PROJETS/Metabohub/Bruker/Tlse_BPASourisCerveau/"
## upper <- 9.5
## lower <- 0.8
## bucket.width <- 0.01
## exclusion <- TRUE
## exclusion.zone <- list(c(5.1,4.5))
## graphique <- "Overlay"
## nomFichier <- "Tlse_BPASourisCerveau_NmrBucketing_graph.pdf"
## tlse_cerveaupnd21.bucket <- NmrBucketing(StudyDir,upper,lower,bucket.width,exclusion,exclusion.zone,graphique,nomFichier)
## write.table(tlse_cerveaupnd21.bucket,file=paste(StudyDir,"Tlse_BPASourisCerveau_NmrBucketing_dataMatrix.tsv",sep=""),
##             quote=FALSE,row.nmaes=FALSE,sep="\t")
#################################################################################################################
