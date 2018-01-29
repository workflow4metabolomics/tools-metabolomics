#!/usr/bin/env Rscript

library(batch) ## parseCommandArgs
library(lme4)     ## mixed model computing
library(Matrix)
library(lmerTest) ## computing pvalue and lsmeans from results of lme4 package 
library(multtest) ## multiple testing

setwd("//Tls-tox-nas/home$/jfmartin/Mes documents/PROJETS/PROG/mixed_model/galaxy/mixmodel/mixmodel4repeated_measures/")
source("mixmodel_script_trycatch.R")
source("diagmfl.R")

options(stringsAsFactors = FALSE)

## loading
##--------
setwd("//Tls-tox-nas/home$/jfmartin/Mes documents/PROJETS/PROG/mixed_model/dataMTF2")

datMN <- t(as.matrix(read.table("Choo_G1G3.Alignment.Bucketing.PQN_datamatrix_Bal.tsv",
                                check.names = FALSE,
                                header = TRUE,
                                row.names = 1,
                                sep = "\t")))

samDF <- read.table("Choo_G1G3.Alignment.Bucketing.PQN_samplemetadata_Bal.tsv",
                    check.names = FALSE,
                    header = TRUE,
                    sep = "\t")
row.names(samDF) <- samDF$Sample

varDF <- read.table("Choo_G1G3.Alignment.Bucketing.PQN_variablemetadata.tsv",
                    check.names = FALSE,
                    header = TRUE,
                    sep = "\t")
row.names(varDF) <- varDF$Bucket

## checking
##---------

##------------------------------
## Computation
##------------------------------


varDF <- lmixedm(datMN = datMN,
                     samDF = samDF,
                     varDF = varDF,
                     fixfact    = "Treatment",
                     time       = "Date",
                     subject    = "Individu",
                     dffOption  = "Kenward-Roger",
					           visu		= "yes", 
                     least.confounded = FALSE,
                     outlier.limit = 3,
                     logtr ="none",
					           pdfC     = "resuAllraw.pdf"
					      )


## affectation des variables pour test
samDF <- samDF 
varDF <- varDF 
fixfact    <- "Treatment" 
time       <- "Date" 
subject    <- "Individu" 
dffOption  <- "Kenward-Roger" 
visu		<- "yes"  
least.confounded <- FALSE 
outlier.limit <- 3 
logtr <-"none" 
pdfC     <- "resuraw.pdf"


##------------------------------
## Ending


##------------------------------


## saving
##--------

varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                          varDF)

write.table(varDF,
            file = "test_raw_variable_out.txt",
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

## closing
##--------

cat("\nEnd of '", modNamC, "' Galaxy module call: ",
    as.character(Sys.time()), "\n", sep = "")

sink()

options(stringsAsFactors = strAsFacL)

rm(list = ls())
