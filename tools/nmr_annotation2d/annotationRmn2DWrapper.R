#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 201919016 2DNmrAnnotation_1.0.0.R
## Marie Tremblay-Franco
## MetaboHUB: The French Infrastructure for Metabolomics and Fluxomics
## www.metabohub.fr/en
## marie.tremblay-franco@toulouse.inra.fr

runExampleL <- FALSE

if (runExampleL) {
    ## ------------------------------
    ## Example of arguments
    ## ------------------------------
}


## ------------------------------
## Options
## ------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## ------------------------------
## Constants
## ------------------------------
topEnvC <- environment()
flagC <- "\n"


## -------------------------
## Input parameters reading
## -------------------------

## ------------------------------
## R libraries laoding
## ------------------------------
library(batch)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(tidyr)

if (!runExampleL) {
    argLs <- parseCommandArgs(evaluate = FALSE)
}
logFile <- argLs[["logOut"]]
sink(logFile)

cat("\tPACKAGE INFO\n")
sessionInfo()

## ------------------------------
## Functions
## ------------------------------
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
# Import the different functions
source_local("annotationRmn2D.R")
source_local("annotationRmn2DGlobale.R")
source_local("viridis.R")

## Input parameter values
fileToAnnotate <- argLs[[1]]
# Chosen sequence(s)
cosy <- 0
hmbc <- 0
hsqc <- 0
jres <- 0
tocsy <- 0
## sequences <- str_split(argLs[[2]], ",")[[1]]
## for (s in 1:length(sequences))
## {
##   argv <- commandArgs(trailingOnly = FALSE)
##   currentDir <- dirname(substring(argv[grep("--file=", argv)], 8))
##   if (sequences[s]=="cosy"){
##         cosy <- 1
##         load(paste(currentDir, "BdDReference_COSY.RData", sep="/"))
##   }else if(sequences[s]=="hmbc"){
##         hmbc <- 1
##         load(paste(currentDir, "BdDReference_HMBC.RData", sep="/"))
##   }else if(sequences[s]=="hsqc"){
##         hsqc <- 1
##         load(paste(currentDir, "BdDReference_HSQC.RData", sep="/"))
##   }else if(sequences[s]=="jres"){
##         jres <- 1
##         load(paste(currentDir, "BdDReference_JRES.RData", sep="/"))
##   }else if(sequences[s]=="tocsy"){
##         tocsy <- 1
##         load(paste(currentDir, "BdDReference_TOCSY.RData", sep="/"))
##   }else
##     stop("No chosen sequence", call.=FALSE)
## }

if (argLs[[2]] == "yes") {
    argv <- commandArgs(trailingOnly = FALSE)
    currentDir <- dirname(substring(argv[grep("--file=", argv)], 8))
    cosy <- 1
    load(paste(currentDir, "BdDReference_COSY.RData", sep = "/"))
}

if (argLs[[3]] == "yes") {
    argv <- commandArgs(trailingOnly = FALSE)
    currentDir <- dirname(substring(argv[grep("--file=", argv)], 8))
    jres <- 1
    load(paste(currentDir, "BdDReference_JRES.RData", sep = "/"))
}

if (argLs[[4]] == "yes") {
    argv <- commandArgs(trailingOnly = FALSE)
    currentDir <- dirname(substring(argv[grep("--file=", argv)], 8))
    hmbc <- 1
    load(paste(currentDir, "BdDReference_HMBC.RData", sep = "/"))
}

if (argLs[[5]] == "yes") {
    argv <- commandArgs(trailingOnly = FALSE)
    currentDir <- dirname(substring(argv[grep("--file=", argv)], 8))
    hsqc <- 1
    load(paste(currentDir, "BdDReference_HSQC.RData", sep = "/"))
}

if (argLs[[6]] == "yes") {
    argv <- commandArgs(trailingOnly = FALSE)
    currentDir <- dirname(substring(argv[grep("--file=", argv)], 8))
    tocsy <- 1
    load(paste(currentDir, "BdDReference_TOCSY.RData", sep = "/"))
}

if (argLs[[2]] == "no" & argLs[[3]] == "no" & argLs[[4]] == "no" & argLs[[5]] == "no" & argLs[[6]] == "no") {
    stop("No chosen sequence", call. = FALSE)
}


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
st0 <- Sys.time()
pdf(AnnotationGraph, onefile = TRUE)
annotationMelange <- annotationRmn2DGlobale(fileToAnnotate,
    tolPpm1 = tolPpm1, tolPpm2HJRes = tolPpm2HJRes,
    tolPpm2C = tolPpm2C, cosy = cosy, hmbc = hmbc, hsqc = hsqc,
    jres = jres, tocsy = tocsy, seuil = seuil, unicite = unicite
)
dev.off()

if (cosy == 1) {
    write.table(annotationMelange$COSY$liste_resultat,
        file = argLs[["annotationCOSY"]], quote = FALSE,
        row.names = FALSE, sep = "\t"
    )
    if (nrow(annotationMelange$COSY$listing_ppm_commun) != 0) {
        write.table(annotationMelange$COSY$listing_ppm_commun,
            file = argLs[["ppmCommunCOSY"]], quote = FALSE,
            row.names = FALSE, sep = "\t"
        )
    }
}

if (hmbc == 1) {
    write.table(annotationMelange$HMBC$liste_resultat,
        file = argLs[["annotationHMBC"]], quote = FALSE,
        row.names = FALSE, sep = "\t"
    )
    if (nrow(annotationMelange$HMBC$listing_ppm_commun) != 0) {
        write.table(annotationMelange$HMBC$listing_ppm_commun,
            file = argLs[["ppmCommunHMBC"]], quote = FALSE,
            row.names = FALSE, sep = "\t"
        )
    }
}

if (hsqc == 1) {
    write.table(annotationMelange$HSQC$liste_resultat,
        file = argLs[["annotationHSQC"]], quote = FALSE,
        row.names = FALSE, sep = "\t"
    )
    if (nrow(annotationMelange$HSQC$listing_ppm_commun) != 0) {
        write.table(annotationMelange$HSQC$listing_ppm_commun,
            file = argLs[["ppmCommunHSQC"]], quote = FALSE,
            row.names = FALSE, sep = "\t"
        )
    }
}

if (jres == 1) {
    write.table(annotationMelange$JRES$liste_resultat,
        file = argLs[["annotationJRES"]], quote = FALSE,
        row.names = FALSE, sep = "\t"
    )
    if (nrow(annotationMelange$JRES$listing_ppm_commun) != 0) {
        write.table(annotationMelange$JRES$listing_ppm_commun,
            file = argLs[["ppmCommunJRES"]], quote = FALSE,
            row.names = FALSE, sep = "\t"
        )
    }
}

if (tocsy == 1) {
    write.table(annotationMelange$TOCSY$liste_resultat,
        file = argLs[["annotationTOCSY"]], quote = FALSE,
        row.names = FALSE, sep = "\t"
    )
    if (nrow(annotationMelange$TOCSY$listing_ppm_commun) != 0) {
        write.table(annotationMelange$TOCSY$listing_ppm_commun,
            file = argLs[["ppmCommunTOCSY"]], quote = FALSE,
            row.names = FALSE, sep = "\t"
        )
    }
}

## Combinaison de sequences
if (cosy + jres + hmbc + hsqc + tocsy > 1) {
    write.table(annotationMelange$combination,
        file = argLs[["annotationCombination"]], quote = FALSE,
        row.names = FALSE, sep = "\t"
    )
}
st1 <- Sys.time()
print(st1 - st0)

## Ending
## --------
cat("\nEnd of '2D NMR annotation' Galaxy module call: ", as.character(Sys.time()), sep = "")
sink()
options(stringsAsFactors = strAsFacL)
rm(list = ls())
