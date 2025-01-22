#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 29122017_asics_wrapper.R
## Remi Servien, Patrick Tardivel, Marie Tremblay-Franco and Gaelle Lefort
## marie.tremblay-franco@inra.fr

runExampleL <- FALSE

## ------------------------------
## Options
## ------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)


## ------------------------------
## Libraries loading
## ------------------------------
# ParseCommandArgs function
library(batch)
library(ASICS)



# R script call
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
# Import the different functions
source_local("DrawSpec.R")


## ------------------------------
## Errors ?????????????????????
## ------------------------------


## ------------------------------
## Constants
## ------------------------------
topEnvC <- environment()
flagC <- "\n"


## ------------------------------
## Script
## ------------------------------
if (!runExampleL) {
    argLs <- parseCommandArgs(evaluate = FALSE)
}

# Standards loading
load(argLs[["standards"]])

## Parameters Loading
## -------------------
# Inputs
## Spectrum to annotate
zipfile <- argLs[["zipfile"]]
directory <- unzip(zipfile, list = F)
directory <- paste(getwd(), strsplit(directory[1], "/")[[1]][2], sep = "/")


## Exclusion zone(s)
exclusionZones <- argLs[["zone_exclusion_choices.choice"]]
exclusionZonesBorders <- NULL
if (!is.null(argLs$zone_exclusion_left)) {
    for (i in which(names(argLs) == "zone_exclusion_left"))
    {
        #     exclusionZonesBorders <- c(exclusionZonesBorders,list(c(argLs[[i]],argLs[[i+1]])))
        exclusionZonesBorders <- c(exclusionZonesBorders, argLs[[i]], argLs[[i + 1]])
    }
    exclusionZonesBorders <- matrix(exclusionZonesBorders, byrow = T, ncol = 2)
}
if (is.null(argLs$zone_exclusion_left)) {
    exclusionZonesBorders <- matrix(c(0, 0), byrow = T, ncol = 2)
}
## Maximal allowed shift
shift <- argLs[["shift"]]

## Graphical zone(s)
graphicalZones <- argLs[["zone_graphical_choices.choice"]]
graphicalZonesBorders <- NULL
if (!is.null(argLs$zone_exclusion_left)) {
    for (i in which(names(argLs) == "zone_graphical_left"))
    {
        graphicalZonesBorders <- c(graphicalZonesBorders, list(c(argLs[[i]], argLs[[i + 1]])))
    }
}

# Outputs
logOut <- argLs[["logOut"]]
proportionEstimation <- argLs[["proportionEstimation"]]
graphOut <- argLs[["graphOut"]]

sink(logOut)
cat("\tPACKAGE INFO\n")
# pkgs=c("batch", "ASICS")
pkgs <- c("batch", "ASICS")
for (pkg in pkgs) {
    suppressPackageStartupMessages(stopifnot(library(pkg, quietly = TRUE, logical.return = TRUE, character.only = TRUE)))
    cat(pkg, "\t", as.character(packageVersion(pkg)), "\n", sep = "")
}
cat("\n")


## Checking arguments
## -------------------
error.stock <- "\n"
if (length(error.stock) > 1) {
    stop(error.stock)
}


## Computation
## ------------
# annotation.Asics <- ASICS(directory, exclusion.areas=matrix(exclusionZonesBorders, byrow=T, ncol=2),
#                            max.shift=shift, which.spectra="last", library.metabolites=NULL,
#                            threshold.noise=0.02, seed=1234, nb.iter.signif=400)
annotation.Asics <- ASICS(directory,
    exclusion.areas = exclusionZonesBorders,
    max.shift = shift, which.spectra = "last", library.metabolites = NULL,
    threshold.noise = 0.02, seed = 1234, nb.iter.signif = 400
)


## Saving
## -------
# Identified metabolites
metabolites.estimation <- present_metabolites(annotation.Asics)
colnames(metabolites.estimation) <- c("Metabolite", colnames(metabolites.estimation)[-1])
write.table(metabolites.estimation, file = argLs$proportionEstimation, row.names = FALSE, quote = FALSE, sep = "\t")


## Graphical display
## ------------------
# Raw and annotated spectra comparison
pdf(graphOut, onefile = TRUE)

## Graphical output: overlay of raw and estimated spectra
ppm.metabolites.estimation <- data.frame(
    round(ppm_grid(annotation.Asics), 3),
    original_mixture(annotation.Asics)
)
colnames(ppm.metabolites.estimation) <- c("PPM", "EstimatedProportion")
ppm.metabolites.estimation <- ppm.metabolites.estimation[order(ppm.metabolites.estimation[, 1], decreasing = T), ]

mix <- data.frame(t(ppm.metabolites.estimation[, 2]))
colnames(mix) <- ppm.metabolites.estimation[, 1]
ppm <- ppm.metabolites.estimation[, 1]

estimatedMix <- data.frame(round(ppm_grid(annotation.Asics), 3), reconstituted_mixture(annotation.Asics))
colnames(estimatedMix) <- c("PPM", "EstimatedProportion")
estimatedMix <- estimatedMix[order(estimatedMix[, 1], decreasing = T), ]
estimatedMix <- estimatedMix[, 2]

## Whole spectra
GraphRange <- 1:ncol(mix)
tempVal <- trunc(length(GraphRange) / 10)
xPos <- c(10:0) * tempVal
plot(1:ncol(mix), mix, type = "l", xlab = "", main = "", xaxt = "n", ylab = "")
axis(1, at = xPos, labels = colnames(mix)[xPos + 1])
lines(estimatedMix, col = "red")
legend("topleft", legend = c("Real Mixture", "Estimated Composition"), lty = c(1, 1), col = c("black", "red"))

## Zoomed spectral window depending on user-selected zone(s)
graphical.zone.length <- length(graphicalZonesBorders)
if (graphical.zone.length != 0) {
    #   par(mfrow=c(2,1))
    for (g in 1:graphical.zone.length)
    {
        print(g)
        plot(1:length((which(round(as.numeric(colnames(mix)), 2) == graphicalZonesBorders[[g]][1])[1]):(which(round(as.numeric(colnames(mix)), 2) == max(graphicalZonesBorders[[g]][2], 0.5))[1])),
            mix[(which(round(as.numeric(colnames(mix)), 2) == graphicalZonesBorders[[g]][1])[1]):(which(round(as.numeric(colnames(mix)), 2) == max(graphicalZonesBorders[[g]][2], 0.5))[1])],
            type = "l", xlab = "", ylab = "Intensity", main = "", xaxt = "n"
        )
        lines(estimatedMix[(which(round(as.numeric(colnames(mix)), 2) == graphicalZonesBorders[[g]][1])[1]):(which(round(as.numeric(colnames(mix)), 2) == max(graphicalZonesBorders[[g]][2], 0.5))[1])], col = "red")

        xPos <- 1
        nAxisPos <- 4
        startP <- length(nAxisPos)
        endP <- length((which(round(as.numeric(colnames(mix)), 2) == graphicalZonesBorders[[g]][1])[1]):(which(round(as.numeric(colnames(mix)), 2) == max(graphicalZonesBorders[[g]][2], 0.5))[1]))
        GraphRange <- c(startP:endP)
        tempVal <- trunc(length(GraphRange) / nAxisPos)
        xPos <- c(0:nAxisPos) * tempVal
        noms <- ppm.metabolites.estimation[xPos + which(ppm == round(graphicalZonesBorders[[g]][1], 1))[1], 1]
        axis(1, at = xPos, labels = noms)
    }
}

invisible(dev.off())


## Ending
## ---------------------
cat("\nEnd of 'NMR annotation' Galaxy module call: ", as.character(Sys.time()), sep = "")
options(stringsAsFactors = strAsFacL)
rm(list = ls())
sink()
