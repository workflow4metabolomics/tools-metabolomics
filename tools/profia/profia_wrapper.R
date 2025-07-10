#!/usr/bin/env Rscript

library(batch) ## parseCommandArgs

argVc <- unlist(parseCommandArgs(evaluate=FALSE))

##------------------------------
## Initializing
##------------------------------


## libraries
##----------

suppressMessages(library(proFIA))


## constants
##----------

modNamC <- "proFIA" ## module name


## log file
##---------

sink(argVc["information"])

cat("\nStart of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")


## arguments
##----------


if("zipfile" %in% names(argVc)) {
    
    zipfile <- argVc["zipfile"]

    ## We unzip automatically the raw files from the zip file
    
    if(exists("zipfile") && (zipfile!="")) {
        if(!file.exists(zipfile)){
            error_message=paste("Cannot access the Zip file:", zipfile)
            print(error_message)
            stop(error_message)
        }

        ## unzip

        suppressWarnings(unzip(zipfile, unzip="unzip"))
        
        ## get the directory name
        
        filesInZip=unzip(zipfile, list=T);
        directories=unique(unlist(lapply(strsplit(filesInZip$Name,"/"), function(x) x[1])));
        directories=directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
        directory = "."
        if (length(directories) == 1) directory = directories
        
        cat("files_root_directory\t",directory,"\n")
        
    }

} else if ("library" %in% names(argVc)) {

    directory <- argVc["library"]
    
    if(!file.exists(directory)) {
        
        error_message=paste("Cannot access the directory:", directory,". Please check that the directory really exists.")
        print(error_message)
        stop(error_message)
        
    }
    
} else {

    error_message <- "No zipfile nor input library available"
    print(error_message)
    stop(error_message)

}

##------------------------------
## Computations
##------------------------------


optWrnN <- options()$warn
options(warn = -1)

stpI <- 1

cat("\n", stpI, ") Peak detection step ('proFIAset'):\n", sep = "")

if("sizeMinN" %in% names(argVc) && argVc["sizeMinN"] != "none") {
    if("scanMaxI" %in% names(argVc) && argVc["scanMaxI"] != "none") {
        fiaset <- proFIAset(directory,
                            ppm = as.numeric(argVc["ppmN"]),
                            dmz = as.numeric(argVc["dmzN"]),
                            bandCoverage = as.numeric(ifelse("bandCoverageN" %in% names(argVc), argVc["bandCoverageN"], "0.3")),
                            sizeMin = as.numeric(argVc["sizeMinN"]),
                            scanmin = as.numeric(ifelse("scanMinI" %in% names(argVc), argVc["scanMinI"], "1")),
                            scanmax = as.numeric(argVc["scanMaxI"]),
                            parallel = TRUE)
    } else {
        fiaset <- proFIAset(directory,
                            ppm = as.numeric(argVc["ppmN"]),
                            dmz = as.numeric(argVc["dmzN"]),
                            bandCoverage = as.numeric(ifelse("bandCoverageN" %in% names(argVc), argVc["bandCoverageN"], "0.3")),
                            sizeMin = as.numeric(argVc["sizeMinN"]),
                            scanmin = as.numeric(ifelse("scanMinI" %in% names(argVc), argVc["scanMinI"], "1")),
                            parallel = TRUE)
    }
} else {
    if("scanMaxI" %in% names(argVc) && argVc["scanMaxI"] != "none") {
        fiaset <- proFIAset(directory,
                            ppm = as.numeric(argVc["ppmN"]),
                            dmz = as.numeric(argVc["dmzN"]),
                            bandCoverage = as.numeric(ifelse("bandCoverageN" %in% names(argVc), argVc["bandCoverageN"], "0.3")),
                            scanmin = as.numeric(ifelse("scanMinI" %in% names(argVc), argVc["scanMinI"], "1")),
                            scanmax = as.numeric(argVc["scanMaxI"]),
                            parallel = TRUE)
    } else {
        fiaset <- proFIAset(directory,
                            ppm = as.numeric(argVc["ppmN"]),
                            dmz = as.numeric(argVc["dmzN"]),
                            bandCoverage = as.numeric(ifelse("bandCoverageN" %in% names(argVc), argVc["bandCoverageN"], "0.3")),
                            scanmin = as.numeric(ifelse("scanMinI" %in% names(argVc), argVc["scanMinI"], "1")),
                            parallel = TRUE)
    }
}

stpI <- stpI + 1

cat("\n", stpI, ") Peak alignment ('group.FIA'):\n", sep = "")

fiaset <- group.FIA(fiaset,
                    ppmGroup = as.numeric(argVc["ppmGroupN"]),
                    dmz = as.numeric(argVc["dmzGroupN"]),
                    fracGroup = as.numeric(argVc["fracGroupN"]))

stpI <- stpI + 1

cat("\n", stpI, ") Creating the peak table ('makeDataMatrix'):\n", sep = "")

fiaset <- makeDataMatrix(fiaset,
                         maxo = FALSE)

stpI <- stpI + 1


if(argVc["imputeC"] == "randomForest") {

    cat("\n", stpI, ") Imputing missing values ('impute.randomForest'):\n", sep = "")

    set.seed(123)

    fiaset <- impute.randomForest(fiaset,
                                  parallel = FALSE)

    set.seed(NULL)

    stpI <- stpI + 1
}
###TODO Fix the error message to set 'parallel' to TRUE: You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.
###TODO add the impute.KNN_TN method


options(warn = optWrnN)


##------------------------------
## Ending
##------------------------------


## Plotting
##---------

cat("\n", stpI, ") Plotting ('plot'):\n", sep = "")

pdf(argVc["figure"])

plot(fiaset)

dev.off()

stpI <- stpI + 1

## Printing
##---------

cat("\n", stpI, ") Printing ('show'):\n", sep = "")

fiaset

stpI <- stpI + 1

## Exporting
##----------

cat("\n", stpI, ") Exporting ('exportDataMatrix', 'exportSampleMetadata', 'exportVariableMetadata'):\n", sep = "")

datMN <- exportDataMatrix(fiaset)
samDF <- exportSampleMetadata(fiaset)
varDF <- exportVariableMetadata(fiaset)

if(nrow(datMN) == nrow(samDF) && ncol(datMN) == nrow(varDF)) {
    datDF <- as.data.frame(t(datMN))
} else {
    datDF <- as.data.frame(datMN)
}
rownames(varDF) <- rownames(datDF)

datDF <- cbind.data.frame(dataMatrix = rownames(datDF),
                          datDF)
write.table(datDF,
            file = argVc["dataMatrix_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

samDF <- cbind.data.frame(sampleMetadata = samDF[, "sampleID"],
                          class = samDF[, "class"])
write.table(samDF,
            file = argVc["sampleMetadata_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                          varDF)
write.table(varDF,
            file = argVc["variableMetadata_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")


## Closing
##--------

cat("\nEnd of '", modNamC, "' Galaxy module call: ",
    as.character(Sys.time()), "\n", sep = "")

cat("\n\n\n============================================================================")
cat("\nAdditional information about the call:\n")
cat("\n1) Parameters:\n")
print(cbind(value = argVc))

cat("\n2) Session Info:\n")
sessioninfo <- sessionInfo()
cat(sessioninfo$R.version$version.string,"\n")
cat("Main packages:\n")
for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
cat("Other loaded packages:\n")
for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")

cat("============================================================================\n")

sink()

rm(list = ls())
