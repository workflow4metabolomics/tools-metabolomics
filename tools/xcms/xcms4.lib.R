# ==============================================================================
# Project: Workflow4Metabolomics / Workflow4Exposomic PARC project
# File: xcms4.lib.r
#
# Description:
# Utility functions used by Galaxy XCMS tools
# This file provides helper functions for:
#   - command-line argument parsing
#   - package loading and session reporting
#   - raw file handling and validation
#   - metadata generation
#   - chromatogram visualization
#   - MsExperiment merging and compatibility utilities
#
# Authors:
#   - ABiMS Team
#   - LABERCA - PARC project founding
#   - Gildas Le Corguille
#   - Misharl Monsoor
#   - Camille Trottier
#
# License:
#   GPL-3.0-or-later
#
# ============================================================================== 

#' Parse command-line arguments and convert boolean strings
#' Solve an issue with batch if arguments are logical TRUE/FALSE
#' Wrapper around `batch::parseCommandArgs()` that converts character values
#' equal to `"TRUE"` or `"FALSE"` into logical values.
#'
#' @param ... Arguments passed to `batch::parseCommandArgs()`.
#'
#' @return A named list of parsed command-line arguments.
#'
#' @author Gildas Le Corguille
parseCommandArgs <- function(...) {
    args <- batch::parseCommandArgs(...)
    for (key in names(args)) {
        if (args[key] %in% c("TRUE", "FALSE")) {
            args[key] <- as.logical(args[key])
        }
    }
    return(args)
}

#------------------------------------------------------------------------
#' Load packages and display session information
#'
#' Loads the requested packages and prints R session information, including
#' versions of attached and loaded packages.
#'
#' @param pkgs Character vector of package names.
#'
#' @return No return value. Information is printed to the console.
#'
#' @author Gildas Le Corguille
loadAndDisplayPackages <- function(pkgs) {
    for (pkg in pkgs) suppressPackageStartupMessages(stopifnot(library(pkg, quietly = TRUE, logical.return = TRUE, character.only = TRUE)))

    sessioninfo <- sessionInfo()
    cat(sessioninfo$R.version$version.string, "\n")
    cat("Main packages:\n")
    for (pkg in names(sessioninfo$otherPkgs)) {
        cat(paste(pkg, packageVersion(pkg)), "\t")
    }
    cat("\n")
    cat("Other loaded packages:\n")
    for (pkg in names(sessioninfo$loadedOnly)) {
        cat(paste(pkg, packageVersion(pkg)), "\t")
    }
    cat("\n")
}

#------------------------------------------------------------------------
#' Read a tabular file with automatic delimiter detection
#'
#' Attempts to read a file using semicolon, tabulation, or comma separators.
#'
#' @param filename Path to the input file.
#' @param header Logical indicating whether the file contains column names.
#'
#' @return A data.frame.
#'
#' @author Gildas Le Corguille
getDataFrameFromFile <- function(filename, header = TRUE) {
    myDataFrame <- read.table(filename, header = header, sep = ";", stringsAsFactors = FALSE)
    if (ncol(myDataFrame) < 2) myDataFrame <- read.table(filename, header = header, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(myDataFrame) < 2) myDataFrame <- read.table(filename, header = header, sep = ",", stringsAsFactors = FALSE)
    if (ncol(myDataFrame) < 2) {
        error_message <- "Your tabular file seems not well formatted. The column separators accepted are ; , and tabulation"
        print(error_message)
        stop(error_message)
    }
    return(myDataFrame)
}

#------------------------------------------------------------------------
#' Compute MD5 checksums
#'
#' Computes MD5 hashes for a set of files to verify data integrity.
#'
#' @param files Character vector of file paths.
#'
#' @return A matrix containing MD5 checksums.
#'
#' @author Gildas Le Corguille
getMd5sum <- function(files) {
    cat("Compute md5 checksum...\n")
    library(tools)
    return(as.matrix(md5sum(files)))
}

#------------------------------------------------------------------------
#' Retrieve raw files into the working directory
#'
#' Imports raw data files from individual inputs or ZIP archives and prepares
#' them for downstream processing.
#'
#' @param singlefile Named list of input files.
#' @param zipfile Path to a ZIP archive.
#' @param args Parsed command-line arguments.
#' @param prefix Acquisition mode prefix.
#'
#' @return A list containing imported files and metadata.
#'
#' @author Gildas Le Corguille
retrieveRawfileInTheWorkingDir <- function(singlefile, zipfile, args, prefix = "") {
    if (!(prefix %in% c("", "Positive", "Negative", "MS1", "MS2"))) stop("prefix must be either '', 'Positive', 'Negative', 'MS1' or 'MS2'")

    # single - if the file are passed in the command arguments -> refresh singlefile
    if (!is.null(args[[paste0("singlefile_galaxyPath", prefix)]])) {
        singlefile_galaxyPaths <- unlist(strsplit(args[[paste0("singlefile_galaxyPath", prefix)]], "\\|"))
        singlefile_sampleNames <- unlist(strsplit(args[[paste0("singlefile_sampleName", prefix)]], "\\|"))

        singlefile <- NULL
        for (singlefile_galaxyPath_i in seq_len(length(singlefile_galaxyPaths))) {
            singlefile_galaxyPath <- singlefile_galaxyPaths[singlefile_galaxyPath_i]
            singlefile_sampleName <- singlefile_sampleNames[singlefile_galaxyPath_i]
            # In case, an url is used to import data within Galaxy
            singlefile_sampleName <- tail(unlist(strsplit(singlefile_sampleName, "/")), n = 1)
            singlefile[[singlefile_sampleName]] <- singlefile_galaxyPath
        }
    }
    # zipfile - if the file are passed in the command arguments -> refresh zipfile
    if (!is.null(args[[paste0("zipfile", prefix)]])) {
        zipfile <- args[[paste0("zipfile", prefix)]]
    }

    # single
    if (!is.null(singlefile) && (length("singlefile") > 0)) {
        files <- vector()
        for (singlefile_sampleName in names(singlefile)) {
            singlefile_galaxyPath <- singlefile[[singlefile_sampleName]]
            if (!file.exists(singlefile_galaxyPath)) {
                error_message <- paste("Cannot access the sample:", singlefile_sampleName, "located:", singlefile_galaxyPath, ". Please, contact your administrator ... if you have one!")
                print(error_message)
                stop(error_message)
            }

            if (!suppressWarnings(try(file.link(singlefile_galaxyPath, singlefile_sampleName), silent = TRUE))) {
                file.copy(singlefile_galaxyPath, singlefile_sampleName)
            }
            files <- c(files, singlefile_sampleName)
        }
    }
    # zipfile
    if (!is.null(zipfile) && (zipfile != "")) {
        if (!file.exists(zipfile)) {
            error_message <- paste("Cannot access the Zip file:", zipfile, ". Please, contact your administrator ... if you have one!")
            print(error_message)
            stop(error_message)
        }
        suppressWarnings(unzip(zipfile, unzip = "unzip"))

        # get the directory name
        suppressWarnings(filesInZip <- unzip(zipfile, list = TRUE))
        directories <- unique(unlist(lapply(strsplit(filesInZip$Name, "/"), function(x) x[1])))
        directories <- directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
        directory <- "."
        if (length(directories) == 1) directory <- directories

        cat("files_root_directory\t", directory, "\n")

        filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
        filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
        info <- file.info(directory)
        listed <- list.files(directory[info$isdir], pattern = filepattern, recursive = TRUE, full.names = TRUE)
        files <- c(directory[!info$isdir], listed)
        exists <- file.exists(files)
        files <- files[exists]
    }
    return(list(zipfile = zipfile, singlefile = singlefile, files = files))
}

#------------------------------------------------------------------------
#' Merge chromatogram objects
#'
#' Combines chromatogram intensity matrices into a single object.
#'
#' @param chrom_merged Existing merged chromatogram.
#' @param chrom Chromatogram to append.
#'
#' @return A merged chromatogram object.
#'
#' @author Gildas Le Corguille
mergeChrom <- function(chrom_merged, chrom) {
    if (is.null(chrom_merged)) {
        return(NULL)
    }
    chrom_merged@.Data <- cbind(chrom_merged@.Data, chrom@.Data)
    return(chrom_merged)
}

#------------------------------------------------------------------------
#' Merge multiple MsExperiment objects
#'
#' Loads and combines several serialized experiment objects and associated
#' chromatographic information.
#'
#' @param args Parsed command-line arguments.
#'
#' @return A list containing merged experiment data and chromatograms.
#'
#' @author Gildas Le Corguille
#' @author Camille Trottier
mergeXData <- function(args) {
    chromTIC <- NULL
    chromBPI <- NULL
    chromTIC_adjusted <- NULL
    chromBPI_adjusted <- NULL
    md5sumList <- NULL
    for (image in args$images) {
        load(image)
        # Handle infiles
        if (!exists("singlefile")) singlefile <- NULL
        if (!exists("zipfile")) zipfile <- NULL
        rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
        zipfile <- rawFilePath$zipfile
        singlefile <- rawFilePath$singlefile

        if (exists("raw_data")) xdata <- raw_data
        if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS 3.*")

        # cat(sampleNamesList$sampleNamesOrigin, "\n")

        if (!exists("xdata_merged")) {
            xdata_merged <- xdata
            singlefile_merged <- singlefile
            md5sumList_merged <- md5sumList
            # sampleNamesList_merged <- sampleNamesList
            chromTIC_merged <- chromTIC
            chromBPI_merged <- chromBPI
            chromTIC_adjusted_merged <- chromTIC_adjusted
            chromBPI_adjusted_merged <- chromBPI_adjusted
        } else {
            if (is(xdata, "MsExperiment")) {
                xdata_merged <- c(xdata_merged, xdata)
                
            # } else if (is(xdata, "OnDiskMSnExp")) {
            #     # xdata_merged <- xcms:::.concatenate_OnDiskMSnExp(xdata_merged, xdata)
            } else {
                stop("\n\nERROR: The RData wrong format. Please use MsExperiment RData")
            }

            singlefile_merged <- c(singlefile_merged, singlefile)
            md5sumList_merged$origin <- rbind(md5sumList_merged$origin, md5sumList$origin)
            # sampleNamesList_merged$sampleNamesOrigin <- c(sampleNamesList_merged$sampleNamesOrigin, sampleNamesList$sampleNamesOrigin)
            # sampleNamesList_merged$sampleNamesMakeNames <- c(sampleNamesList_merged$sampleNamesMakeNames, sampleNamesList$sampleNamesMakeNames)
            chromTIC_merged <- mergeChrom(chromTIC_merged, chromTIC)
            chromBPI_merged <- mergeChrom(chromBPI_merged, chromBPI)
            chromTIC_adjusted_merged <- mergeChrom(chromTIC_adjusted_merged, chromTIC_adjusted)
            chromBPI_adjusted_merged <- mergeChrom(chromBPI_adjusted_merged, chromBPI_adjusted)
        }
    }
    rm(image)
    #Rebuild a merge MsExperiment object
    all_spectra <- do.call(c, lapply(xdata_merged, spectra))
    all_sampleData <- do.call(rbind, lapply(xdata_merged, sampleData))
    xdata <- MsExperiment(spectra = all_spectra, sampleData = all_sampleData)
    #Don't forget to link spectra and sample datasets  
    xdata <- linkSampleData(
        xdata,
        with = "sampleData.spectraOrigin = spectra.dataOrigin")
    rm(xdata_merged)
    singlefile <- singlefile_merged
    rm(singlefile_merged)
    md5sumList <- md5sumList_merged
    rm(md5sumList_merged)
    # sampleNamesList <- sampleNamesList_merged
    # rm(sampleNamesList_merged)

    if (!is.null(args$sampleMetadata)) {
        cat("\tXSET METADATA SETTING...\n")
        sampleMetadataFile <- args$sampleMetadata
        sampleMetadata <- getDataFrameFromFile(sampleMetadataFile, header = FALSE)
        xdata@sampleData$sample_group <- sampleMetadata$V2[match(xdata@sampleData$sample_name, sampleMetadata$V1)]

        if (any(is.na(sampleData(xdata)$sample_group))) {
            sample_missing <- sampleData(xdata)$sample_name[is.na(sampleData(xdata)$sample_group)]
            error_message <- paste("Those samples are missing in your sampleMetadata:", paste(sample_missing, collapse = " "))
            print(error_message)
            stop(error_message)
        }
    }

    if (!is.null(chromTIC_merged)) {
        chromTIC <- chromTIC_merged
        chromTIC@sampleData <- xdata@sampleData
    }
    if (!is.null(chromBPI_merged)) {
        chromBPI <- chromBPI_merged
        chromBPI@sampleData <- xdata@sampleData
    }
    if (!is.null(chromTIC_adjusted_merged)) {
        chromTIC_adjusted <- chromTIC_adjusted_merged
        chromTIC_adjusted@sampleData <- xdata@sampleData
    }
    if (!is.null(chromBPI_adjusted_merged)) {
        chromBPI_adjusted <- chromBPI_adjusted_merged
        chromBPI_adjusted@sampleData <- xdata@sampleData
    }

    # return(list("xdata" = xdata, "singlefile" = singlefile, "md5sumList" = md5sumList, "sampleNamesList" = sampleNamesList, "chromTIC" = chromTIC, "chromBPI" = chromBPI, "chromTIC_adjusted" = chromTIC_adjusted, "chromBPI_adjusted" = chromBPI_adjusted))
    return(list("xdata" = xdata, "singlefile" = singlefile, "md5sumList" = md5sumList, "chromTIC" = chromTIC, "chromBPI" = chromBPI, "chromTIC_adjusted" = chromTIC_adjusted, "chromBPI_adjusted" = chromBPI_adjusted))
}

#------------------------------------------------------------------------
#' Generate interactive chromatogram plots
#'
#' Creates Plotly visualizations of TIC or BPI chromatograms and exports them
#' as a self-contained HTML file.
#'
#' @param chrom Chromatogram object.
#' @param xdata Experiment object.
#' @param htmlFile Output HTML file.
#' @param aggregationFun Aggregation function used to generate chromatograms.
#'
#' @return No return value. An HTML file is written to disk.
#'
#' @author Camille Trottier
getPlotChromHTML <- function(chrom, xdata, htmlFile = "Chromatogram.html", aggregationFun = "max") {
    if (aggregationFun == "sum") {
        type <- "Total Ion Chromatograms"
    } else {
        type <- "Base Peak Intensity Chromatograms"
    }

    adjusted <- "Raw"
    if (hasAdjustedRtime(xdata)) {
        adjusted <- "Adjusted"
    }

    main <- paste(type, ":", adjusted, "data")

    # pdf(pdfname, width = 16, height = 10)

    # Color by sample
    plots <- lapply(seq_len(ncol(chrom)), function(i) {
        chr_i <- chrom[1, i]
        data.frame(
            rt = rtime(chr_i),
            intensity = intensity(chr_i),
            sample = xdata@sampleData$sample_name[i],
            group = xdata@sampleData$sample_group[i]
            )
    })

    df_all <- do.call(rbind, plots)
    p_sample <- plot_ly(df_all, x = ~rt, y = ~intensity, color = ~sample, 
                        type = "scatter",  mode = "lines", line = list(width = 0.6), 
                        legendgroup = "group", legendgrouptitle = list(text = "Samples")) 

    # Color by group
    p_group <- plot_ly(df_all, x = ~rt, y = ~intensity, color = ~group, 
                        type = "scatter",  mode = "lines", line = list(width = 0.6), 
                        legendgroup = "sample", legendgrouptitle = list(text = "Groups")) 

    # Combine plots
    combined_plot <- subplot(
        p_group,
        p_sample,
        nrows = 2,
        shareX = TRUE,
        titleY = TRUE
    ) %>%
    layout(legend = list(x = 1.02, y = 1))
    
    saveWidget(combined_plot, htmlFile, selfcontained = TRUE) 
}

#------------------------------------------------------------------------
#' Generate sample metadata
#'
#' Extracts sample information and polarity information from an experiment and
#' writes the metadata table to disk.
#'
#' @param xdata An MsExperiment object.
#' @param sampleMetadataOutput Output TSV file.
#'
#' @return A list containing original and sanitized sample names.
#'
#' @author Misharl Monsoor
#' @author Gildas Le Corguille
#' @author Camille Trottier
getSampleMetadata <- function(xdata = NULL, sampleMetadataOutput = "sampleMetadata.tsv") {
    cat("Creating the sampleMetadata file...\n")

    # Create the sampleMetada dataframe
    sampleMetadata <- xdata@sampleData
    rownames(sampleMetadata) <- NULL
    colnames(sampleMetadata) <- c("sample_name", "sample_group", "polarity")

    sampleNamesOrigin <- sampleMetadata$sample_name
    sampleNamesMakeNames <- make.names(sampleNamesOrigin)

    if (any(duplicated(sampleNamesMakeNames))) {
        write("\n\nERROR: Usually, R has trouble to deal with special characters in its column names, so it rename them using make.names().\nIn your case, at least two columns after the renaming obtain the same name, thus XCMS will collapse those columns per name.", stderr())
        for (sampleName in sampleNamesOrigin) {
            write(paste(sampleName, "\t->\t", make.names(sampleName)), stderr())
        }
        stop("\n\nERROR: One or more of your files will not be imported. It may due to bad characters in their filenames.")
    }

    if (!all(sampleNamesOrigin == sampleNamesMakeNames)) {
        cat("\n\nWARNING: Usually, R has trouble to deal with special characters in its column names, so it rename them using make.names()\nIn your case, one or more sample names will be renamed in the sampleMetadata and dataMatrix files:\n")
        for (sampleName in sampleNamesOrigin) {
            cat(paste(sampleName, "\t->\t", make.names(sampleName), "\n"))
        }
    }

    sampleMetadata$sample_name <- sampleNamesMakeNames
    #Initialisation 
    sp <- spectra(xdata)
    files <- xdata@sampleData$spectraOrigin

    # For each sample file, the following actions are done
    for (fileIdx in seq_along(files)) {
        # Check if the file is in the CDF format
        if (!mzR:::netCDFIsFile(files[fileIdx])) {
            # If the column isn't exist, with add one filled with NA
            if (is.null(sampleMetadata$polarity)) sampleMetadata$polarity <- NA

            # Extract the polarity (a list of polarities)
            pol <- polarity(sp[dataOrigin(sp) == files[fileIdx]])

            uniq_list <- unique(pol[!is.na(pol)])

            # Verify if all the scans have the same polarity
            sampleMetadata$polarity[fileIdx] <-
                if (length(uniq_list) > 1) {
                    pol <- "mixed"
                } else {
                    pol <- as.character(uniq_list)
                }
                # Set the polarity attribute
            sampleMetadata$polarity[fileIdx] <- pol
        }
    }

    write.table(sampleMetadata, sep = "\t", quote = FALSE, row.names = FALSE, file = sampleMetadataOutput)

    return(list("sampleNamesOrigin" = sampleNamesOrigin, "sampleNamesMakeNames" = sampleNamesMakeNames))
}

#------------------------------------------------------------------------
#' Retrieve a legacy xcmsSet object
#'
#' Converts modern XCMS objects into the legacy `xcmsSet` format when needed.
#'
#' @param xobject An XCMS object.
#'
#' @return An `xcmsSet` object.
#'
#' @author Gildas Le Corguille
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    # if (class(xobject) == "xcmsSet") {
    if (inherits(xobject, "xcmsSet")) {
        return(xobject)
    }
    # XCMS >= 3.x
    # if (class(xobject) == "XCMSnExp") {
    if (inherits(xobject, "XCMSnExp")) {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, "xcmsSet"))
        if (!is.null(xset@phenoData$sample_group)) {
            sampclass(xset) <- xset@phenoData$sample_group
        } else {
            sampclass(xset) <- "."
        }
        return(xset)
    }
}
