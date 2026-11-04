# ==============================================================================
# Project: Workflow4Metabolomics / Workflow4Exposomic PARC project
# File: lib.r
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
retrieveRawfileInTheWorkingDir <- function(singlefile, zipfile, args, prefix = "") { # nolint
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
#' Merge chromatogram objects
#'
#' Extract data from MsExperiment object.
#'
#' @param xdata MsExperiment R object.
#'
#' @return a file list.
#'
#' @author Camille Trottier
getMsExperimentSampleFiles <- function(xdata) {
    sp <- MsExperiment::spectra(xdata)
    sampleIdx <- MsExperiment::spectraSampleIndex(xdata)
    origins <- Spectra::dataOrigin(sp)
    nSamples <- nrow(MsExperiment::sampleData(xdata))
    vapply(seq_len(nSamples), function(i) {
        f <- unique(origins[sampleIdx == i])
        if (length(f) != 1) {
            stop("\n\nERROR: Unable to determine a single raw data file for sample ", i, " of a 'MsExperiment' object (found ", length(f), "). Merging currently only supports one raw data file per sample.")
        }
        f
    }, character(1))
}

#------------------------------------------------------------------------
#' Merge multiple MsExperiment or XcmsExperiment objects
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
    rawFiles_merged <- NULL
    msExperimentSampleData_merged <- NULL
    for (image in args$images) {
        load(image)
        # Handle infiles
        if (!exists("singlefile")) singlefile <- NULL
        if (!exists("zipfile")) zipfile <- NULL
        rawFilePath <- retrieveRawfileInTheWorkingDir(singlefile, zipfile, args)
        # zipfile <- rawFilePath$zipfile
        singlefile <- rawFilePath$singlefile

        if (exists("raw_data")) xdata <- raw_data
        if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS 2.*")

        # XCMS 4.x - MsExperiment/XcmsExperiment objects only.
        if (is(xdata, "XCMSnExp") || is(xdata, "OnDiskMSnExp")) {
            stop("\n\nERROR: The RData contains a legacy '", class(xdata)[1], "' object. This function only supports the 'MsExperiment'/'XcmsExperiment' objects produced by xcms >= 4. Please reprocess your data with a recent version of xcms.")
        }
        if (!is(xdata, "MsExperiment")) {
            stop("\n\nERROR: Unsupported object of class '", paste(class(xdata), collapse = "/"), "'. Expected a 'MsExperiment' or 'XcmsExperiment' object.")
        }

        cat(sampleNamesList$sampleNamesOrigin, "\n")

        if (!exists("xdata_merged")) {
            xdata_merged <- xdata
            singlefile_merged <- singlefile
            md5sumList_merged <- md5sumList
            sampleNamesList_merged <- sampleNamesList
            chromTIC_merged <- chromTIC
            chromBPI_merged <- chromBPI
            chromTIC_adjusted_merged <- chromTIC_adjusted
            chromBPI_adjusted_merged <- chromBPI_adjusted
            if (!is(xdata, "XcmsExperiment")) {
                rawFiles_merged <- getMsExperimentSampleFiles(xdata)
                msExperimentSampleData_merged <- as.data.frame(MsExperiment::sampleData(xdata))
            }
        } else {
            if (!identical(class(xdata)[1], class(xdata_merged)[1])) {
                stop("\n\nERROR: All the RData to merge must contain the same type of object (currently mixing '", class(xdata_merged)[1], "' and '", class(xdata)[1], "').")
            }
            if (is(xdata, "XcmsExperiment")) {
                xdata_merged <- c(xdata_merged, xdata)
            } else {
                # The merged object itself is rebuilt once, after the
                # loop, with MsExperiment::readMsExperiment()
                rawFiles_merged <- c(rawFiles_merged, getMsExperimentSampleFiles(xdata))
                msExperimentSampleData_merged <- rbind(msExperimentSampleData_merged, as.data.frame(MsExperiment::sampleData(xdata)))
            }

            singlefile_merged <- c(singlefile_merged, singlefile)
            md5sumList_merged$origin <- rbind(md5sumList_merged$origin, md5sumList$origin)
            sampleNamesList_merged$sampleNamesOrigin <- c(sampleNamesList_merged$sampleNamesOrigin, sampleNamesList$sampleNamesOrigin)
            sampleNamesList_merged$sampleNamesMakeNames <- c(sampleNamesList_merged$sampleNamesMakeNames, sampleNamesList$sampleNamesMakeNames)
            chromTIC_merged <- mergeChrom(chromTIC_merged, chromTIC)
            chromBPI_merged <- mergeChrom(chromBPI_merged, chromBPI)
            chromTIC_adjusted_merged <- mergeChrom(chromTIC_adjusted_merged, chromTIC_adjusted)
            chromBPI_adjusted_merged <- mergeChrom(chromBPI_adjusted_merged, chromBPI_adjusted)
        }
    }
    rm(image)
    xdata <- xdata_merged
    rm(xdata_merged)
    singlefile <- singlefile_merged
    rm(singlefile_merged)
    md5sumList <- md5sumList_merged
    rm(md5sumList_merged)
    sampleNamesList <- sampleNamesList_merged
    rm(sampleNamesList_merged)

    # Plain MsExperiment objects could not be combined with c() : rebuild the full merged object
    # from the raw files and sample metadata
    if (is(xdata, "MsExperiment") && !is(xdata, "XcmsExperiment") && length(args$images) > 1) {
        xdata <- MsExperiment::readMsExperiment(spectraFiles = rawFiles_merged, sampleData = msExperimentSampleData_merged)
    }
    if (!is.null(args$sampleMetadata)) {
        cat("\tXSET METADATA SETTING...\n")
        sampleMetadataFile <- args$sampleMetadata
        sampleMetadata <- getDataFrameFromFile(sampleMetadataFile, header = FALSE)

        sd <- as.data.frame(MsExperiment::sampleData(xdata))
        sd$sample_group <- sampleMetadata$V2[match(sd$sample_name, sampleMetadata$V1)]

        if (any(is.na(sd$sample_group))) {
            sample_missing <- sd$sample_name[is.na(sd$sample_group)]
            error_message <- paste("Those samples are missing in your sampleMetadata:", paste(sample_missing, collapse = " "))
            print(error_message)
            stop(error_message)
        }

        sd_x <- MsExperiment::sampleData(xdata)
        sd_x$sample_group <- sd$sample_group
        MsExperiment::sampleData(xdata) <- sd_x
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

    return(list("xdata" = xdata, "md5sumList" = md5sumList, "sampleNamesList" = sampleNamesList, "chromTIC" = chromTIC, "chromBPI" = chromBPI, "chromTIC_adjusted" = chromTIC_adjusted, "chromBPI_adjusted" = chromBPI_adjusted))
}

#------------------------------------------------------------------------
#' Generate interactive chromatogram plots
#'
#' Creates Plotly visualizations of TIC or BPI chromatograms and exports them
#' as a self-contained HTML file.
#'
#' @param chrom Chromatogram object.
#' @param xdata MsExperiment object.
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
    p_sample <- plot_ly(df_all,
        x = ~rt, y = ~intensity, color = ~sample,
        type = "scatter", mode = "lines", line = list(width = 0.6),
        legendgroup = "group", legendgrouptitle = list(text = "Samples")
    )

    # Color by group
    p_group <- plot_ly(df_all,
        x = ~rt, y = ~intensity, color = ~group,
        type = "scatter", mode = "lines", line = list(width = 0.6),
        legendgroup = "sample", legendgrouptitle = list(text = "Groups")
    )

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
    # Initialisation
    sp <- spectra(xdata)
    files <- xdata@sampleData$spectraOrigin

    sp_origin <- normalizePath(Spectra::dataOrigin(sp), mustWork = FALSE)

    # For each sample file, the following actions are done
    for (fileIdx in seq_along(files)) {
        # Check if the file is in the CDF format
        if (!mzR:::netCDFIsFile(files[fileIdx])) {
            # If the column isn't exist, with add one filled with NA
            if (is.null(sampleMetadata$polarity)) sampleMetadata$polarity <- NA

            # Extract the polarity (a list of polarities)
            normFile <- normalizePath(files[fileIdx], mustWork = FALSE)
            pol <- Spectra::polarity(sp[sp_origin == normFile])
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


#------------------------------------------------------------------------
#' Generate chromatographic peak density plots for all features
#' @param xdata An XcmsExperiment object.
#' @param param Parameters for xcms::plotChromPeakDensity() function
#' @param mzdigit Number (integer) of digits for round() function
#'
#' @return A PDF file named `"plotChromPeakDensity.pdf"`
#'
#' @author Gildas Le Corguille
#' @author Camille Trottier

getPlotChromPeakDensity <- function(xdata, param = NULL, mzdigit = 4) {
    pdf(file = "plotChromPeakDensity.pdf", width = 16, height = 10)

    par(
        mfrow = c(3, 1), mar = c(5, 5, 3, 1),
        cex.main = 1.8, # title
        cex.lab = 1.6, # axes
        cex.axis = 1.1
    ) # ticks

    sample_groups <- MsExperiment::sampleData(xdata)$sample_group
    n_groups <- length(unique(MsExperiment::sampleData(xdata)$sample_group))

    if (n_groups <= 9) {
        group_colors <- brewer.pal(max(n_groups, 3), "Set1")[1:n_groups]
    } else {
        group_colors <- hcl.colors(n_groups, palette = "Dark 3")
    }
    names(group_colors) <- unique(sample_groups)
    col_per_samp <- unname(group_colors[as.character(sample_groups)])

    xlim <- c(min(featureDefinitions(xdata)$rtmin), max(featureDefinitions(xdata)$rtmax))

    for (i in seq_len(nrow(featureDefinitions(xdata)))) {
        # for (i in 1:10) { #for test only
        mzmin <- featureDefinitions(xdata)[i, ]$mzmin
        mzmax <- featureDefinitions(xdata)[i, ]$mzmax
        chr <- chromatogram(xdata, mz = c(mzmin, mzmax))
        # xlim cannot be used in plotChromPeakDensity in xcms v4
        # so filtering is used
        chr <- filterRt(chr, rt = xlim)

        peak_sample_idx <- chromPeaks(chr)[, "sample"]
        peak_colors <- col_per_samp[peak_sample_idx]
        title <- paste("Feature", i, "- m/z :", round(mzmin, mzdigit), "-", round(mzmax, mzdigit))

        plotChromPeakDensity(chr,
            param = param,
            col = col_per_samp,
            peakBg = adjustcolor(peak_colors, alpha.f = 0.3),
            peakCol = peak_colors,
            peakPch = 16,
            main = title
        )
        legend("topright", legend = names(group_colors), col = group_colors, cex = 1, lty = 1, lwd = 2)
    }

    dev.off()
}


#------------------------------------------------------------------------
#' If required, convert Retention Time (in seconds) to minutes
#'
#' @param variableMetadata variable data table obtained at group step
#' @param convertRTMinute logical converts retention times from seconds to minutes.
#'
#' @return variableMetadata with converted RTime in minutes
#'
#' @author Gildas Le Corguille
#' @author Camille Trottier
RTSecondToMinute <- function(variableMetadata, convertRTMinute) {
    if (convertRTMinute) {
        # converting the retention times (seconds) into minutes
        print("converting the retention times into minutes in the variableMetadata")
        variableMetadata[, "rt"] <- variableMetadata[, "rt"] / 60
        variableMetadata[, "rtmin"] <- variableMetadata[, "rtmin"] / 60
        variableMetadata[, "rtmax"] <- variableMetadata[, "rtmax"] / 60
    }
    return(variableMetadata)
}


#------------------------------------------------------------------------
#' This function format ions identifiers

#' @param variableMetadata variable data table obtained at group step
#' @param numDigitsRT
#' @param numDigitsMZ
#'
#' @return variableMetadata formated ion identifiers
#'
#' @author Gildas Le Corguille
#' @author Camille Trottier
formatIonIdentifiers <- function(variableMetadata, numDigitsRT = 0, numDigitsMZ = 0) {
    splitDeco <- strsplit(as.character(rownames(variableMetadata)), "_")
    idsDeco <- sapply(
        splitDeco,
        function(x) {
            deco <- unlist(x)[2]
            if (is.na(deco)) {
                return("")
            } else {
                return(paste0("_", deco))
            }
        }
    )
    namecustom <- make.unique(paste0("M", round(variableMetadata[, "mz"], numDigitsMZ), "T", round(variableMetadata[, "rt"], numDigitsRT), idsDeco))
    variableMetadata <- cbind(name = rownames(variableMetadata), namecustom = namecustom, variableMetadata[, !(colnames(variableMetadata) %in% c("name"))])

    return(variableMetadata)
}


#------------------------------------------------------------------------
#' This function convert the remain NA to 0 in the dataMatrix

#' @param dataMatrix data matrix obtained at the group step
#' @param naTOzero boolean
#'
#' @return dataMatrix coorected with 0 replacing NA values
#'
#' @author Gildas Le Corguille
#' @author Camille Trottier
naTOzeroDataMatrix <- function(dataMatrix, naTOzero) {
    if (naTOzero) {
        dataMatrix[is.na(dataMatrix)] <- 0
    }
    return(dataMatrix)
}


#------------------------------------------------------------------------
#' This function built W4M files for variableMetadata and dataMatrix from XcmsExperiment object

#' @param xdata An XcmsExperiment object.
#' @param intval intensity value to extract for each feature
#' @param convertRTMinute logical converts retention times from seconds to minutes.
#' @param numDigitsMZ number of decimal digits to keep for m/z
#' @param numDigitsRT number of decimal digits to keep for retention time
#' @param naTOzero logical replace NA values in the dataMatrix with 0
#' @param variableMetadataOutput Output file path for the variableMetadata table
#' @param dataMatrixOutput Output file path for the dataMatrix table
#' @param sampleNamesList sample name list
#'
#' @return variableMetadata table
#' @return dataMatrix table
#'
#' @author Gildas Le Corguille
#' @author Camille Trottier
getPeaklistW4M <- function(xdata, intval = "into", convertRTMinute = FALSE, numDigitsMZ = 4, numDigitsRT = 0, naTOzero = TRUE, variableMetadataOutput, dataMatrixOutput, sampleNamesList) {
    dataMatrix <- featureValues(xdata, method = "medret", value = intval)
    colnames(dataMatrix) <- make.names(tools::file_path_sans_ext(colnames(dataMatrix)))

    variableMetadata <- as.data.frame(featureDefinitions(xdata))
    colnames(variableMetadata)[1] <- "mz"
    colnames(variableMetadata)[4] <- "rt"

    feat_names <- rownames(featureDefinitions(xdata))
    dataMatrix <- cbind(name = feat_names, dataMatrix)
    variableMetadata <- RTSecondToMinute(variableMetadata, convertRTMinute)
    variableMetadata <- formatIonIdentifiers(variableMetadata, numDigitsRT = numDigitsRT, numDigitsMZ = numDigitsMZ)
    dataMatrix <- naTOzeroDataMatrix(dataMatrix, naTOzero)

    # FIX: issue when the vector at peakidx is too long and is written in a new line during the export
    variableMetadata[, "peakidx"] <- vapply(variableMetadata[, "peakidx"], FUN = paste, FUN.VALUE = character(1), collapse = ",")

    write.table(variableMetadata, file = variableMetadataOutput, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(dataMatrix, file = dataMatrixOutput, sep = "\t", quote = FALSE, row.names = FALSE)
}


#------------------------------------------------------------------------
#' Generates a PDF with plots of adjusted retention time deviation across samples from an XCMS object
#'
#' @param xdata XcmsExperiment object
#'
#' @return Write PDF file named `raw_vs_adjusted_rt.pdf`
#'
#' @author Gildas Le Corguille
#' @author Camille Trottier
getPlotAdjustedRtime <- function(xdata) {
    pdf(file = "raw_vs_adjusted_rt.pdf", width = 16, height = 12)

    sample_groups <- MsExperiment::sampleData(xdata)$sample_group
    n_groups <- length(unique(MsExperiment::sampleData(xdata)$sample_group))

    if (n_groups <= 9) {
        group_colors <- brewer.pal(max(n_groups, 3), "Set1")[1:n_groups]
    } else {
        group_colors <- hcl.colors(n_groups, palette = "Dark 3")
    }
    # Color by group
    if (length(group_colors) > 1) {
        names(group_colors) <- unique(sample_groups)
        col_per_samp <- unname(group_colors[as.character(sample_groups)])
        plotAdjustedRtime(xdata, col = col_per_samp)
        legend("topright", legend = names(group_colors), col = group_colors, cex = 0.8, lty = 1)
    }

    # Color by sample
    plotAdjustedRtime(xdata, col = rainbow(length(xdata@sampleData$sample_name)))
    legend("topright", legend = xdata@sampleData$sample_name, col = rainbow(length(xdata@sampleData$sample_name)), cex = 0.8, lty = 1)

    dev.off()
}
