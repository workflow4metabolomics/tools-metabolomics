#!/usr/bin/env Rscript
# metaMS_runGC.r version="3.0.0"
# created by Yann GUITTON and updated by Julien SAINT-VANNE
# use RI options + add try on plotUnknown add session Info
# use make.names in sampleMetadata to avoid issues with files names


# ----- LOG FILE -----
# log_file=file("log.txt", open = "wt")
# sink(log_file)
# sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n\n")

# Import the different functions and packages
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep = "/"))
}
source_local("lib_metams.r")

pkgs <- c("metaMS", "stringr", "batch", "CAMERA") # "batch" necessary for parseCommandArgs function
loadAndDisplayPackages(pkgs)

cat("\n\n")

modNamC <- "metaMS:runGC" ## module name
cat("\nStart of the '", modNamC, "' Galaxy module call: ", format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")


# ----- PROCESSING INFILE -----
cat("\n\n\tARGUMENTS PROCESSING INFO\n\n")
args <- parseCommandArgs(evaluate = FALSE) # interpretation of arguments given in command line as an R list of objects
# write.table(as.matrix(args), col.names=F, quote=F, sep='\t\t')
print(cbind(value = unlist(args)))

# ----- PROCESSING INFILE -----
cat("\n\n\tARGUMENTS PROCESSING INFO\n\n")

# Saving the specific parameters
# RI parameter
if (args$ri != "NULL") {
  RIarg <- read.table(args$ri)
  if (ncol(RIarg) < 2) RIarg <- read.table(args$ri, h = T, sep = ";")
  if (ncol(RIarg) < 2) RIarg <- read.table(args$ri, h = T, sep = "\t")
  if (ncol(RIarg) < 2) RIarg <- read.table(args$ri, h = T, sep = ",")
  if (ncol(RIarg) < 2) {
    error_message <- "Your RI file seems not well formatted. The column separators accepted are ; , and tabulation"
    print(error_message)
    stop(error_message)
  }
  # to do check real column names
  colnames(RIarg) <- c("rt", "RI")
} else {
  RIarg <- NULL
}

# RIshift parameter
if (args$rishift != "none") {
  RIshift <- args$rishift
} else {
  RIshift <- "none"
}

# Personal database parameter
if (args$db != "NULL") {
  DBgc <- args$db
} else {
  DBgc <- NULL
}

# settings process
if (args$settings == "default") {
  cat("Using default parameters")
  data(FEMsettings)
  if (args$rtrange[1] != "NULL") {
    rtrange <- args$rtrange
  } else {
    rtrange <- NULL
  }

  if (!is.null(DBgc)) {
    manual <- read.msp(DBgc)
    DBgc <- createSTDdbGC(stdInfo = NULL, settings = TSQXLS.GC, manualDB = manual)
  }

  # use RI instead of rt for time comparison vs DB
  if (RIshift != "none") {
    TSQXLS.GC@match2DB.timeComparison <- "RI"
    TSQXLS.GC@match2DB.RIdiff <- as.numeric(RIshift)
    TSQXLS.GC@betweenSamples.timeComparison <- "RI"
    TSQXLS.GC@betweenSamples.RIdiff <- as.numeric(RIshift)
  }
  nSlaves <- args$nSlaves
}

if (args$settings == "User_defined") {
  cat("Using user's parameters\n")
  rtdiffparam <- args$rtdiff
  minfeatparam <- args$minfeat
  simthreshparam <- args$simthreshold
  minclassfractionparam <- args$minclassfraction
  minclasssizeparam <- args$minclasssize

  if (args$rtrange != "NULL") {
    rtrange <- args$rtrange
    cat("rtrange= ", rtrange)
  } else {
    rtrange <- NULL
    cat("rtrange= ", rtrange)
  }

  nSlaves <- args$nSlaves
  # peakpicking parameters are not needed anymore xcms does it before runGC but keep it for protocolName and chrom type (2020-06-18 Yann Guitton)
  GALAXY.GC <-
    metaMSsettings(
      "protocolName" = "GALAXY.GC",
      "chrom" = "GC",
      PeakPicking = list(
        method = "matchedFilter",
        step = 0.5,
        steps = 2,
        mzdiff = .5,
        fwhm = 5,
        snthresh = 2,
        max = 500
      ),
      CAMERA = list(perfwhm = 1)
    )

  metaSetting(GALAXY.GC, "DBconstruction") <- list(
    minintens = 0.0,
    rttol = rtdiffparam,
    intensityMeasure = "maxo",
    DBthreshold = .80,
    minfeat = minfeatparam
  )
  metaSetting(GALAXY.GC, "match2DB") <- list(
    simthresh = simthreshparam,
    timeComparison = "rt",
    rtdiff = rtdiffparam,
    RIdiff = 5,
    minfeat = minfeatparam
  )

  # to used if contaminant filter

  # metaSetting(GALAXY.GC, "matchIrrelevants") <- list(
  # irrelevantClasses = c("Bleeding", "Plasticizers"),
  # timeComparison = "RI",
  # RIdiff = RIdiffparam,
  # rtdiff = rtdiffparam,
  # simthresh = simthreshparam)

  metaSetting(GALAXY.GC, "betweenSamples") <- list(
    min.class.fraction = minclassfractionparam,
    min.class.size = minclasssizeparam,
    timeComparison = "rt",
    rtdiff = rtdiffparam,
    RIdiff = 2,
    simthresh = simthreshparam
  )

  # ONLY use RI instead of rt for time comparison vs DB or samples
  if (RIshift != "none") {
    GALAXY.GC@match2DB.timeComparison <- "RI"
    GALAXY.GC@match2DB.RIdiff <- as.numeric(RIshift)
    GALAXY.GC@betweenSamples.timeComparison <- "RI"
    GALAXY.GC@betweenSamples.RIdiff <- as.numeric(RIshift)
  }

  if (!is.null(DBgc)) {
    manual <- read.msp(DBgc)
    DBgc <- createSTDdbGC(stdInfo = NULL, settings = GALAXY.GC, manualDB = manual)
  }
}


# ----- INFILE PROCESSING -----
cat("\n\n\n\tINFILE PROCESSING INFO\n\n")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
directory <- retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)


# ----- MAIN PROCESSING INFO -----
cat("\n\tMAIN PROCESSING INFO\n")

cat("\t\tCOMPUTE\n\n")

# runGC accept either a list of files a zip folder or an xset object from xcms.xcmsSet tool
# From xset is an .RData file necessary to use the xcmsSet object generated by xcms.xcmsSet given by previous tools
if (!is.null(args$singlefile_galaxyPath)) {
  cat("Loading datas from XCMS files...\n")
  load(args$singlefile_galaxyPath)

  # Transform XCMS object if needed
  if (!exists("xset")) {
    if (exists("xdata")) {
      xset <- getxcmsSetObject(xdata)
    } else {
      error_message <- "no xset and no xdata... Probably a problem"
      print(error_message)
      stop(error_message)
    }
  }
  # Verify that there are more than 1 file (can't run metaMS on only 1 file)
  if (length(rownames(xdata@phenoData)) < 2) {
    error_message <- "You need more than 1 file to be able to run metaMS"
    print(error_message)
    stop(error_message)
  }

  # xset from xcms.xcmsSet is not well formatted for metaMS this function do the formatting
  if (class(xset) == "xcmsSet") {
    if (length(xset@rt$raw) > 1) {
      # create an exceptable list of xset for metaMS
      xset.l <- vector("list", length(xset@rt$raw))
      for (i in 1:length(xset@rt$raw)) {
        xset.l[[i]] <- new("xcmsSet")
        xset.l[[i]]@peaks <- xset@peaks[which(xset@peaks[, "sample"] == i), ]
        df <- data.frame(class = xset@phenoData[i, ])
        rownames(df) <- rownames(xset@phenoData)[i]
        xset.l[[i]]@phenoData <- df
        xset.l[[i]]@rt$raw <- xset@rt$raw[[i]]
        xset.l[[i]]@rt$corrected <- xset@rt$corrected[[i]]
        xset.l[[i]]@filepaths <- xset@filepaths[i]
        xset.l[[i]]@profinfo <- xset@profinfo
      }
    } else {
      xset.l <- xset
    }

    # create sampleMetadata, get sampleMetadata and class
    sampleMetadata <- xset@phenoData
    colnames(sampleMetadata) <- c("sampleMetadata", "sample_group", "class")
    sampleMetadata <- sampleMetadata[, -2]
    row.names(sampleMetadata) <- NULL
    samples <- xset@filepaths
  } else {
    xset <- NULL
  }
  if (args$settings == "default") {
    settingslist <- TSQXLS.GC
    if (class(xset.l[[1]]) != "xsAnnotate") {
      cat("Process xsAnnotate with CAMERA package...\n")
      xsetCAM <- lapply(
        xset.l,
        function(x) {
          y <- xsAnnotate(x, sample = 1)
          capture.output(z <- groupFWHM(y, perfwhm = settingslist@CAMERA$perfwhm),
            file = NULL
          )
          z
        }
      )
    } else {
      xsetCAM <- xset.l
    }

    # default settings for GC from Wehrens et al
    cat("Process runGC with metaMS package...\n\n")
    print(str(TSQXLS.GC))
    resGC <- runGC(
      xset = xsetCAM, settings = TSQXLS.GC, rtrange = rtrange, DB = DBgc, removeArtefacts = TRUE,
      findUnknowns = TRUE, returnXset = TRUE, RIstandards = RIarg, nSlaves = nSlaves
    )
  } else {
    if (args$settings == "User_defined") {
      settingslist <- GALAXY.GC
      if (class(xset.l[[1]]) != "xsAnnotate") {
        cat("Process xsAnnotate with CAMERA package...\n")
        xsetCAM <- lapply(
          xset.l,
          function(x) {
            y <- xsAnnotate(x, sample = 1)
            capture.output(z <- groupFWHM(y, perfwhm = settingslist@CAMERA$perfwhm),
              file = NULL
            )
            z
          }
        )
      } else {
        xsetCAM <- xset.l
      }

      # user settings for GC
      cat("Process runGC with metaMS package...\n\n")
      print(str(GALAXY.GC))
      resGC <- runGC(
        xset = xsetCAM, settings = GALAXY.GC, rtrange = rtrange, DB = DBgc, removeArtefacts = TRUE,
        findUnknowns = TRUE, returnXset = TRUE, RIstandards = RIarg, nSlaves = nSlaves
      )
    } else {
      error_message <- "There is no xset"
      print(error_message)
      stop(error_message)
    }
  }
} else {
  error_message <- "No galaxy path entered"
  print(error_message)
  stop(error_message)
}


# ----- EXPORT -----
# peakTable ordered by rt
cat("\nGenerating peakTable file")
peaktable <- getCorrectFileName(resGC$PeakTable, sampleMetadata)
cat("\t.\t.")
write.table(peaktable, file = "peaktable.tsv", sep = "\t", row.names = FALSE)
cat("\t.\tOK")

# peakTable for PCA
# dataMatrix
cat("\nGenerating dataMatrix file")
dataMatrix <- cbind(Name = peaktable[, "Name"], peaktable[, (colnames(peaktable) %in% sampleMetadata[, 1])])
rownames(dataMatrix) <- NULL
cat("\t.\t.")
write.table(dataMatrix, file = "dataMatrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\t.\tOK")

# variableMetadata
cat("\nGenerating variableMetadata file")
variableMetadata <- peaktable[, !(colnames(peaktable) %in% sampleMetadata[, 1])]
rownames(variableMetadata) <- NULL
cat("\t.")
write.table(variableMetadata, file = "variableMetadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\t.\tOK")

# sampleMetadata
cat("\nGenerating sampleMetadata file")
cat("\t.\t.")
write.table(sampleMetadata, file = "sampleMetadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\t.\tOK")

# peak spectrum as MSP for DB search
cat("\nGenerating", length(resGC$PseudoSpectra), "peakspectra in peakspectra.msp file\n")
write.msp(resGC$PseudoSpectra, file = "peakspectra.msp", newFile = TRUE)

# saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("resGC", "xset", "singlefile", "zipfile", "DBgc")
save(list = objects2save[objects2save %in% ls()], file = "runGC.RData")

cat("\nEnd of '", modNamC, "' Galaxy module call: ", as.character(Sys.time()), "\n", sep = "")

# WARNING if user has CDF files (not yet good for plotting)
files <- paste("./", names(singlefile), sep = "")
if (MSnbase:::isCdfFile(files)) {
  warning_message <- "You have CDF files, for the moment you can't obtain plot after runGC! A new update will follow with the good correction"
  warning(warning_message)
  cat(paste("\n", "/!\\Warning/!\\", warning_message, sep = "\n"))
}
