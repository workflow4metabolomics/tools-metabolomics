# lib_metams.r version 2.1.1
# R function for metaMS runGC under W4M
# author Yann GUITTON CNRS IRISA/LINA Idealg project 2014-2015
# author Yann GUITTON Oniris Laberca 2015-2017


# @author G. Le Corguille
# This function will
# - load the packages
# - display the sessionInfo
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

# This function list the compatible files within the directory as xcms did
# @author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
getMSFiles <- function(directory) {
  filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
  filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
  info <- file.info(directory)
  listed <- list.files(directory[info$isdir], pattern = filepattern, recursive = TRUE, full.names = TRUE)
  files <- c(directory[!info$isdir], listed)
  exists <- file.exists(files)
  files <- files[exists]
  return(files)
}

# This function retrieve a xset like object
# @author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
  # XCMS 1.x
  if (class(xobject) == "xcmsSet") {
    return(xobject)
  }
  # XCMS 3.x
  if (class(xobject) == "XCMSnExp") {
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

# J.Saint-Vanne
# Function to correct the file names which can be written like "..alg8.mzData" and we just want "alg8"
getCorrectFileName <- function(peaktable, sampleMetadata) {
  # Try to start for the first sample, avoid description of line with colnamesdontwant
  i <- 1
  while (!(sampleMetadata[1, 1] %in% strsplit(colnames(peaktable)[i], "\\.")[[1]])) {
    if (i < length(peaktable)) {
      i <- i + 1
    } else {
      break
    }
  }
  # i now correspond to the first column with a sample
  for (j in seq_len(nrow(sampleMetadata))) {
    col <- j + i - 1 # minus 1 cause i is the good column to start and j start at 1
    if (col <= length(colnames(peaktable))) {
      newname <- gsub("(^.*)(\\..*$)", "\\1", colnames(peaktable)[col])
      if (newname != sampleMetadata[j, 1]) {
        # Correction for 2 points starting the name (I don't know why they are here...)
        if (".." == gsub("(^\\.+)(.*)", "\\1", newname)) {
          newname <- sub("(^\\.+)(.*)", "\\2", newname)
        }
      }
      colnames(peaktable)[col] <- newname
    }
  }
  return(peaktable)
}


# @author Gildas Le Corguille lecorguille@sb-roscoff.fr
# This function get the raw file path from the arguments
getRawfilePathFromArguments <- function(singlefile, zipfile, listArguments) {
  if (!is.null(listArguments[["zipfile"]])) zipfile <- listArguments[["zipfile"]]
  if (!is.null(listArguments[["zipfilePositive"]])) zipfile <- listArguments[["zipfilePositive"]]
  if (!is.null(listArguments[["zipfileNegative"]])) zipfile <- listArguments[["zipfileNegative"]]

  if (!is.null(listArguments[["singlefile_galaxyPath"]])) {
    singlefile_galaxyPaths <- listArguments[["singlefile_galaxyPath"]]
    singlefile_sampleNames <- listArguments[["singlefile_sampleName"]]
  }
  if (!is.null(listArguments[["singlefile_galaxyPathPositive"]])) {
    singlefile_galaxyPaths <- listArguments[["singlefile_galaxyPathPositive"]]
    singlefile_sampleNames <- listArguments[["singlefile_sampleNamePositive"]]
  }
  if (!is.null(listArguments[["singlefile_galaxyPathNegative"]])) {
    singlefile_galaxyPaths <- listArguments[["singlefile_galaxyPathNegative"]]
    singlefile_sampleNames <- listArguments[["singlefile_sampleNameNegative"]]
  }
  if (exists("singlefile_galaxyPaths")) {
    singlefile_galaxyPaths <- unlist(strsplit(singlefile_galaxyPaths, ","))
    singlefile_sampleNames <- unlist(strsplit(singlefile_sampleNames, ","))

    singlefile <- NULL
    for (singlefile_galaxyPath_i in seq_along(singlefile_galaxyPaths)) {
      singlefile_galaxyPath <- singlefile_galaxyPaths[singlefile_galaxyPath_i]
      singlefile_sampleName <- singlefile_sampleNames[singlefile_galaxyPath_i]
      singlefile[[singlefile_sampleName]] <- singlefile_galaxyPath
    }
  }
  for (argument in c("zipfile", "zipfilePositive", "zipfileNegative", "singlefile_galaxyPath", "singlefile_sampleName", "singlefile_galaxyPathPositive", "singlefile_sampleNamePositive", "singlefile_galaxyPathNegative", "singlefile_sampleNameNegative")) {
    listArguments[[argument]] <- NULL
  }
  return(list(zipfile = zipfile, singlefile = singlefile, listArguments = listArguments))
}


# @author Gildas Le Corguille lecorguille@sb-roscoff.fr
# This function retrieve the raw file in the working directory
#   - if zipfile: unzip the file with its directory tree
#   - if singlefiles: set symlink with the good filename
retrieveRawfile <- function(singlefile, zipfile) {
  if (!is.null(singlefile) && (length("singlefile") > 0)) {
    for (singlefile_sampleName in names(singlefile)) {
      singlefile_galaxyPath <- singlefile[[singlefile_sampleName]]
      if (!file.exists(singlefile_galaxyPath)) {
        error_message <- paste("Cannot access the sample:", singlefile_sampleName, "located:", singlefile_galaxyPath, ". Please, contact your administrator ... if you have one!")
        print(error_message)
        stop(error_message)
      }
      file.symlink(singlefile_galaxyPath, singlefile_sampleName)
    }
    directory <- "."
  }
  if (!is.null(zipfile) && (zipfile != "")) {
    if (!file.exists(zipfile)) {
      error_message <- paste("Cannot access the Zip file:", zipfile, ". Please, contact your administrator ... if you have one!")
      print(error_message)
      stop(error_message)
    }

    # unzip
    suppressWarnings(unzip(zipfile, unzip = "unzip"))

    # get the directory name
    filesInZip <- unzip(zipfile, list = TRUE)
    directories <- unique(unlist(lapply(strsplit(filesInZip$Name, "/"), function(x) x[1])))
    directories <- directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
    directory <- "."
    if (length(directories) == 1) directory <- directories

    cat("files_root_directory\t", directory, "\n")
  }
  return(directory)
}

## ADDITIONS FROM Y. Guitton
getBPC <- function(file, rtcor = NULL, ...) {
  object <- xcmsRaw(file)
  sel <- profRange(object, ...)
  cbind(if (is.null(rtcor)) object@scantime[sel$scanidx] else rtcor, xcms:::colMax(object@env$profile[sel$massidx, sel$scanidx, drop = FALSE]))
}

getBPC2s <- function(files, xset = NULL, pdfname = "BPCs.pdf", rt = c("raw", "corrected"), scanrange = NULL) {
  require(xcms)

  # create sampleMetadata, get sampleMetadata and class
  if (!is.null(xset)) {
    # When files come from XCMS3 directly before metaMS
    sampleMetadata <- xset@phenoData
  } else {
    # When files come from a zip directory with raw files before metaMS
    sampleMetadata <- xcms:::phenoDataFromPaths(files)
  }
  class <- unique(sampleMetadata[, "class"]) # create phenoData like table
  classnames <- vector("list", length(class))
  for (i in seq_along(class)) {
    classnames[[i]] <- which(sampleMetadata[, "class"] == class[i])
  }

  N <- dim(sampleMetadata)[1]
  BPC <- vector("list", N)

  for (j in 1:N) {
    BPC[[j]] <- getBPC(files[j])
    # good for raw
    # seems strange for corrected
    # errors if scanrange used in xcmsSetgeneration
    if (!is.null(xcmsSet) && rt == "corrected") {
      rtcor <- xcmsSet@rt$corrected[[j]]
    } else {
      rtcor <- NULL
    }
    BPC[[j]] <- getBPC(files[j], rtcor = rtcor)
  }

  pdf(pdfname, w = 16, h = 10)
  cols <- rainbow(N)
  lty <- 1:N
  pch <- 1:N
  # search for max x and max y in BPCs

  xlim <- range(sapply(BPC, function(x) range(x[, 1])))
  ylim <- range(sapply(BPC, function(x) range(x[, 2])))

  ylim <- c(-ylim[2], ylim[2])

  ## plot start
  if (length(class) > 2) {
    for (k in 1:(length(class) - 1)) {
      for (l in (k + 1):length(class)) {
        cat(paste(class[k], "vs", class[l], sep = " ", "\n"))
        plot(0, 0, type = "n", xlim = xlim / 60, ylim = ylim, main = paste("Base Peak Chromatograms \n", "BPCs_", class[k], " vs ", class[l], sep = ""), xlab = "Retention Time (min)", ylab = "BPC")
        colvect <- NULL
        for (j in seq_along(classnames[[k]])) {
          bpc <- BPC[[classnames[[k]][j]]]

          points(bpc[, 1] / 60, bpc[, 2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
          colvect <- append(colvect, cols[classnames[[k]][j]])
        }
        for (j in seq_along(classnames[[l]])) {

          bpc <- BPC[[classnames[[l]][j]]]
          points(bpc[, 1] / 60, -bpc[, 2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type = "l")
          colvect <- append(colvect, cols[classnames[[l]][j]])
        }
        legend("topright", paste(gsub("(^.+)\\..*$", "\\1", basename(files[c(classnames[[k]], classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
      }
    }
  } # end if length >2

  if (length(class) == 2) {
    k <- 1
    l <- 2
    colvect <- NULL
    plot(0, 0, type = "n", xlim = xlim / 60, ylim = ylim, main = paste("Base Peak Chromatograms \n", "BPCs_", class[k], "vs", class[l], sep = ""), xlab = "Retention Time (min)", ylab = "BPC")

    for (j in seq_along(classnames[[k]])) {
      bpc <- BPC[[classnames[[k]][j]]]

      points(bpc[, 1] / 60, bpc[, 2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
      colvect <- append(colvect, cols[classnames[[k]][j]])
    }
    for (j in seq_along(classnames[[l]])) {

      bpc <- BPC[[classnames[[l]][j]]]
      points(bpc[, 1] / 60, -bpc[, 2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type = "l")
      colvect <- append(colvect, cols[classnames[[l]][j]])
    }
    legend("topright", paste(gsub("(^.+)\\..*$", "\\1", basename(files[c(classnames[[k]], classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
  } # end length ==2

  if (length(class) == 1) {
    k <- 1

    ylim <- range(sapply(BPC, function(x) range(x[, 2])))

    colvect <- NULL
    plot(0, 0, type = "n", xlim = xlim / 60, ylim = ylim, main = paste("Base Peak Chromatograms \n", "BPCs_", class[k], sep = ""), xlab = "Retention Time (min)", ylab = "BPC")

    for (j in seq_along(classnames[[k]])) {
      bpc <- BPC[[classnames[[k]][j]]]

      points(bpc[, 1] / 60, bpc[, 2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
      colvect <- append(colvect, cols[classnames[[k]][j]])
    }
    legend("topright", paste(gsub("(^.+)\\..*$", "\\1", basename(files[c(classnames[[k]])]))), col = colvect, lty = lty, pch = pch)
  } # end length ==1
  dev.off()
}

getTIC <- function(file, rtcor = NULL) {
  object <- xcmsRaw(file)
  cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object, mzrange = range(object@env$mz))$intensity)
}

##  overlay TIC from all files in current folder or from xcmsSet, create pdf
getTIC2s <- function(files, xset = NULL, pdfname = "TICs.pdf", rt = c("raw", "corrected")) {
  require(xcms)

  # create sampleMetadata, get sampleMetadata and class
  if (!is.null(xset)) {
    # When files come from XCMS3 before metaMS treatment
    sampleMetadata <- xset@phenoData
  } else {
    # When files come from a zip directory with raw files before metaMS
    sampleMetadata <- xcms:::phenoDataFromPaths(files)
  }
  class <- as.vector(levels(sampleMetadata[, "class"])) # create phenoData like table
  classnames <- vector("list", length(class))
  for (i in seq_along(class)) {
    classnames[[i]] <- which(sampleMetadata[, "class"] == class[i])
  }

  N <- dim(sampleMetadata)[1]
  TIC <- vector("list", N)

  for (i in 1:N) {
    if (!is.null(xcmsSet) && rt == "corrected") {
      rtcor <- xcmsSet@rt$corrected[[i]]
    } else {
      rtcor <- NULL
    }
    TIC[[i]] <- getTIC(files[i], rtcor = rtcor)
  }

  pdf(pdfname, w = 16, h = 10)
  cols <- rainbow(N)
  lty <- 1:N
  pch <- 1:N
  # search for max x and max y in TICs
  xlim <- range(sapply(TIC, function(x) range(x[, 1])))
  ylim <- range(sapply(TIC, function(x) range(x[, 2])))
  ylim <- c(-ylim[2], ylim[2])

  ## plot start
  if (length(class) > 2) {
    for (k in 1:(length(class) - 1)) {
      for (l in (k + 1):length(class)) {
        cat(paste(class[k], "vs", class[l], "\n", sep = " "))
        plot(0, 0, type = "n", xlim = xlim / 60, ylim = ylim, main = paste("Total Ion Chromatograms \n", "TICs_", class[k], " vs ", class[l], sep = ""), xlab = "Retention Time (min)", ylab = "TIC")
        colvect <- NULL
        for (j in seq_along(classnames[[k]])) {
          tic <- TIC[[classnames[[k]][j]]]

          points(tic[, 1] / 60, tic[, 2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
          colvect <- append(colvect, cols[classnames[[k]][j]])
        }
        for (j in seq_along(classnames[[l]])) {

          tic <- TIC[[classnames[[l]][j]]]
          points(tic[, 1] / 60, -tic[, 2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type = "l")
          colvect <- append(colvect, cols[classnames[[l]][j]])
        }
        legend("topright", paste(gsub("(^.+)\\..*$", "\\1", basename(files[c(classnames[[k]], classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
      }
    }
  } # end if length >2

  if (length(class) == 2) {
    k <- 1
    l <- 2
    plot(0, 0, type = "n", xlim = xlim / 60, ylim = ylim, main = paste("Total Ion Chromatograms \n", "TICs_", class[k], "vs", class[l], sep = ""), xlab = "Retention Time (min)", ylab = "TIC")
    colvect <- NULL
    for (j in seq_along(classnames[[k]])) {
      tic <- TIC[[classnames[[k]][j]]]

      points(tic[, 1] / 60, tic[, 2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
      colvect <- append(colvect, cols[classnames[[k]][j]])
    }
    for (j in seq_along(classnames[[l]])) {

      tic <- TIC[[classnames[[l]][j]]]
      points(tic[, 1] / 60, -tic[, 2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type = "l")
      colvect <- append(colvect, cols[classnames[[l]][j]])
    }
    legend("topright", paste(gsub("(^.+)\\..*$", "\\1", basename(files[c(classnames[[k]], classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
  } # end length ==2

  if (length(class) == 1) {
    k <- 1
    ylim <- range(sapply(TIC, function(x) range(x[, 2])))
    plot(0, 0, type = "n", xlim = xlim / 60, ylim = ylim, main = paste("Total Ion Chromatograms \n", "TICs_", class[k], sep = ""), xlab = "Retention Time (min)", ylab = "TIC")
    colvect <- NULL
    for (j in seq_along(classnames[[k]])) {
      tic <- TIC[[classnames[[k]][j]]]

      points(tic[, 1] / 60, tic[, 2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
      colvect <- append(colvect, cols[classnames[[k]][j]])
    }
    legend("topright", paste(gsub("(^.+)\\..*$", "\\1", basename(files[c(classnames[[k]])]))), col = colvect, lty = lty, pch = pch)
  } # end length ==1
  dev.off()
}


# Update by J.Saint-Vanne
## addition for quality control of peak picking
# metaMS EIC and pspectra plotting option
# version 20190520
# only for Galaxy
plotUnknowns <- function(resGC, unkn = "", DB = NULL, fileFrom = NULL) {
  ## Annotation table each value is a pcgrp associated to the unknown
  ## NOTE pcgrp index are different between xcmsSet and resGC due to filtering steps in metaMS
  ## R. Wehrens give me some clues on that and we found a correction

  # correction of annotation matrix due to pcgrp removal by quality check in runGCresult
  # matrix of correspondance between an@pspectra and filtered pspectra from runGC
  # Select only pspectra which correpond to them select in metaMS
  # col1 = filtered spectra from runGC and col2 = an@spectra
  allPCGRPs <- lapply(
    seq_along(resGC$xset),
    function(i) {
      an <- resGC$xset[[i]]
      huhn <- an@pspectra[which(sapply(an@pspectra, length) >=
        metaSetting(
          resGC$settings,
          "DBconstruction.minfeat"
        ))]
      matCORR <- cbind(seq_along(huhn), match(huhn, an@pspectra))
    }
  )
  # Build a new annotation list with sampnames and pseudospectra number from xset
  helpannotation <- list()
  for (j in seq_along(resGC$xset)) {
    helpannotation[[j]] <- resGC$annotation[[j]][1:2]
    pspvector <- vector()
    for (i in seq_len(nrow(helpannotation[[j]]))) {
      # Find corresponding pspec
      psplink <- allPCGRPs[[j]][match(helpannotation[[j]][i, 1], allPCGRPs[[j]]), 2]
      pspvector <- c(pspvector, psplink)
      # Change the annotation column into sampname column
      if (helpannotation[[j]][i, 2] < 0) {
        # It's an unknown
        new_name <- paste("Unknown", abs(as.integer(helpannotation[[j]][i, 2])))
        helpannotation[[j]][i, 2] <- new_name
      } else {
        # It has been found in local database
        for (k in seq_along(DB)) {
          if (helpannotation[[j]][i, 2] == k) {
            helpannotation[[j]][i, 2] <- DB[[k]]$Name
            break
          }
        }
      }
    }
    helpannotation[[j]] <- cbind(helpannotation[[j]], pspvector)
    names(helpannotation)[j] <- names(resGC$annotation[j])
  }
  peaktable <- resGC$PeakTable

  par(mar = c(5, 4, 4, 2) + 0.1)
  # For each unknown
  for (l in seq_along(unkn)) {
    # recordPlot
    perpage <- 3 # if change change layout also!
    dev.new(width = 21 / 2.54, height = 29.7 / 2.54, file = paste("Unknown_", unkn[l], ".pdf", sep = "")) # A4 pdf

    layout(matrix(c(1, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9), 6, 2, byrow = TRUE), widths = rep(c(1, 1), perpage), heights = rep(c(1, 5), perpage))

    oma.saved <- par("oma")
    par(oma = rep.int(0, 4))
    par(oma = oma.saved)
    o.par <- par(mar = rep.int(0, 4))
    on.exit(par(o.par))
    # For each sample
    for (c in seq_along(resGC$xset)) {
      # get sample name
      sampname <- basename(resGC$xset[[c]]@xcmsSet@filepaths)
      # remove .cdf, .mzXML filepattern
      filepattern <- c(
        "[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
        "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]"
      )
      filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
      sampname <- gsub(filepattern, "", sampname)
      title1 <- paste(peaktable[unkn[l], 1], "from", sampname, sep = " ")
      an <- resGC$xset[[c]]
      if (fileFrom == "zipfile") {
        an@xcmsSet@filepaths <- paste0("./", an@xcmsSet@phenoData[, "class"], "/", basename(an@xcmsSet@filepaths))
      }



      # Find the good annotation for this sample
      for (a in seq_len(length(helpannotation))) {
        if (gsub(filepattern, "", names(helpannotation)[a]) == paste0("./", sampname)) {
          # Find the unkn or the matched std in this sample
          findunkn <- FALSE
          for (r in seq_len(nrow(helpannotation[[a]]))) {
            if (helpannotation[[a]][r, "annotation"] == peaktable[unkn[l], 1]) {
              findunkn <- TRUE
              pcgrp <- helpannotation[[a]][r, "pspvector"]
              par(mar = c(0, 0, 0, 0) + 0.1)
              # Write title
              plot.new()
              box()
              text(0.5, 0.5, title1, cex = 2)
              par(mar = c(3, 2.5, 3, 1.5) + 0.1)
              # Window for EIC
              plotEICs(an, pspec = pcgrp, maxlabel = 2)
              # Window for pseudospectra
              plotPsSpectrum(an, pspec = pcgrp, maxlabel = 2)
            }
          }
          if (!findunkn) {
            par(mar = c(0, 0, 0, 0) + 0.1)
            # Write title
            plot.new()
            box()
            text(0.5, 0.5, title1, cex = 2)
            # Window for EIC
            plot.new()
            box()
            text(0.5, 0.5, "NOT FOUND", cex = 2)
            # Window for pseudospectra
            plot.new()
            box()
            text(0.5, 0.5, "NOT FOUND", cex = 2)
          }
          break
        }
      }
    }
    graphics.off()
  } # end  for unkn[l]
} # end function
