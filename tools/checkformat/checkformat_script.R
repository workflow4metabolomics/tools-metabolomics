## Etienne Thevenot
## CEA, MetaboHUB Paris
## etienne.thevenot@cea.fr



## Reads the dataMatrix, sampleMetadata, and variableMetadata .tsv files
## and checks the formats
readAndCheckF <- function(datFilC = "dataMatrix.tsv",
                          samFilC = "sampleMetadata.tsv",
                          varFilC = "variableMetadata.tsv",
                          makNamL) {
  ## options

  optStrAsFacL <- options()[["stringsAsFactors"]]
  options(stringsAsFactors = FALSE)


  ## checking that the tables have no duplicated row or column names

  for (tabC in c("dat", "sam", "var")) {
    tabNamC <- switch(tabC,
      dat = "dataMatrix",
      sam = "sampleMetadata",
      var = "variableMetadata"
    )

    rowVc <- read.table(eval(parse(text = paste0(tabC, "FilC"))),
      check.names = FALSE,
      header = TRUE,
      sep = "\t"
    )[, 1]

    colVc <- unlist(read.table(eval(parse(text = paste0(tabC, "FilC"))),
      check.names = FALSE,
      nrow = 1,
      sep = "\t"
    ))[-1]

    if (any(duplicated(rowVc))) {
      stop("The following row name(s) is/are duplicated in the ",
        tabNamC,
        " table: '",
        paste(rowVc[duplicated(rowVc)], collapse = "', '"), "'",
        call. = FALSE
      )
    }

    if (any(duplicated(colVc))) {
      stop("The following column name(s) is/are duplicated in the ",
        tabNamC,
        " table: '",
        paste(colVc[duplicated(colVc)], collapse = "', '"), "'",
        call. = FALSE
      )
    }
  }


  ## reading tables

  datMN <- t(as.matrix(read.table(datFilC,
    check.names = FALSE,
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )))

  samDF <- read.table(samFilC,
    check.names = FALSE,
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )

  varDF <- read.table(varFilC,
    check.names = FALSE,
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )


  ## checking that dataMatrix is numeric and that the sample and variable numbers are coherent

  if (mode(datMN) != "numeric") {
    stop("The dataMatrix is not of the 'numeric' type",
      call. = FALSE
    )
  }

  if (nrow(datMN) != nrow(samDF)) {
    if (nrow(datMN) > nrow(samDF)) {
      print(setdiff(rownames(datMN), rownames(samDF)))
      stop("The sample names above from dataMatrix were not found in sampleMetadata",
        call. = FALSE
      )
    } else {
      print(setdiff(rownames(samDF), rownames(datMN)))
      stop("The sample names above from sampleMetadata were not found in dataMatrix",
        call. = FALSE
      )
    }
  }

  if (ncol(datMN) != nrow(varDF)) {
    if (ncol(datMN) > nrow(varDF)) {
      print(setdiff(colnames(datMN), rownames(varDF)))
      stop("The variable names above from dataMatrix were not found in variableMetadata",
        call. = FALSE
      )
    } else {
      print(setdiff(rownames(varDF), colnames(datMN)))
      stop("The variable names above from variableMetadata were not found in dataMatrix",
        call. = FALSE
      )
    }
  }


  ## making sample and variable names (optional)

  newL <- FALSE

  if (makNamL) {
    cat("\n\nMessage: Converting sample and variable names to the standard R format\n")

    rownames(datMN) <- make.names(rownames(datMN), unique = TRUE)
    colnames(datMN) <- make.names(colnames(datMN), unique = TRUE)
    rownames(samDF) <- make.names(rownames(samDF), unique = TRUE)
    rownames(varDF) <- make.names(rownames(varDF), unique = TRUE)

    newL <- TRUE
  }


  ## checking sample and variable names

  chkL <- TRUE

  if (!identical(rownames(datMN), rownames(samDF))) {
    if (identical(sort(rownames(datMN)), sort(rownames(samDF)))) {
      cat("\n\nMessage: Re-ordering dataMatrix sample names to match sampleMetadata\n")
      datMN <- datMN[rownames(samDF), , drop = FALSE]

      stopifnot(identical(sort(rownames(datMN)), sort(rownames(samDF))))

      newL <- TRUE
    } else {
      cat("\n\nStop: The sample names of dataMatrix and sampleMetadata do not match:\n")
      print(cbind.data.frame(
        indice = 1:nrow(datMN),
        dataMatrix = rownames(datMN),
        sampleMetadata = rownames(samDF)
      )[rownames(datMN) != rownames(samDF), , drop = FALSE])
      chkL <- FALSE
    }
  }

  if (!identical(colnames(datMN), rownames(varDF))) {
    if (identical(sort(colnames(datMN)), sort(rownames(varDF)))) {
      cat("\n\nMessage: Re-ordering dataMatrix variable names to match variableMetadata\n")
      datMN <- datMN[, rownames(varDF), drop = FALSE]

      stopifnot(identical(sort(colnames(datMN)), sort(rownames(varDF))))

      newL <- TRUE
    } else {
      cat("\n\nStop: The variable names of dataMatrix and variableMetadata do not match:\n")
      print(cbind.data.frame(
        indice = 1:ncol(datMN),
        dataMatrix = colnames(datMN),
        variableMetadata = rownames(varDF)
      )[colnames(datMN) != rownames(varDF), , drop = FALSE])
      chkL <- FALSE
    }
  }


  options(stringsAsFactors = optStrAsFacL)

  resLs <- list(
    chkL = chkL,
    newL = newL,
    datMN = datMN,
    samDF = samDF,
    varDF = varDF
  )

  return(resLs)
} ## end of checkAndReadF
