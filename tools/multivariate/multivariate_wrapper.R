#!/usr/bin/env Rscript

library(batch) ## parseCommandArgs

# Constants
argv <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=","",argv[grep("--file=",argv)])
prog.name <- basename(script.path)

# Print help
if (length(grep('-h', argv)) >0) {
	cat("Usage:", prog.name,
	    "dataMatrix_in myDataMatrix.tsv",
	    "sampleMetadata_in mySampleData.tsv",
	    "variableMetadata_in myVariableMetadata.tsv",
		"respC ...",
		"predI ...",
		"orthoI ...",
		"testL ...",
		"typeC ...",
		"parAsColC ...",
		"parCexN ...",
		"parPc1I ...",
		"parPc2I ...",
		"parMahalC ...",
		"parLabVc ...",
		"algoC ...",
		"crossvalI ...",
		"log10L ...",
		"permI ...",
		"scaleC ...",
	    "sampleMetadata_out mySampleMetadata_out.tsv",
	    "variableMetadata_out myVariableMetadata_out.tsv",
	    "figure figure.pdf",
	    "information information.txt",
		"\n")
	quit(status = 0)
}

########
# MAIN #
########

argVc <- unlist(parseCommandArgs(evaluate=FALSE))

##------------------------------
## Initializing
##------------------------------

## options
##--------

strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## libraries
##----------

suppressMessages(library(ropls))

if(packageVersion("ropls") < "1.4.0")
    stop("Please use 'ropls' versions of 1.4.0 and above")

## constants
##----------

modNamC <- "Multivariate" ## module name

topEnvC <- environment()
flgC <- "\n"

## functions
##----------

flgF <- function(tesC,
                 envC = topEnvC,
                 txtC = NA) { ## management of warning and error messages

    tesL <- eval(parse(text = tesC), envir = envC)

    if(!tesL) {

        sink()
        stpTxtC <- ifelse(is.na(txtC),
                          paste0(tesC, " is FALSE"),
                          txtC)

        stop(stpTxtC,
             call. = FALSE)

    }

} ## flgF


## log file
##---------

sink(argVc["information"])

cat("\nStart of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")


## arguments
##----------

xMN <- t(as.matrix(read.table(argVc["dataMatrix_in"],
                              check.names = FALSE,
                              header = TRUE,
                              row.names = 1,
                              sep = "\t",
                              comment.char = "")))

samDF <- read.table(argVc["sampleMetadata_in"],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t",
                    comment.char = "")
flgF("identical(rownames(xMN), rownames(samDF))", txtC = "Sample names (or number) in the data matrix (first row) and sample metadata (first column) are not identical; use the 'Check Format' module in the 'Quality Control' section")

varDF <- read.table(argVc["variableMetadata_in"],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t",
                    comment.char = "")
flgF("identical(colnames(xMN), rownames(varDF))", txtC = "Variable names (or number) in the data matrix (first column) and sample metadata (first column) are not identical; use the 'Check Format' module in the 'Quality Control' section")

flgF("argVc['respC'] == 'none' || (argVc['respC'] %in% colnames(samDF))",
     txtC = paste0("Y Response argument (", argVc['respC'], ") must be either none or one of the column names (first row) of your sample metadata"))
if(argVc["respC"] != "none") {
    yMCN <- matrix(samDF[, argVc['respC']], ncol = 1, dimnames = list(rownames(xMN), argVc['respC']))
} else
    yMCN <- NULL

if(argVc["testL"] == "TRUE") {
    flgF("!is.null(yMCN)",
         txtC = "Predictions cannot be peformed with PCA models")
    flgF("'test.' %in% colnames(samDF)",
         txtC = "No 'test.' column found in the sample metadata")
    flgF("identical(sort(unique(samDF[, 'test.'])), c('no', 'yes'))",
         txtC = "'test.' column of sample metadata must contain both 'yes' (tested samples) and 'no' (samples to be used for model training) values, and nothing else")
    flgF("identical(sort(unique(samDF[, 'test.'])), c('no', 'yes'))",
         txtC = "'test.' column of sample metadata must contain both 'yes' (tested samples) and 'no' (samples to be used for model training) values, and nothing else")
    flgF("!any(is.na(yMCN[samDF[, 'test.'] == 'no', ]))",
         txtC = "samples for model training (i.e. 'no' value in the 'test.' column) should not contain NA in the response")
    tesVl <- samDF[, "test."] == "yes"
    xTesMN <- xMN[tesVl, , drop = FALSE]
    xMN <- xMN[!tesVl, , drop = FALSE]
    yMCN <- yMCN[!tesVl, , drop = FALSE]
} else
    tesVl <- NULL

if(!('parAsColC' %in% names(argVc)))
    argVc["parAsColC"] <- "none"
flgF("argVc['parAsColC'] == 'none' || argVc['parAsColC'] %in% colnames(samDF)", txtC = paste0("Sample color argument (", argVc['parAsColC'], ") must be either none or one of the column names (first row) of your sample metadata"))
if(argVc["parAsColC"] != "none") {
    parAsColFcVn <- samDF[, argVc['parAsColC']]
    if(!is.null(tesVl)) {
      parAsColFcVn <- parAsColFcVn[!tesVl]
    }
    if(is.character(parAsColFcVn))
        parAsColFcVn <- factor(parAsColFcVn)
} else
    parAsColFcVn <- NA

if(!('parMahalC' %in% names(argVc)) || argVc["parMahalC"] == "NA") {
    if(!is.null(yMCN) && ncol(yMCN) == 1 && mode(yMCN) == "character")
        argVc["parMahalC"] <- argVc["respC"]
    else
        argVc["parMahalC"] <- "none"
}
flgF("argVc['parMahalC'] == 'none' || (argVc['parMahalC'] %in% colnames(samDF))",
     txtC = paste0("Mahalanobis argument (", argVc['parMahalC'], ") must be either 'NA', 'none' or one of the column names (first row) of your sample metadata"))
if(argVc["parMahalC"] == "none") {
    parEllipsesL <- FALSE
} else {
    if(is.null(yMCN)) { ## PCA case
        flgF("mode(samDF[, argVc['parMahalC']]) == 'character'",
             txtC = paste0("Mahalanobis argument (", argVc['parMahalC'], ") must correspond to a column of characters in your sampleMetadata"))
        parAsColFcVn <- factor(samDF[, argVc["parMahalC"]])
        parEllipsesL <- TRUE
    } else { ## (O)PLS-DA case
        flgF("identical(as.character(argVc['respC']), as.character(argVc['parMahalC']))",
             txtC = paste0("The Mahalanobis argument (", argVc['parMahalC'], ") must be identical to the Y response argument (", argVc['respC'], ")"))
        parEllipsesL <- TRUE
    }
}

if(!('parLabVc' %in% names(argVc)))
    argVc["parLabVc"] <- "none"
flgF("argVc['parLabVc'] == 'none' || (argVc['parLabVc'] %in% colnames(samDF))",
     txtC = paste0("Sample labels argument (", argVc['parLabVc'], ") must be either none or one of the column names (first row) of your sample metadata"))
if('parLabVc' %in% names(argVc))
    if(argVc["parLabVc"] != "none") {
        flgF("mode(samDF[, argVc['parLabVc']]) == 'character'",
             txtC = paste0("The sample label argument (", argVc['parLabVc'], ") must correspond to a sample metadata column of characters (not numerics)"))
        parLabVc <- samDF[, argVc['parLabVc']]
        if(!is.null(tesVl)) {
          parLabVc <- parLabVc[!tesVl]
        }
    } else
        parLabVc <- NA

if('parPc1I' %in% names(argVc)) {
    parCompVi <-  as.numeric(c(argVc["parPc1I"], argVc["parPc2I"]))
} else
    parCompVi <- c(1, 2)


## checking
##---------


flgF("argVc['predI'] == 'NA' || argVc['orthoI'] == 'NA' || as.numeric(argVc['orthoI']) > 0 || parCompVi[2] <=  as.numeric(argVc['predI'])",
     txtC = paste0("The highest component to display (", parCompVi[2], ") must not exceed the number of predictive components of the model (", argVc['predI'], ")"))


if(argVc["orthoI"] == "NA" || argVc["orthoI"] != "0")
    if(argVc["predI"] == "NA" || argVc["predI"] != "0") {
        argVc["predI"] <- "1"
        cat("\nWarning: OPLS: number of predictive components ('predI' argument) set to 1\n", sep = "")
    }

if(argVc["predI"] != "NA")
    if(as.numeric(argVc["predI"]) > min(nrow(xMN), ncol(xMN))) {
        argVc["predI"] <- as.character(min(nrow(xMN), ncol(xMN)))
        cat("\nWarning: 'predI' set to the minimum of the dataMatrix dimensions: ", as.numeric(argVc["predI"]), "\n", sep = "")
    }

if("algoC" %in% names(argVc) && argVc["algoC"] == "svd" && length(which(is.na(c(xMN)))) > 0) {
    minN <- min(c(xMN[!is.na(xMN)])) / 2
    cat("\nWarning: Missing values set to ", round(minN, 1), " (half minimum value) for 'svd' algorithm to be used\n", sep = "")
}


##------------------------------
## Computation and plot
##------------------------------


sink()

optWrnN <- options()$warn
options(warn = -1)

ropLs <- opls(x = xMN,
              y = yMCN,
              predI = ifelse(argVc["predI"] == "NA", NA, as.numeric(argVc["predI"])),
              orthoI = ifelse(argVc["orthoI"] == "NA", NA, as.numeric(argVc["orthoI"])),
              algoC = ifelse('algoC' %in% names(argVc), argVc["algoC"], "default"),
              crossvalI = ifelse('crossvalI' %in% names(argVc), as.numeric(argVc["crossvalI"]), 7),
              log10L = ifelse('log10L' %in% names(argVc), as.logical(argVc["log10L"]), FALSE),
              permI = ifelse('permI' %in% names(argVc), as.numeric(argVc["permI"]), 20),
              scaleC = ifelse('scaleC' %in% names(argVc), argVc["scaleC"], "standard"),
              subset = NULL,
              printL = FALSE,
              plotL = FALSE,
              .sinkC = argVc['information'])

modC <- ropLs@typeC
sumDF <- getSummaryDF(ropLs)
desMC <- ropLs@descriptionMC
scoreMN <- getScoreMN(ropLs)
loadingMN <- getLoadingMN(ropLs)

vipVn <- coeMN <- orthoScoreMN <- orthoLoadingMN <- orthoVipVn <- NULL

if(grepl("PLS", modC)) {

    vipVn <- getVipVn(ropLs)
    coeMN <- coef(ropLs)

    if(grepl("OPLS", modC)) {
        orthoScoreMN <- getScoreMN(ropLs, orthoL = TRUE)
        orthoLoadingMN <- getLoadingMN(ropLs, orthoL = TRUE)
        orthoVipVn <- getVipVn(ropLs, orthoL = TRUE)
    }

}

ploC <- ifelse('typeC' %in% names(argVc), argVc["typeC"], "summary")

if(sumDF[, "pre"] + sumDF[, "ort"] < 2) {
    if(!(ploC %in% c("permutation", "overview"))) {
        ploC <- "summary"
        plotWarnL <- TRUE
    }
} else
    plotWarnL <- FALSE

plot(ropLs,
     typeVc = ploC,
     parAsColFcVn = parAsColFcVn,
     parCexN = ifelse('parCexN' %in% names(argVc), as.numeric(argVc["parCexN"]), 0.8),
     parCompVi = parCompVi,
     parEllipsesL = parEllipsesL,
     parLabVc = parLabVc,
     file.pdfC = argVc['figure'],
     .sinkC = argVc['information'])

options(warn = optWrnN)


##------------------------------
## Print
##------------------------------


sink(argVc["information"], append = TRUE)

if(plotWarnL)
    cat("\nWarning: For single component models, only 'overview' (and 'permutation' in case of single response (O)PLS(-DA)) plot(s) are available\n", sep = "")


cat("\n", modC, "\n", sep = "")

cat("\n", desMC["samples", ],
    " samples x ",
    desMC["X_variables", ],
    " variables",
    ifelse(modC != "PCA",
           " and 1 response",
           ""),
    "\n", sep = "")

cat("\n", ropLs@suppLs[["scaleC"]], " scaling of dataMatrix",
            ifelse(modC == "PCA",
                   "",
                   paste0(" and ",
                          ifelse(mode(ropLs@suppLs[["yMCN"]]) == "character" && ropLs@suppLs[["scaleC"]] != "standard",
                                 "standard scaling of ",
                                 ""),
                          "response\n")), sep = "")

if(substr(desMC["missing_values", ], 1, 1) != "0")
    cat("\n", desMC["missing_values", ], " NAs\n", sep = "")

if(substr(desMC["near_zero_excluded_X_variables", ], 1, 1) != "0")
    cat("\n", desMC["near_zero_excluded_X_variables", ],
        " excluded variables during model building (because of near zero variance)\n", sep = "")

cat("\n")

optDigN <- options()[["digits"]]
options(digits = 3)
print(ropLs@modelDF)
options(digits = optDigN)


##------------------------------
## Ending
##------------------------------


## Saving
##-------


rspModC <- gsub("-", "", modC)
if(rspModC != "PCA")
    rspModC <- paste0(make.names(argVc['respC']), "_", rspModC)

if(sumDF[, "pre"] + sumDF[, "ort"] < 2) {

    tCompMN <- scoreMN
    pCompMN <- loadingMN

} else {

    if(sumDF[, "ort"] > 0) {
        if(parCompVi[2] > sumDF[, "ort"] + 1)
            stop("Selected orthogonal component for plotting (ordinate) exceeds the total number of orthogonal components of the model", call. = FALSE)
        tCompMN <- cbind(scoreMN[, 1], orthoScoreMN[, parCompVi[2] - 1])
        pCompMN <- cbind(loadingMN[, 1], orthoLoadingMN[, parCompVi[2] - 1])
        colnames(pCompMN) <- colnames(tCompMN) <- c("h1", paste("o", parCompVi[2] - 1, sep = ""))
    } else {
        if(max(parCompVi) > sumDF[, "pre"])
            stop("Selected component for plotting as ordinate exceeds the total number of predictive components of the model", call. = FALSE)
        tCompMN <- scoreMN[, parCompVi, drop = FALSE]
        pCompMN <- loadingMN[, parCompVi, drop = FALSE]
    }

}

## x-scores and prediction

colnames(tCompMN) <- paste0(rspModC, "_XSCOR-", colnames(tCompMN))
tCompDF <- as.data.frame(tCompMN)[rownames(samDF), , drop = FALSE]

if(modC != "PCA") {

    if(!is.null(tesVl)) {
        tCompFulMN <- matrix(NA,
                             nrow = nrow(samDF),
                             ncol = ncol(tCompMN),
                             dimnames = list(rownames(samDF), colnames(tCompMN)))
        mode(tCompFulMN) <- "numeric"
        tCompFulMN[rownames(tCompMN), ] <- tCompMN
        tCompMN <- tCompFulMN

        fitMCN <- fitted(ropLs)
        fitFulMCN <- matrix(NA,
                            nrow = nrow(samDF),
                            ncol = 1,
                            dimnames = list(rownames(samDF), NULL))
        mode(fitFulMCN) <- mode(fitMCN)
        fitFulMCN[rownames(fitMCN), ] <- fitMCN
        yPreMCN <- predict(ropLs, newdata = as.data.frame(xTesMN))
        fitFulMCN[rownames(yPreMCN), ] <- yPreMCN
        fitMCN <- fitFulMCN

    } else
        fitMCN <- fitted(ropLs)

    colnames(fitMCN) <- paste0(rspModC,
                               "_predictions")
    fitDF <- as.data.frame(fitMCN)[rownames(samDF), , drop = FALSE]

    tCompDF <- cbind.data.frame(tCompDF, fitDF)
}

samDF <- cbind.data.frame(samDF, tCompDF)

## x-loadings and VIP

colnames(pCompMN) <- paste0(rspModC, "_XLOAD-", colnames(pCompMN))
if(!is.null(vipVn)) {
    pCompMN <- cbind(pCompMN, vipVn)
    colnames(pCompMN)[ncol(pCompMN)] <- paste0(rspModC,
                                               "_VIP",
                                               ifelse(!is.null(orthoVipVn),
                                                      "_pred",
                                                      ""))
    if(!is.null(orthoVipVn)) {
        pCompMN <- cbind(pCompMN, orthoVipVn)
        colnames(pCompMN)[ncol(pCompMN)] <- paste0(rspModC,
                                                   "_VIP_ortho")
    }
}
if(!is.null(coeMN)) {
    pCompMN <- cbind(pCompMN, coeMN)
    if(ncol(coeMN) == 1)
        colnames(pCompMN)[ncol(pCompMN)] <- paste0(rspModC, "_COEFF")
    else
        colnames(pCompMN)[(ncol(pCompMN) - ncol(coeMN) + 1):ncol(pCompMN)] <- paste0(rspModC, "_", colnames(coeMN), "-COEFF")
}
pCompDF <- as.data.frame(pCompMN)[rownames(varDF), , drop = FALSE]
varDF <- cbind.data.frame(varDF, pCompDF)

## sampleMetadata

samDF <- cbind.data.frame(sampleMetadata = rownames(samDF),
                          samDF)
write.table(samDF,
            file = argVc["sampleMetadata_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

## variableMetadata

varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                          varDF)
write.table(varDF,
            file = argVc["variableMetadata_out"],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

# Output ropLs
if (!is.null(argVc['ropls_out']) && !is.na(argVc['ropls_out']))
    save(ropLs, file = argVc['ropls_out'])

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

options(stringsAsFactors = strAsFacL)

rm(list = ls())
