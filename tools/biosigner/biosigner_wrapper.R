#!/usr/bin/env Rscript


library(batch) ## parseCommandArgs

# Constants
argv <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=","",argv[grep("--file=",argv)])
prog.name <- basename(script.path)

# Print help
if (length(grep('-h', argv)) >0) {
	cat("Usage:", prog.name,
	    "data_matrix myDataMatrix.tsv",
	    "sample_metadata mySampleData.tsv",
	    "variable_metadata myVariableMetadata.tsv",
		"response_column ...",
  		"classification_method ...",
  		"num_bootstraps ...",
  		"selection_tier ...",
  		"pvalue_threshold ...",
  		"random_seed ...",
	    "variable_metadata_out myVariableMetadata_out.tsv",
	    "tier_figure tier_figure.pdf",
	    "boxplot_figure boxplot_figure.pdf",
	    "log_information log_information.txt",
		"\n")
	quit(status = 0)
}

# Parse all arguments
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

suppressMessages(library(biosigner))

if(packageVersion("biosigner") < "1.0.0")
    stop("Please use 'biosigner' versions of 1.0.0 and above")
if(packageVersion("ropls") < "1.4.0")
    stop("Please use 'ropls' versions of 1.4.0 and above")

## constants
##----------

modNamC <- "Biosigner" ## module name

topEnvC <- environment()
flgC <- "\n"

## functions
##----------

flgF <- function(tesC,
                 envC = topEnvC,
                 txtC = NA) { ## management of warning and error messages

    tesL <- eval(parse(text = tesC), envir = envC)

    if(!tesL) {

        sink(NULL)
        stpTxtC <- ifelse(is.na(txtC),
                          paste0(tesC, " is FALSE"),
                          txtC)

        stop(stpTxtC,
             call. = FALSE)

    }

} ## flgF


## log file
##---------

sink(argVc["log_information"])

cat("\nStart of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")


## arguments
##----------

xMN <- t(as.matrix(read.table(argVc["data_matrix"],
                              check.names = FALSE,
                              comment.char = '',
                              header = TRUE,
                              row.names = 1,
                              sep = "\t")))

samDF <- read.table(argVc["sample_metadata"],
                    check.names = FALSE,
                    comment.char = '',
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")
flgF("identical(rownames(xMN), rownames(samDF))", txtC = "Sample names (or number) in the data matrix (first row) and sample metadata (first column) are not identical; use the 'Check Format' module in the 'Quality Control' section")

varDF <- read.table(argVc["variable_metadata"],
                    check.names = FALSE,
                    comment.char = '',
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")
flgF("identical(colnames(xMN), rownames(varDF))", txtC = "Variable names (or number) in the data matrix (first column) and sample metadata (first column) are not identical; use the 'Check Format' module in the 'Quality Control' section")

flgF("argVc['response_column'] %in% colnames(samDF)",
     txtC = paste0("Class argument (", argVc['response_column'], ") must be either none or one of the column names (first row) of your sample metadata"))
respVc <- samDF[, argVc["response_column"]]
flgF("mode(respVc) == 'character'",
     txtC = paste0("'", argVc['response_column'], "' column of sampleMetadata does not contain only characters"))
respFc <- factor(respVc)
flgF("length(levels(respFc)) == 2",
     txtC = paste0("'", argVc['response_column'], "' column of sampleMetadata does not contain only 2 types of characters (e.g., 'case' and 'control')"))
tierMaxC <- ifelse("selection_tier" %in% names(argVc), argVc["selection_tier"], "S")
pvalN <- ifelse("pvalue_threshold" %in% names(argVc), as.numeric(argVc["pvalue_threshold"]), 0.05)


##------------------------------
## Computation and plot
##------------------------------


sink()

optWrnN <- options()$warn
options(warn = -1)

if("random_seed" %in% names(argVc) && argVc["random_seed"] != "0")
    set.seed(as.integer(argVc["random_seed"]))

bsnLs <- biosign(xMN,
                 respFc,
                 methodVc = ifelse("classification_method" %in% names(argVc), argVc["classification_method"], "all"),
                 bootI = ifelse("num_bootstraps" %in% names(argVc), as.numeric(argVc["num_bootstraps"]), 50),
                 pvalN = pvalN)

if("random_seed" %in% names(argVc) && argVc["random_seed"] != "0")
    set.seed(NULL)

tierMC <- bsnLs@tierMC

if(!is.null(tierMC)) {
    pdf(argVc["tier_figure"])
    plot(bsnLs, typeC = "tier")
    dev.off()
    pdf(argVc["boxplot_figure"])
    plot(bsnLs, typeC = "boxplot")
    dev.off()
} else {
    pdf(argVc["tier_figure"])
    plot(1, bty = "n", type = "n",
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),
         labels = "No significant variable to display")
    dev.off()
    pdf(argVc["boxplot_figure"])
    plot(1, bty = "n", type = "n",
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),
         labels = "No significant variable to display")
    dev.off()
}


options(warn = optWrnN)


##------------------------------
## Print
##------------------------------

sink(argVc["log_information"], append = TRUE)

tierFullVc <- c("S", LETTERS[1:5])
tierVc <- tierFullVc[1:which(tierFullVc == tierMaxC)]

if(sum(tierMC %in% tierVc)) {
    cat("\nSignificant features from '", paste(tierVc, collapse = "', '"), "' tiers:\n", sep = "")
    print(tierMC[apply(tierMC, 1, function(rowVc) sum(rowVc %in% tierVc) > 0), ,
                         drop = FALSE])
    cat("\nAccuracy:\n")
    print(round(getAccuracyMN(bsnLs), 3))
} else
    cat("\nNo significant variable found for any classifier\n")


##------------------------------
## Ending
##------------------------------

## Saving
##-------

if(!is.null(tierMC)) {
    tierDF <- data.frame(tier = sapply(rownames(varDF),
                             function(varC) {
                                 varTirVc <- tierMC[varC, ]
                                 varTirVc <- names(varTirVc)[varTirVc %in% tierVc]
                                 paste(varTirVc, collapse = "|")
                             }),
                         stringsAsFactors = FALSE)
    colnames(tierDF) <- paste(argVc["response_column"],
                              colnames(tierDF),
                              paste(tierVc, collapse = ""),
                              sep = "_")
    varDF <- cbind.data.frame(varDF, tierDF)
}

## variableMetadata

varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                          varDF)
write.table(varDF,
            file = argVc["variable_metadata_out"],
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

options(stringsAsFactors = strAsFacL)

rm(list = ls())
