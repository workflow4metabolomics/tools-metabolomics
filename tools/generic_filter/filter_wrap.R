#############################################################################
# Wrapper for the generic_filter() function from the W4MRUtils R package
# XML associated to the tool: generic_filter.xml
#
# Input data: Data Matrix, Variable Metadata, Sample Metadata
# Output data: Data Matrix, Variable Metadata, Sample Metadata
#############################################################################

suppressPackageStartupMessages(library(W4MRUtils))

# Constants
argv <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=", "", argv[grep("--file=", argv)])
prog.name <- basename(script.path)

# Help
if (length(grep("-h", argv)) > 0) {
    cat(
        "Usage:", prog.name,
        "dataMatrix_in myDataMatrix.tsv",
        "sampleMetadata_in mySampleMetadata.tsv",
        "variableMetadata_in myVariableMetadata.tsv",
        "...",
        "\n"
    )
    quit(status = 0)
}

# Parameter formating -----------------------------------------------------

args <- parse_args() # interpretation of arguments given in command line as an R list of objects

if (length(args) < 8) {
    stop("NOT enough argument!!!")
}

cat(
    "\nJob starting time:\n", format(Sys.time(), "%a %d %b %Y %X"),
    "\n\n--------------------------------------------------------------------",
    "\nParameters used in 'Generic Filter':\n\n"
)
print(args)
cat("--------------------------------------------------------------------\n\n")

list_num <- NULL
if (!is.null(args$parm_col)) {
    for (i in which(names(args) == "num_file")) {
        if (args[[i + 2]] %in% c("lower", "upper")) {
            list_num <- c(list_num, list(c(args[[i]], args[[i + 1]], args[[i + 2]], args[[i + 3]])))
        }
        if (args[[i + 2]] %in% c("between", "extremity")) {
            list_num <- c(list_num, list(c(args[[i]], args[[i + 1]], args[[i + 2]], args[[i + 3]], args[[i + 4]])))
        }
    }
}

list_fact <- NULL
if (!is.null(args$factor_col)) {
    for (i in which(names(args) == "qual_file")) {
        list_fact <- c(list_fact, list(c(args[[i + 1]], args[[i + 2]], args[[i]])))
    }
}

# Begining of processing --------------------------------------------------

data3tables <- W4MRUtils::import3(args$dataMatrix_in, args$sampleMetadata_in, args$variableMetadata_in)

filteredset <- generic_filter(
    data3tables$dataMatrix,
    data3tables$sampleMetadata,
    data3tables$variableMetadata,
    args$Numeric,
    list_num,
    args$Factors,
    list_fact
)

write.table(filteredset$dataMatrix, args$dataMatrix_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filteredset$sampleMetadata, args$sampleMetadata_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filteredset$variableMetadata, args$variableMetadata_out, sep = "\t", quote = FALSE, row.names = FALSE)

# End of processing -------------------------------------------------------

cat(
    "\n--------------------------------------------------------------------",
    "\nInformation about R (version, Operating System, attached or loaded packages):\n\n"
)
sessionInfo()
cat(
    "--------------------------------------------------------------------\n",
    "\nJob ending time:\n", format(Sys.time(), "%a %d %b %Y %X")
)

# delete the parameters to avoid the passage to the next tool in .RData image
rm(args)
