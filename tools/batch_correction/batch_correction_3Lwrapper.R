#!/usr/bin/env Rscript

################################################################################################
# batch_correction_wrapper                                                                     #
#                                                                                              #
# Authors: Marion LANDI / Jean-Francois MARTIN / Melanie Petera                                #
# User: Galaxy                                                                                 #
# Original data: --                                                                            #
# Starting date: 22-07-2014                                                                    #
# Version 1: 22-07-2014                                                                        #
# Version 2: 08-12-2014                                                                        #
# Version 2.1: 09-01-2015 modification in Error message of sample matching                     #
# Version 2.2: 16-03-2015 inclusion of miniTools' functions for special characters             #
# Version 2.90: 18-08-2015 new parameter valnull                                               #
# Version 2.91: 25-08-2016 error message improvment                                            #
# Version 3: 02-10-2020                                                                        #
#            - split of tool-linked code and script-linked one                                 #
#            - addition of args print and sessionInfo()                                        #
#            - adjustment of sample tags' parameters to 3L methods                             #
#            - addition of the min.norm argument in meth3L() call                              #
#                                                                                              #
# Input files: dataMatrix.txt, sampleMetadata.txt, variableMetadata.txt (BC only)              #
# Output files: graph.pdf, corrected table (BC only), diagnostic table (DBC only),             #
#               variableMetadata (BC only)                                                     #
#                                                                                              #
################################################################################################


library(W4MRUtils)

## ------------------------------
## test help option
## ------------------------------

# Prog. constants
argv.help <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=", "", argv.help[grep("--file=", argv.help)])
prog.name <- basename(script.path)

# Test Help
if (length(grep("-h", argv.help)) > 0) {
    cat(
        "Usage: Rscript ",
        prog.name,
        "{args} \n",
        "parameters: \n",
        "\tanalyse {val}: must be set to \"batch_correction\"",
        "\tdataMatrix {file}: set the input data matrix file (mandatory) \n",
        "\tsampleMetadata {file}: set the input sample metadata file (mandatory) \n",
        "\tvariableMetadata {file}: set the input variable metadata file (mandatory) \n",
        "\tmethod {opt}: set the method; can set to \"linear\", \"lowess\" or \"loess\" (mandatory) \n",
        "\tspan {condition}: set the span condition; set to \"none\" if method is set to \"linear\" (mandatory) \n",
        "\tref_factor {value}: set the ref_factor value; (if span value is set to NULL, optional) \n",
        "\tdetail {value}: set the detail value; (if span value is set to NULL, optional) \n",
        "\tdataMatrix_out {file}: set the output data matrix file (mandatory) \n",
        "\tvariableMetadata_out {file}: set the output variable metadata file (mandatory) \n",
        "\tgraph_output {file}: set the output graph file (mandatory) \n",
        "\trdata_output {file}: set the output Rdata file (mandatory) \n",
        "\tbatch_col_name {val}: the column name for batch. Default value is \"batch\".\n",
        "\tinjection_order_col_name {val}: the column name for the injection order. Default value is \"injectionOrder\".\n",
        "\tsample_type_col_name {val}: the column name for the sample types. Default value is \"sampleType\".\n",
        "\tsample_type_tags {val}: the tags used inside the sample type column, defined as key/value pairs separated by commas (example: blank=blank,pool=pool,sample=sample).\n",
        "\n"
    )
    quit(status = 0)
}

## ------------------------------
## init. params
## ------------------------------

args <- parse_args() # interpretation of arguments given in command line as an R list of objects


cat(
    "\nJob starting time:\n", format(Sys.time(), "%a %d %b %Y %X"),
    "\n\n--------------------------------------------------------------------",
    "\nParameters used:\n\n"
)
print(args)
cat("--------------------------------------------------------------------\n\n")


# Set default col names
if (!"batch_col_name" %in% names(args)) {
    args[["batch_col_name"]] <- "batch"
}
if (!"injection_order_col_name" %in% names(args)) {
    args[["injection_order_col_name"]] <- "injectionOrder"
}
if (!"sample_type_col_name" %in% names(args)) {
    args[["sample_type_col_name"]] <- "sampleType"
}
if (!"sample_type_tags" %in% names(args)) {
    args[["sample_type_tags"]] <- "blank=blank,pool=pool,sample=sample"
}

# Parse sample type tags
sample.type.tags <- list()
for (kv in strsplit(strsplit(args$sample_type_tags, ",")[[1]], "=")) {
    sample.type.tags[[kv[[1]]]] <- kv[-1]
}
if (!all(c("pool", "blank", "sample") %in% names(sample.type.tags))) {
    stop("All tags pool, blank and sample must be defined in option sampleTypeTags.")
}
args$sample_type_tags <- sample.type.tags

## ------------------------------
## init. functions
## ------------------------------

source_local <- function(...) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    for (i in 1:length(list(...))) {
        source(paste(base_dir, list(...)[[i]], sep = "/"))
    }
}
# Import the different functions
source_local("batch_correction_3Lfct.R", "batch_correction_3Llauncher.R")

# Specificities of BC and DBC
if (args$analyse == "batch_correction") {
    args$out_graph_pdf <- NULL
    args$out_preNormSummary <- NULL
} else {
    args$variableMetadata <- NULL
    args$rdata_output <- NULL
    args$dataMatrix_out <- NULL
    args$variableMetadata_out <- NULL
    args$graph_output <- NULL
    args$method <- NULL
    args$detail <- NULL
    args$valnull <- NULL
}

# Launch tool
meth3L(
    idsample = args$sampleMetadata, iddata = args$dataMatrix, sample_type_col_name = args$sample_type_col_name, injection_order_col_name = args$injection_order_col_name,
    batch_col_name = args$batch_col_name, sample_type_tags = args$sample_type_tags, factbio = args$ref_factor, analyse = args$analyse, metaion = args$variableMetadata,
    detail = args$detail, method = args$method, outlog = args$graph_output, span = args$span, valnull = args$valnull, rdata_output = args$rdata_output,
    dataMatrix_out = args$dataMatrix_out, variableMetadata_out = args$variableMetadata_out, out_graph_pdf = args$out_graph_pdf, out_preNormSummary = args$out_preNormSummary,
    min.norm = 1
)


cat(
    "\n\n--------------------------------------------------------------------",
    "\nInformation about R (version, Operating System, attached or loaded packages):\n\n"
)
sessionInfo()
cat(
    "--------------------------------------------------------------------\n",
    "\nJob ending time:\n", format(Sys.time(), "%a %d %b %Y %X")
)

rm(args)
