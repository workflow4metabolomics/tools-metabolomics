#############################################################################
# WRAPPER for INTENSITY CHECK                                               #
#                                                                           #
# Core function used: intens_check() from the W4MRUtils package             #
# XML associated to the tool: xml_intensity_check.xml                       #
#                                                                           #
# Input: Data Matrix, VariableMetadata, SampleMetadata                      #
# Output: VariableMetadata, Graphics                                        #
#                                                                           #
#############################################################################

suppressPackageStartupMessages(library(W4MRUtils))

args <- parse_args() # interpretation of arguments given in command line as an R list of objects

if (length(args) < 7) {
    stop("NOT enough argument !!!")
}

cat(
    "\nJob starting time:\n",
    format(Sys.time(), "%a %d %b %Y %X"),
    "\n\n--------------------------------------------------------------------",
    "\nIntensity Check parameters:\n\n"
)
print(args)
cat("--------------------------------------------------------------------\n\n")

class_col <- NULL
test_fold <- NULL
class1 <- NULL
fold_frac <- NULL
logarithm <- NULL
if (args$method == "each_class") {
    class_col <- args$class_col
    test_fold <- args$test_fold
    if (args$test_fold == "Yes") {
        logarithm <- args$logarithm
    }
}
if (args$method == "one_class") {
    class_col <- args$class_col
    class1 <- args$class1
    test_fold <- args$test_fold
    if (args$test_fold == "Yes") {
        fold_frac <- args$fold_frac
        logarithm <- args$logarithm
    }
}

err_no_option <- NULL

if (
    ((args$method == "no_class") && (args$chosen_stat == "None")) ||
        ((args$method != "no_class") &&
            (args$chosen_stat == "None") &&
            (test_fold == "No"))
) {
    err_no_option <- "You did not select any computational option. Program can not be executed."
    stop("\n- - - - - - - - -\n", err_no_option, "\n- - - - - - - - -\n")
}

# Begining of processing ----------------------------------------------------------------------------------------------

data3tables <- W4MRUtils::import3(args$dataMatrix_in, args$sampleMetadata_in, args$variableMetadata_in)
dm <- data3tables$dataMatrix
sm <- data3tables$sampleMetadata
vm <- data3tables$variableMetadata

if (is.null(err_no_option)) {
    VM.output <- intens_check(
        dm,
        sm,
        vm,
        args$method,
        args$chosen_stat,
        class_col,
        test_fold,
        class1,
        fold_frac,
        logarithm,
        args$graphs_out
    )
    write.table(
        VM.output,
        args$variableMetadata_out,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

# End of processing ---------------------------------------------------------------------------------------------------

cat(
    "\n--------------------------------------------------------------------",
    "\nInformation about R (version, Operating System, attached or loaded packages):\n\n"
)
sessionInfo()
cat(
    "--------------------------------------------------------------------\n",
    "\nJob ending time:\n",
    format(Sys.time(), "%a %d %b %Y %X")
)


# delete the parameters to avoid the passage to the next tool in .RData image
rm(args)
