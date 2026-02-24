#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(optparse)
library(dispersionIndicators)

#### ---- Define command-line options ----
option_list <- list(
  make_option(c("-q", "--dataMatrix"),
              type = "character",
              help = "dataMatrix containing the data",
              metavar = "FILE"),
  make_option(c("-s", "--sampleMetadata"),
              type = "character",
              help = "sampleMetadata containing the data",
              metavar = "FILE"),
  make_option(c("-v", "--variableMetadata"),
              type = "character",
              help = "variableMetadata containing the data",
              metavar = "FILE"),
  make_option(c("-m", "--metabolites"),
            type = "character",
            default = "",
            help = "variableMetadata containing the data",
            metavar = "FILE"),
  make_option(c("-g", "--global"),
              type = "logical",
              default = FALSE,
              help = "Injection Order Global used",
              metavar = "BOOL"),
  make_option(c("-x", "--batch_col"),
              type = "character",
              default = "batch",
              help = "Name of the column with batch information",
              metavar = "BOOL"),
  make_option(c("-t", "--sample_order_col"),
              type = "character",
              default = "injectionOrder",
              help = "Name of the columns with sample order information (e.g. injectionOrder)",
              metavar = "NAME"),
  make_option(c("-p", "--points"),
              type = "logical",
              default = TRUE,
              help = "Display points (TRUE/FALSE)",
              metavar = "BOOL"),
  make_option(c("-o", "--output"),
              type = "character",
              default = "plot_convex_hull.pdf",
              help = "Name of the output file",
              metavar = "OUTPUT")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#### ---- Read Data ---- 
read_data_file <- function(file, description) {
  if (!file.exists(file)) {
    stop(paste(description, "file does not exist:", file))
  }
  df <- tryCatch({
    read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  }, error = function(e) {
    stop(paste("Error reading", description, "file:", conditionMessage(e)))
  })
  return(df)
}

dataMatrix <- read_data_file(opt$dataMatrix, "dataMatrix")
sampleMetadata <- read_data_file(opt$sampleMetadata, "sampleMetadata")
variableMetadata <- read_data_file(opt$variableMetadata, "variableMetadata")

#### ---- Verification for tests ----
# Check the colname
if (as.integer(opt$sample_order_col) > length(names(sampleMetadata))) {
  stop("Error : Specified column number for sample order in the sampleMetadata File is higher that the number of
  columns ")
}
if (as.integer(opt$batch) > length(names(sampleMetadata))) {
  stop("Error : Specified column number for the batch in the sampleMetadata File is higher that the number of columns ")
}
batch_col <- names(sampleMetadata)[as.integer(opt$batch)]
order_col <- names(sampleMetadata)[as.integer(opt$sample_order_col)]
#### ---- Create data ----
variableMetadata <- data.frame(variableMetadata, row.names = 1)
dataMatrix_t <- as.data.frame(t(dataMatrix[-1]))
colnames(dataMatrix_t) <- dataMatrix[[1]]         
dataMatrix_t$sampleMetadata <- rownames(dataMatrix_t)
dataMatrix_t <- dataMatrix_t[, c("sampleMetadata", setdiff(names(dataMatrix_t), "sampleMetadata"))]
pool_s <- merge(sampleMetadata, dataMatrix_t, by = "sampleMetadata")
pool_s <- pool_s[order(pool_s[[order_col]]), ]

# # Use global injection order or not
mode <- "batchwise"
if (opt$global) {
  mode <- "global"
}
# if (!opt$global) {
#   pool_s[[opt$sample_order_col]] <- ave(
#     pool_s[[opt$sample_order_col]],
#     pool_s[[opt$batch_col]],
#     FUN = function(x) seq_along(x)
#   )
# }

#### ---- Call plotting function ----
result <- convex_analysis_of_variables(
  pool_s,
  variable_columns=rownames(variableMetadata),
  batch_col=batch_col,
  sample_order_col=order_col,
  impute_if_needed="median",
  mode=mode
)
tryCatch(
  {
    plot_all_convex_hulls(
    target_file_path = opt$output,
    convex_analysis_res = result,
    show_points = opt$points,
    mode = mode
    )
    cat("Plot saved as", opt$output, "\n")
  },
  warning = function (war) {
    print("Caught warning:")
    warning(war$message)
  },
  error = function (err) {
    print("Caught exception:")
    stop(err$message)
  }
)


