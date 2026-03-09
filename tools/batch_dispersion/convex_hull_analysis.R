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
  make_option(c("-m", "--variables"),
            type = "character",
            default = "",
            help = "variableMetadata containing the data",
            metavar = "NAME"),
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
  make_option(c("--output_plot"),
              type = "character",
              default = "plot_convex_hull.pdf",
              help = "Name of the output pdf file plot",
              metavar = "OUTPUT"),
  make_option(c("--output_vm"),
              type = "character",
              default = "variableMetadata.tsv",
              help = "Name of the output tabular file variableMetadata",
              metavar = "OUTPUT")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#### ---- Read data function ----
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

#### ---- Load data and extract script parameters
dataMatrix <- read_data_file(opt$dataMatrix, "dataMatrix")
sampleMetadata <- read_data_file(opt$sampleMetadata, "sampleMetadata")
variableMetadata <- read_data_file(opt$variableMetadata, "variableMetadata")
variables <- strsplit(opt$variables, ",")[[1]]
batch_col <- names(sampleMetadata)[as.integer(opt$batch)]
order_col <- names(sampleMetadata)[as.integer(opt$sample_order_col)]
mode <- "batchwise"
if (opt$global) {
  mode <- "global"
}

#### ---- Create data ----
variableMetadata <- data.frame(variableMetadata, row.names = 1)
dataMatrix_t <- as.data.frame(t(dataMatrix[-1]))
colnames(dataMatrix_t) <- dataMatrix[[1]]         
dataMatrix_t$sampleMetadata <- rownames(dataMatrix_t)
dataMatrix_t <- dataMatrix_t[, c("sampleMetadata", setdiff(names(dataMatrix_t), "sampleMetadata"))]
pool_s <- merge(sampleMetadata, dataMatrix_t, by = "sampleMetadata")
pool_s <- pool_s[order(pool_s[[order_col]]), ]
if (length(variables) == 0) {
  variable_columns <- rownames(variableMetadata)
} else {
  variable_columns <- variables
}
#### ---- Call function for convex analysis ----
result <- convex_analysis_of_variables(
  pool_s,
  variable_columns=variable_columns,
  batch_col=batch_col,
  sample_order_col=order_col,
  impute_if_needed="median",
  mode=mode
)
#### ---- Adding dispersion indicators to variable metadata ----
variableMetadata[c(
  "IntraBatchConvexDispersion",
  "InterBatchConvexDispersion",
  "IntraInterBatchDispersonRatio"
)] <- result$indicators[c("IntraB", "InterB", "Ratio")]
variableMetadata <- data.frame(variableMetadata=rownames(variableMetadata), variableMetadata)
#### ---- Creating new variable metadata file ----
write.table(
  variableMetadata[variable_columns, ],
  file = opt$output_vm,
  sep="\t",
  row.names = FALSE
)
cat("VM saved as", opt$output_vm, "\n")
#### ---- Call to ploting function ----
tryCatch(
  {
    plot_all_convex_hulls(
    target_file_path = opt$output_plot,
    convex_analysis_res = result,
    show_points = opt$points,
    mode = mode
    )
    cat("Plot saved as", opt$output_plot, "\n")
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


