#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(optparse)
library(tools)
library(dplyr)

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
  make_option(c("-g", "--global"),
              type = "logical",
              default = FALSE,
              help = "Injection Order Global used",
              metavar = "BOOL"),
  make_option(c("-x", "--multiple"),
              type = "logical",
              default = FALSE,
              help = "multiple plot depend of the number of Intensity",
              metavar = "BOOL"),
  make_option(c("-m", "--metabolite"),
              type = "character",
              default = "Metabolite X",
              help = "Name of the metabolite for the plot title (e.g. Intensity1)",
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
# Check the colnames
required_cols_data <- character(ncol(dataMatrix) - 1)
for (i in 1:(ncol(dataMatrix) - 1)) {
  required_cols_data[i] <- paste0("S", i)
}
missing_cols <- setdiff(required_cols_data, names(dataMatrix))
if (length(missing_cols) > 0) {
  stop(paste("Error : Missing columns in the dataMatrix File :", paste(missing_cols, collapse = ", ")))
}

required_cols_sample <- c("InjectionOrder", "Batch")
missing_cols <- setdiff(required_cols_sample, names(sampleMetadata))
if (length(missing_cols) > 0) {
  stop(paste("Error : Missing columns in the sampleMetadata File :", paste(missing_cols, collapse = ", ")))
}

required_cols_variable <- c("compoundname")
missing_cols <- setdiff(required_cols_variable, names(variableMetadata))
if (length(missing_cols) > 0) {
  stop(paste("Error : Missing columns in the sampleMetadata File :", paste(missing_cols, collapse = ", ")))
}

# Check if the value have the good format
numeric_cols <- c(required_cols_data)
non_numeric <- numeric_cols[!sapply(dataMatrix[numeric_cols], is.numeric)]
if (length(non_numeric) > 0) {
  stop(paste("The following columns msut be numeric :", paste(non_numeric, collapse = ", ")))
}


#### ---- Create data ----
dataMatrix_t <- as.data.frame(t(dataMatrix[-1]))
colnames(dataMatrix_t) <- dataMatrix[[1]]         
dataMatrix_t$sampleMetadata <- rownames(dataMatrix_t)
dataMatrix_t <- dataMatrix_t[, c("sampleMetadata", setdiff(names(dataMatrix_t), "sampleMetadata"))]
pool_s <- merge(sampleMetadata, dataMatrix_t, by = "sampleMetadata")
pool_s <- pool_s[order(pool_s$InjectionOrder), ]

# Use global injection order or not
if (!opt$global) {
  pool_s$InjectionOrder <- ave(
    pool_s$InjectionOrder,
    pool_s$Batch,
    FUN = function(x) seq_along(x)
  )
}

#### ---- Function to plot the convex hull ----
plot_convex_hull <- function(data, metabolite_name, output_file, show_points) {
  
  # Base plot
  plot <- ggplot(data, aes(x = .data$InjectionOrder,
                           y = .data$Intensity1,
                           color = .data$Batch)) +
    labs(title = paste("Convex Hull by Batch for", metabolite_name),
         x = "Injection Order", y = "Intensity") +
    theme_minimal()
  
  # Add points if requested
  if (show_points) {
    plot <- plot + geom_point(size = 2)
  }
  
  # Add convex hulls
  plot <- plot + geom_polygon(
    data = do.call(rbind, by(data, data$Batch, function(batch) {
      hull_indices <- chull(batch$InjectionOrder, batch$Intensity1)
      hull_indices <- c(hull_indices, hull_indices[1])
      batch[hull_indices, ]
    })),
    aes(x = .data$InjectionOrder,
        y = .data$Intensity1,
        fill = .data$Batch),
    alpha = 0.2
  )

  ## ---- Compute intraD ----
  intra_distances <- data %>%
    group_by(Batch) %>%
    summarise(
      intraD = if (n() > 1) {
        mean(dist(cbind(InjectionOrder, Intensity1)))
      } else {
        NA
      }
    )
  
  intraD <- median(intra_distances$intraD, na.rm = TRUE)
  
  ## ---- Compute interD ----
  centroids <- data %>%
    group_by(Batch) %>%
    summarise(
      x = mean(InjectionOrder),
      y = mean(Intensity1)
    )
  
  if (nrow(centroids) > 1) {
    interD <- mean(dist(cbind(centroids$x, centroids$y)))
  } else {
    interD <- NA
  }
  
  ## ---- Compute Ratio ----
  ratio <- intraD / (1 + interD)
  
  # ---- Format indicators for annotation ----
  indicator_text <- paste0(
    "intraD: ", round(intraD, 3), "\n",
    "interD: ", round(interD, 3), "\n",
    "Ratio: ", round(ratio, 3)
  )
  
  # Add annotation with the indicators
  plot <- plot + annotate("text",
                          x = Inf, y = -Inf,
                          hjust = 1.1, vjust = -0.2,
                          label = indicator_text,
                          size = 3.5, color = "black")
  
  # ---- Save plot as PDF ----
  ggsave(filename = output_file, plot = plot, device = "pdf")
}

#### ---- Function for multiple metabolite ----

plot_convex_hull_metabolite <- function(data, injection_order_col, intensity_col, metabolite_name) {
  
  # ---- Compute intraD ----
  intra_distances <- data %>%
    group_by(Batch) %>%
    summarise(
      intraD = if (n() > 1) {
        mean(dist(cbind(.data[[injection_order_col]], .data[[intensity_col]])))
      } else {
        NA
      }
    )
  
  intraD <- median(intra_distances$intraD, na.rm = TRUE)
  
  # ---- Compute interD ----
  centroids <- data %>%
    group_by(Batch) %>%
    summarise(
      x = mean(.data[[injection_order_col]]),
      y = mean(.data[[intensity_col]])
    )
  
  if (nrow(centroids) > 1) {
    interD <- mean(dist(cbind(centroids$x, centroids$y)))
  } else {
    interD <- NA
  }
  
  # ---- Compute Ratio ----
  ratio <- if (!is.na(intraD) && !is.na(interD) && interD != 0) {
    intraD / (1 + interD)
  } else {
    NA
  }
  
  # ---- Format indicator text ----
  indicator_text <- paste0(
    "intraD: ", ifelse(!is.na(intraD), round(intraD, 3), "NA"), "\n",
    "interD: ", ifelse(!is.na(interD), round(interD, 3), "NA"), "\n",
    "Ratio: ", ifelse(!is.na(ratio), round(ratio, 3), "NA")
  )
  
  # ---- Create base plot ----
  plot <- ggplot(data, aes_string(x = injection_order_col, y = intensity_col, color = "Batch")) +
    geom_point(size = 2) +
    geom_polygon(
      data = do.call(rbind, by(data, data$Batch, function(batch) {
        if (nrow(batch) >= 3) {
          hull_indices <- chull(batch[[injection_order_col]], batch[[intensity_col]])
          hull_indices <- c(hull_indices, hull_indices[1])
          batch[hull_indices, ]
        } else {
          batch[NULL, ]
        }
      })),
      aes_string(x = injection_order_col, y = intensity_col, fill = "Batch"),
      alpha = 0.5, inherit.aes = FALSE
    ) +
    labs(title = paste("Convex Hull for", metabolite_name),
         x = "Injection Order", y = "Intensity") +
    annotate("text",
             x = Inf, y = -Inf,
             hjust = 1.1, vjust = -0.2,
             label = indicator_text,
             size = 3.5, color = "black") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(plot)
}

plot_convex_metabolites <- function(data, injection_order_col, intensity_cols, output_pdf = "convex_hulls.pdf") {
  
  pdf(output_pdf, width = 8, height = 6)
  
  for (intensity_col in intensity_cols) {
    p <- plot_convex_hull_metabolite(data, injection_order_col, intensity_col, intensity_col)
    print(p)
  }
  
  dev.off()
  cat("All plots saved in", output_pdf, "\n")
}
#### ---- Call plotting function ----
if (opt$multiple) {
  intensity_cols <- grep("^Intensity", colnames(pool_s), value = TRUE)
  injectionorder = 'InjectionOrder'
  plot_convex_metabolites(pool_s, injectionorder, intensity_cols)
}else{
  plot_convex_hull(pool_s, opt$metabolite, opt$output, opt$points)
  if (!file.exists(opt$output)) {
    stop("Plot file was not created: ", opt$output)
  }
  
  cat("Plot saved as", opt$output, "\n")
}


