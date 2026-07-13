# ========================================================.
# Tool R Wrapper
#
# Original tool wrapping: Anthony SLEIMAN (Internship student, PFEM, INRAE, MetaboHub)
# Maintainer: Melanie PETERA (PFEM, INRAE, MetaboHub)
#
# Creation date: 2026-03-09 - v0.1.0
# Last modification: 2026-05-26 - v1.0.0
#
# Versions:
# version 0.1.0 - 2026-03-09 - Initial version
# Version 1.0.0: version 1.0.0 - 2026-06-29 - addition of new features
# ========================================================.


# ====================== Prepare environment =============.
# Clean R environment
rm(list = ls())

# Import libraries
library("W4MRUtils")
# ========================================================.


# ====================== Establish XML connections =======.
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) == 1 && idx < length(args)) {
        return(args[idx + 1])
    }
    return(NULL)
}

tool_dir <- get_arg("--tool_dir")

if (is.null(tool_dir)) {
    stop("tool_dir not provided")
}

source(file.path(tool_dir, "ImputeMSmv_function.R"))
# ========================================================.


# ====================== Get XML parameters ==============.
para <- W4MRUtils::parse_args(args = commandArgs())

cat("Job starting time:", format(Sys.time(), "%a %d %b %Y %X"), "\n", strrep("-", 60),
    "\nWelcome to MissMDA, a Galaxy tool dedicated to the handling of missing values from Mass Spectrometry data. \n\n",
    "\nParameters used by this tool:\n")
cat("Paramters set: \n")
print(para)
cat("\n", strrep("-", 60), "\n")
# ========================================================.


# ====================== Get tables ======================.
# Retrieve parameters
raw_data <- W4MRUtils::import3(pathDM = para$DM_in,
                                pathSM = para$SM_in,
                                pathVM = para$VM_in)
# Find input data
d <- raw_data$dataMatrix
s <- raw_data$sampleMetadata
v <- raw_data$variableMetadata

# row.names = FALSE by default. Fix and save table name:
# DM
dName <- colnames(d)[1]
d <- fix_head(d)
# SM
sName <- colnames(s)[1]
s <- fix_head(s)
# VM
vName <- colnames(v)[1]
v <- fix_head(v)
# ========================================================.


# ====================== Quality control =================.
# Check that all columns are numeric
if (!all(vapply(d, is.numeric, logical(1)))) {
    warning("Input data matrix contains non-numeric values")
    stop("Please provide a numeric data matrix")
}

# Convert NaN to NA
if (any(sapply(d, function(x) any(is.nan(x))))) {
    cat("NaN values detected in the input data matrix. They have been converted to NA.\n")
    for (i in seq_along(d)) {
        d[[i]][is.nan(d[[i]])] <- NA
    }
}

# Check presence of missing values
# In global data matrix
if (sum(is.na(d)) == 0) {
    cat("\n\nNo missing values found in Datamatrix, nothing to impute. Aborting job... \n\n")
    cat("+-----------------------------------------------------------------------+
|                                                                       |
|                NO MISSING DATA FOUND! ABORTING JOB!                   |
|                                                                       |
+-----------------------------------------------------------------------+")
    stop("JOB TERMINATED - CHECK RAW DATA MATRIX - NO MISSING DATA FOUND")
} else {
    cat("Missing values found. Job running with the following parameters:\n")
}
# Per feature (row)
sm_new <- rbind(apply(d, 2, function(x) sum(is.na(x))))
row_na <- rownames(sm_new)[apply(sm_new == 0, 1, all)]
if (length(row_na) != 0) {
    stop("One or more features contain 100% of missing values (missing values across all groups), please check your data.\nThe concerned variables are the following:\n", row_na)
}

# Check the number of decimal digits
DIGITS <- read.table(para$DM_in, header = TRUE, row.names = 1, dec = ".", sep = "\t", colClasses = "character")
stripped_nb <- sapply(strsplit(unlist(DIGITS), "\\."), "[", 2)
char_nb <- nchar(stripped_nb)
digits <- suppressWarnings(max(char_nb, na.rm = TRUE))
if (digits == -Inf) {
    digits <- 0
} else {

}
cat("The maximal number of decimal digits found in the data matrix is:", digits, "\nTherefore, imputation will be rounded to", digits, "digits after the decimal separator.\n\n")

# Check features' order in datamatrix and sample metadata
if (identical(rownames(s), colnames(d))) {
    cat("Variables order check passed.\n")
} else {
    cat("Tables not in the same order . . . ")
    d <- d[, sort(names(d))]
    s <- s[order(rownames(s)), ]
    cat("Sorting done.\n")
}

# Check column name entered by user
if (para$col %in% colnames(s)) {
    cat("Column name check passed. \n")
} else {
    stop("Column name: '", para$col, "' is not found in the sample metadata file.\nMake sure that the entered name matches the column name in the sample metadata.\n")
}
# ========================================================.


# ====================== If Advanced =====================.
if (para$advanced == TRUE) {
    cat("Job running in ADVANCED mode (please remember this is NOT THE RECOMMENDED mode!) with the following parameters:\n")
} else {
    cat("Job running in Recommended mode, with the following parameters:\n - Missing values' threshold is set to: ", para$LimNan, "\n")
}
# ========================================================.


# ====================== If MAR_Method ===================.
cat(" - Imputing missing values of the MAR category with the", para$mar_methods, "of group replicates.\n")
# ========================================================.


# ====================== If MNAR_Method ==================.
if (para$mnar_method == "lod") {
    cat(" - Missing values categorized as MNAR will be imputed using the limit of detection found in the global datamatrix.\n")
    if (para$mnar_computation == "quantile") {
        cat(" - LOD will be computed using the user-defined quantile of the global data matrix.\n")
    } else {
        cat(" - LOD will be computed using the", para$mnar_computation, "of values below the user-defined LOD.\n")
    }
} else if (para$mnar_method == "min") {
    cat(" - Missing values categorized as MNAR will be imputed using", ((para$pmin) * 100), "percent of the lowest value computed per feature across all samples.\n")
}
# ========================================================.


# ====================== If Seed =========================.
if (para$setSeed == TRUE) {
    cat(" - A seed has been set to:", para$seed, "(Seeds are used for reproductibility)\n")
    set.seed(as.numeric(para$seed))
} else {
    cat(" - No seed has been manually set by the user. Consequently, a seed will be randomly assigned, that can be reused for reproducibility. ")
}
# ========================================================.


# ====================== Run imputation ==================.
result_tables <- imputemsmv(DataMatrix = d, sampleMetadata = s, variableMetadata = v,
                            colClass = para$col, LimNan = para$LimNan, setSeed = para$setSeed, seed = para$seed,
                            mnar_method = para$mnar_methods, percLOD = para$percLOD, SDLOD = para$SDLOD, pmin = para$pmin,
                            mar_method = para$mar_methods, SDmean = para$SDmean, mnar_computation = para$mnar_computation,
                            advanced = para$advanced, digits = digits)
# ========================================================.


# ====================== Get tabular outputs =============.
# DM
DM_out <- data.frame(df = rownames(d), result_tables[[1]], check.names = FALSE)
colnames(DM_out)[colnames(DM_out) == "df"] <- dName
# SM
SM_out <- data.frame(df = rownames(s), result_tables[[2]], check.names = FALSE)
colnames(SM_out)[colnames(SM_out) == "df"] <- sName
# VM
VM_out <- data.frame(df = rownames(v), result_tables[[3]], check.names = FALSE)
colnames(VM_out)[colnames(VM_out) == "df"] <- vName

# Export output tables
write.table(DM_out, file = para$DM_out, sep = "\t", row.names = FALSE)
write.table(SM_out, file = para$SM_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(VM_out, file = para$VM_out, sep = "\t", row.names = FALSE, quote = FALSE)
# ========================================================.


# ====================== Plots ===========================.
chosen_plots <- as.vector(strsplit(para$plots, ",")[[1]])

pdf(file = para$plot_out, width = 11, height = 7)

on.exit(dev.off(), add = TRUE)

if ("None" %in% chosen_plots) {
    plot.new()
    text(0.5, 0.5, "NO PLOT SELECTED")

} else {
    cat("\n", strrep("-", 60), "\n")
    cat("\nGraphical representations: plotting", length(chosen_plots), "chart(s):\n")

    if (para$mnar_methods == "min") {
        cat("\nPlotting distribution first:\n")
        plt_distribution(d, 16)
    }

    if ("DP" %in% chosen_plots) {
        cat("\nDensity plot")
        plt_density(d, DM_out)
    }

    if ("HM" %in% chosen_plots) {
        cat("\nHeatmap")
        plt_matrix(d)
    }

    if ("BP1" %in% chosen_plots) {
        cat("\nGroup barplot")
        plt_group(SM_out, VM_out, para$col)
    }

    if ("BP2" %in% chosen_plots) {
        cat("\nType barplot")
        plt_type(SM_out, VM_out, para$col)
    }

    if ("PC" %in% chosen_plots) {
        cat("\nPie chart\n")
        plt_pie(SM_out, VM_out, para$col)
    }
}

dev.off()
# ========================================================.


# ====================== Script finished =================.
cat("\n", strrep("-", 60), "\nInformation about R (version, Operating System, attached or loaded packages):\n")
sessionInfo()
cat("\n", strrep("-", 60), "\nJob finished at:\n", format(Sys.time(), "%a %d %b %Y %X"))
# ========================================================.
