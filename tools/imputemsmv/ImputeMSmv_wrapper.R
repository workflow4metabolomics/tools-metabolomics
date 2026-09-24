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
# Import libraries
library("W4MRUtils")
# ========================================================.


# ====================== Load the tool functions =========.
# Resolved against the script's own directory by W4MRUtils, the same way
# w4mconcatenate does it. Hand-rolling this -- recovering --tool_dir from the
# command line and sourcing the result -- duplicated source_local(), and put
# the wrapper's install path into para, which print(para) then echoed to the
# job log as if it were a tool parameter. A missing file now errors from
# source() itself, naming the path it tried.
source_local("ImputeMSmv_function.R")
# ========================================================.


# ====================== Get XML parameters ==============.
para <- W4MRUtils::parse_args(args = commandArgs())

cat(
    "Job starting time:", format(Sys.time(), "%a %d %b %Y %X"), "\n", strrep("-", 60),
    "\nWelcome to MissMDA, a Galaxy tool dedicated to the handling of missing values from Mass Spectrometry data. \n\n",
    "\nParameters used by this tool:\n"
)
cat("Paramters set: \n")
print(para)
cat("\n", strrep("-", 60), "\n")
# ========================================================.


# ====================== Get tables ======================.
# Retrieve parameters
raw_data <- W4MRUtils::import3(
    pathDM = para$DM_in,
    pathSM = para$SM_in,
    pathVM = para$VM_in
)
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
# All preconditions are checked here, before any computation starts, so that a
# bad job fails fast instead of half-way through the imputation.
validate_inputs <- function(d, s, v, col, LimNan, percLOD, pmin) {
    problems <- character(0)
    require_ <- function(cond, msg) {
        if (!isTRUE(cond)) problems <<- c(problems, msg)
    }

    require_(nrow(d) > 0, "The data matrix contains no features.")
    require_(ncol(d) > 0, "The data matrix contains no samples.")
    require_(nrow(v) == nrow(d),
              "The data matrix and the variable metadata describe a different number of features.")
    require_(nrow(s) == ncol(d),
              "The data matrix and the sample metadata describe a different number of samples.")
    require_(setequal(rownames(v), rownames(d)),
              "Feature identifiers differ between the data matrix and the variable metadata.")
    require_(setequal(rownames(s), colnames(d)),
              "Sample identifiers differ between the data matrix and the sample metadata.")
    require_(all(vapply(d, is.numeric, logical(1))),
              "The data matrix contains non-numeric values.")
    require_(col %in% colnames(s),
              paste0("Column '", col, "' is not found in the sample metadata."))
    require_(!any(is.na(s[[col]])),
              paste0("Column '", col, "' contains missing values."))

    require_(is.numeric(LimNan) && length(LimNan) == 1L && !is.na(LimNan) &&
                 LimNan >= 0 && LimNan <= 1,
              paste0("The missing value rate threshold must be a single number in [0, 1]; got: ",
                     format(LimNan)))

    # percLOD / pmin arrive only on their own MNAR branch, so NULL is legitimate.
    # percLOD = 0 is allowed on purpose: quantile(x, 0) is the global minimum,
    # which is a defensible LOD.
    if (!is.null(percLOD)) {
        require_(is.numeric(percLOD) && length(percLOD) == 1L && !is.na(percLOD) &&
                     percLOD >= 0 && percLOD <= 1,
                  paste0("The LOD quantile must be a single number in [0, 1]; got: ",
                         format(percLOD)))
    }
    if (!is.null(pmin)) {
        require_(is.numeric(pmin) && length(pmin) == 1L && !is.na(pmin) &&
                     pmin > 0 && pmin <= 1,
                  paste0("The fraction of the lowest detected intensity must be a single ",
                         "number in (0, 1]; got: ", format(pmin),
                         " (imputing literal zeros is what this tool exists to avoid)."))
    }

    if (length(problems) > 0) {
        stop(paste(problems, collapse = "\n"), call. = FALSE)
    }
}

# parse_args applies a numeric converter to every value, so a grouping column
# literally named "1" would arrive as the number 1 and later resolve by
# position. Force the identifier back to a name before it is used.
colGroup <- as.character(para$col)

# Convert NaN to NA (must happen before validation: NaN counts as missing)
if (any(sapply(d, function(x) any(is.nan(x))))) {
    cat("NaN values detected in the input data matrix. They have been converted to NA.\n")
    for (i in seq_along(d)) {
        d[[i]][is.nan(d[[i]])] <- NA
    }
}

LimNan <- suppressWarnings(as.numeric(para$LimNan))
percLOD <- if (is.null(para$percLOD)) NULL else suppressWarnings(as.numeric(para$percLOD))
pmin <- if (is.null(para$pmin)) NULL else suppressWarnings(as.numeric(para$pmin))

validate_inputs(d, s, v, colGroup, LimNan, percLOD, pmin)

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

# Grouping column presence/validity is already enforced by validate_inputs()
cat("Column name check passed. \n")
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
if (para$mnar_methods == "lod") {
    cat(" - Missing values categorized as MNAR will be imputed using the limit of detection found in the global datamatrix.\n")
    if (para$mnar_computation == "quantile") {
        cat(" - LOD will be computed using the user-defined quantile of the global data matrix.\n")
    } else {
        cat(" - LOD will be computed using the", para$mnar_computation, "of values below the user-defined LOD.\n")
    }
} else if (para$mnar_methods == "min") {
    cat(" - Missing values categorized as MNAR will be imputed using", ((para$pmin) * 100), "percent of the lowest value computed per feature across all samples.\n")
}
# ========================================================.


# ====================== If Seed =========================.
# The seed itself is applied inside imputemsmv(), which owns the random draws;
# this block only reports what will happen.
if (isTRUE(para$setSeed)) {
    cat(" - A seed has been set to:", para$seed, "(Seeds are used for reproducibility)\n")
} else {
    cat(" - No seed has been manually set by the user; one will be drawn and logged by the tool, and can be reused for reproducibility.\n")
}
# ========================================================.


# ====================== Run imputation ==================.
result_tables <- imputemsmv(
    DataMatrix = d, sampleMetadata = s, variableMetadata = v,
    colClass = colGroup, LimNan = LimNan, setSeed = para$setSeed, seed = para$seed,
    mnar_method = para$mnar_methods, percLOD = percLOD, SDLOD = para$SDLOD, pmin = pmin,
    mar_method = para$mar_methods, SDmean = para$SDmean, mnar_computation = para$mnar_computation,
    advanced = para$advanced, digits = digits
)

# Post-condition: imputation must leave no gap behind. Without this, a matrix
# that is entirely missing runs to completion and Galaxy reports a green
# dataset full of NAs.
if (any(is.na(result_tables$Datamatrix))) {
    still_na <- rownames(result_tables$Datamatrix)[
        rowSums(is.na(result_tables$Datamatrix)) > 0
    ]
    stop(
        "Imputation finished but ", sum(is.na(result_tables$Datamatrix)),
        " value(s) are still missing across ", length(still_na),
        " feature(s), so no reference value could be found for them.\n",
        "Features left unimputed:\n", paste(utils::head(still_na, 20), collapse = ", "),
        if (length(still_na) > 20) " ..." else "",
        call. = FALSE
    )
}
# ========================================================.


# ====================== Get tabular outputs =============.
# Named access: imputemsmv() returns a named list, and [[1]]/[[2]]/[[3]] silently
# follow any reordering of that return value.
# The identifier column is named by position, not by matching the placeholder:
# matching on the name also renamed any sample, feature or metadata field
# genuinely called "df", yielding two columns of the same name.
DM_out <- data.frame(id = rownames(d), result_tables$Datamatrix, check.names = FALSE)
colnames(DM_out)[1] <- dName
# SM
SM_out <- data.frame(id = rownames(s), result_tables$sampleMetadata, check.names = FALSE)
colnames(SM_out)[1] <- sName
# VM
VM_out <- data.frame(id = rownames(v), result_tables$VariableMetadata, check.names = FALSE)
colnames(VM_out)[1] <- vName

# Export output tables
# quote = FALSE for all three, as the SM/VM writes already did and as
# w4mconcatenate does: quoting left identifiers as "M60T53" in the data matrix.
write.table(DM_out, file = para$DM_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(SM_out, file = para$SM_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(VM_out, file = para$VM_out, sep = "\t", row.names = FALSE, quote = FALSE)
# ========================================================.


# ====================== Plots ===========================.
# Plotting lives in a function so the device close is a real handler rather than
# dead code (C13): a top-level on.exit() registers into R_GlobalContext, which
# Rscript's exit path never unwinds, so the registration at the end of the old
# block could only ever have double-closed the explicit dev.off() below it -- it
# ran in neither case. Inside a function the handler does fire on the error path.
# Measured with graphics::image patched to fail mid-plot: old and new both exit 1
# and leave a well-formed PDF holding the pages completed so far, because
# grDevices finalises the device at process exit. So this is hygiene -- cleanup
# that no longer depends on that finalisation -- not a fix for an observable bug.
make_plots <- function(chosen_plots) {
    pdf(file = para$plot_out, width = 11, height = 7)
    on.exit(dev.off(), add = TRUE)

    # "None" was accepted here but is not one of the tool's option values, so the
    # branch was unreachable and an empty selection fell through to the plotting
    # path, which logged "plotting 0 chart(s)" and left Galaxy holding a valid
    # 0-page PDF (C5). Treat "None" as what it always meant -- no plot -- and say
    # so on the page, so the PDF always has at least one readable page.
    chosen_plots <- setdiff(chosen_plots, c("None", ""))

    if (length(chosen_plots) == 0) {
        plot_placeholder("NO PLOT SELECTED")
        return(invisible(NULL))
    }

    cat("\n", strrep("-", 60), "\n")
    cat("\nGraphical representations: plotting", length(chosen_plots), "chart(s):\n")

    if (para$mnar_methods == "min") {
        cat("\nPlotting distribution first:\n")
        plt_distribution(d, 16)
    }

    if ("DP" %in% chosen_plots) {
        cat("\nDensity plot")
        # the tidy matrix, not DM_out: passing the output table meant the "after"
        # curve depended on whether column 1 happened to hold identifiers (R25)
        plt_density(d, result_tables$Datamatrix)
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

make_plots(as.vector(strsplit(para$plots, ",")[[1]]))
# ========================================================.


# ====================== Script finished =================.
cat("\n", strrep("-", 60), "\nInformation about R (version, Operating System, attached or loaded packages):\n")
sessionInfo()
cat("\n", strrep("-", 60), "\nJob finished at:\n", format(Sys.time(), "%a %d %b %Y %X"))
# ========================================================.
