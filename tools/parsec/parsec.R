options(warn = -1)
# --- IMPORT PACKAGES ---
suppressPackageStartupMessages({
    library(optparse)
    library(lme4)
})

# --- MAIN FUNCTION (CORRECTION) ---
batch_cohort_correction <- function(
    data,
    batch_col,
    injection_order_col,
    intensity_cols) {
    # Checks for mandatory columns
    missing_cols <- setdiff(
        c(batch_col, intensity_cols, injection_order_col),
        colnames(data)
    )
    if (length(missing_cols) > 0) {
        stop(
            paste(
                "❌  missing columns",
                paste(missing_cols, collapse = ", ")
            )
        )
    }

    # 🔧 Cleaning the data and conversion to numeric values
    data[intensity_cols] <- lapply(data[intensity_cols], function(x) {
        x <- gsub("\\s+", "", as.character(x))
        as.numeric(x)
    })

    # 1. Log-transform
    data[intensity_cols] <- lapply(data[intensity_cols], log1p)

    # 2. Batch standardisation
    data <- unsplit(
        lapply(
            split(data, data[[batch_col]]),
            function(df) {
                df[intensity_cols] <- lapply(
                    df[intensity_cols],
                    function(col) as.numeric(scale(col))
                )
                df
            }
        ),
        data[[batch_col]]
    )

    # 3. Linear mixed model
    for (col in intensity_cols) {
        model <- lmer(
            as.formula(paste0(
                col,
                " ~ ", injection_order_col, " + (1|", batch_col, ")"
            )),
            data = data,
            REML = TRUE,
            control = lmerControl(check.conv.singular = "ignore")
        )
        data[[col]] <- residuals(model)
    }

    # 4. Inverse transform
    data[intensity_cols] <- lapply(data[intensity_cols], expm1)

    return(data)
}

# --- CLI ARGUMENTS ---
option_list <- list(
    make_option(
        c("-d", "--dataMatrix"),
        type = "character",
        help = "Data matrix"
    ),
    make_option(
        c("-s", "--sampleMData"),
        type = "character",
        help = "Sample metadata"
    ),
    make_option(
        c("-v", "--variableMData"),
        type = "character",
        help = "Variable metadata"
    ),
    make_option(
        c("-o", "--output"),
        type = "character",
        help = "Output file"
    )
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- FILE LOADINGS ---
if (!all(file.exists(opt$dataMatrice, opt$sampleMData, opt$variableMData))) {
    stop("❌ At least one of the input files could not be found !")
}

data_matrix <- read.csv(
    opt$dataMatrice,
    header = TRUE,
    sep = "\t",
    row.names = 1
)
sample_metadata <- read.csv(
    opt$sampleMData,
    header = TRUE,
    sep = "\t",
    row.names = 1
)
variable_metadata <- read.csv(
    opt$variableMData,
    header = TRUE,
    sep = "\t",
    row.names = 1
)


# --- DATA TRANSPOSE ---
data_t <- as.data.frame(t(data_matrix))

if (ncol(data_t) != nrow(variable_metadata)) {
    stop("❌ Incompatibility: nb of variables ≠ nb of variableMetadata")
}

# --- COLUMNS SUBJECTED TO CORRECTION ---
intensity_cols <- rownames(variable_metadata)

# --- DM+SM FUSION ---
data_set <- transform(
    merge(sample_metadata, data_t, by = 0),
    row.names = Row.names,
    Row.names = NULL
)

# --- CORRECTION---
corrected_data <- batch_cohort_correction(
    data_set,
    batch_col = "batch",
    injection_order_col = "injectionOrder",
    intensity_cols = intensity_cols
)

# --- FINAL EXPORT ---
write.table(
    round(corrected_data[intensity_cols], 10),
    file = opt$output,
    quote = TRUE,
    row.names = TRUE,
    sep = "\t"
)
