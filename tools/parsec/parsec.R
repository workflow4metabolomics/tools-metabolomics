options(warn = -1)

# --- LIBRAIRIES ---
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(lme4)
})

# --- FONCTION DE CORRECTION ---
batch_cohort_correction <- function(data, batch_col, sample_col, intensity_cols) {
    # Vérifie présence des colonnes
    missing_cols <- setdiff(c(batch_col, sample_col, intensity_cols), colnames(data))
    if (length(missing_cols) > 0) {
        stop(paste("❌ Colonnes manquantes :", paste(missing_cols, collapse = ", ")))
    }

    # 🔧 Nettoyage et conversion en numérique
    data[intensity_cols] <- lapply(data[intensity_cols], function(x) {
        x <- gsub("\\s+", "", as.character(x))
        as.numeric(x)
    })

    # 1. Log-transform
    data <- data %>%
        mutate(across(all_of(intensity_cols), log1p))

    # 2. Standardisation par batch
    data <- data %>%
        group_by(!!sym(batch_col)) %>%
        mutate(across(all_of(intensity_cols), ~ scale(.x)[, 1])) %>%
        ungroup()

    # 3. Vérification Injection_Order
    if (!"Injection_Order" %in% colnames(data)) {
        stop("❌ Colonne Injection_Order manquante.")
    }

    # 4. Modèle mixte linéaire
    for (col in intensity_cols) {
        model <- lmer(
            as.formula(paste0(col, " ~ Injection_Order + (1|", batch_col, ")")),
            data = data,
            REML = TRUE,
            control = lmerControl(check.conv.singular = "ignore")
        )
        data[[col]] <- residuals(model)
    }

    # 5. Inverse transform
    data <- data %>%
        mutate(across(all_of(intensity_cols), expm1))

    return(data)
}

# --- ARGUMENTS CLI ---
option_list <- list(
    make_option(c("-d", "--dataMatrice"), type = "character", help = "Data matrix"),
    make_option(c("-s", "--sampleMData"), type = "character", help = "Sample metadata"),
    make_option(c("-v", "--variableMData"), type = "character", help = "Variable metadata"),
    make_option(c("-o", "--output"), type = "character", help = "Output file")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- CHARGEMENT DES FICHIERS ---
if (!all(file.exists(opt$dataMatrice, opt$sampleMData, opt$variableMData))) {
    stop("❌ Un ou plusieurs fichiers d'entrée sont introuvables.")
}

data_matrix <- read.csv(opt$dataMatrice, header = TRUE, sep = "\t")
sample_metadata <- read.csv(opt$sampleMData, header = TRUE, sep = "\t")
variable_metadata <- read.csv(opt$variableMData, header = TRUE, sep = "\t")

# --- RENOMMAGE DES COLONNES ---
colnames(sample_metadata)[1] <- "SampleID"
colnames(sample_metadata) <- sub("^batch$", "Batch", colnames(sample_metadata))
colnames(sample_metadata) <- sub("^injectionOrder$", "Injection_Order", colnames(sample_metadata))

# --- TRANSFORMATION DES DONNÉES ---
data_t <- as.data.frame(t(data_matrix))
data_t$SampleID <- rownames(data_t)

if ((ncol(data_t) - 1) != nrow(variable_metadata)) {
    stop("❌ Incompatibilité : nombre de variables ≠ nombre de métadonnées.")
}

ion_names <- paste0("Ion", seq_len(nrow(variable_metadata)))
colnames(data_t) <- c(ion_names, "SampleID")

# --- FUSION ---
data_set <- merge(sample_metadata, data_t, by = "SampleID")

# --- COLONNES À CORRIGER ---
intensity_cols <- ion_names

# Vérifie que toutes les colonnes nécessaires sont là
required_columns <- c("Batch", "SampleID", "Injection_Order", intensity_cols[1:2])
missing_columns <- setdiff(required_columns, colnames(data_set))
if (length(missing_columns) > 0) {
    stop(paste("❌ Colonnes manquantes dans le fichier fusionné :", paste(missing_columns, collapse = ", ")))
}

# --- APPLICATION DE LA CORRECTION ---
corrected_data <- batch_cohort_correction(
    data_set,
    batch_col = "Batch",
    sample_col = "SampleID",
    intensity_cols = intensity_cols
)

# --- EXPORT FINAL ---
write.csv(
    corrected_data,
    file = opt$output,
    quote = TRUE,
    row.names = FALSE
)
