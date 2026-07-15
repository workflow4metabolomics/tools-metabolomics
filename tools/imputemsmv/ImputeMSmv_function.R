# ========================================================.
# Tool: ImputeMsmv
#
# Original tool code: Mathieu CLADIERE (Université Paris Saclay, INRAE, AgroParisTech, UMR Sayfood, Campus Agro Paris Saclay - Palaiseau, France)
# Original tool wrapping: Anthony SLEIMAN (Internship student, PFEM, INRAE, MetaboHub)
# Maintainer: Melanie PETERA (PFEM, INRAE, MetaboHub, W4M)
#
# Creation date: 2026-03-09 - v0.1.0
# Last modification: 2026-05-26 - v1.0.0
#
# Versions:
# version 0.1.0 - 2026-03-09 - Initial version
# Version 1.0.0: version 1.0.0 - 2026-06-29 - addition of new features
#
# Description:
# This tool performs imputation of missing values in mass spectrometry–based data
# matrices. The primary objective is to handle incomplete metabolomic datasets by
# estimating and replacing missing intensities in a statistically informed manner.
#
# Libraries and dependencies:
# This is an R script
# No external dependencies yet
#
# Input:
# 3 tabular files in standard W4M format:
# Datamatrix containing missing values in NA format
# Sample Metadata containing a column describing samples from same group
# Variable Metadata
#
# Output:
# Datamatrix with no missing values (all imputed)
# Sample Metadata with additional columns describing NA category and proportions in samples
# Variable Metadata with additionam columns describing NA category and proportions in variables
# Graphics and plots in a pdf format
# ========================================================.


# ================= imputemsmv ===========================.
#' Main tool function imputemsmv
#'
#' This function performs missing value analysis and imputation on a
#' feature-by-sample intensity matrix. It first classifies missing values into three categories (No_MV, MAR, MNAR)
#' at the feature level within each biological group, then imputes missing values using different strategies
#' depending on their classification.
#' No_MV: no imputation required
#' MAR: imputed using group-wise replicates mean or median with multiplicative noise
#' MNAR: imputed using either LOD-based or minimum-value-based strategy
#'
#' @param DataMatrix Numeric matrix (features x samples) containing intensities with NA values
#' @param sampleMetadata Data frame containing sample annotations and grouping variable
#' @param variableMetadata Data frame containing feature annotations
#' @param colClass Column name in sampleMetadata defining sample's group
#' @param LimNan Numeric (0–1); threshold of missingness per feature within a group above which values are considered MNAR
#' @param setSeed Whether to set random seed for reproducibility
#' @param seed Random seed value
#' @param mnar_method Method for MNAR imputation ("lod" or "min")
#' @param percLOD Numeric (0–1); quantile used to estimate limit of detection (LOD)
#' @param SDLOD Numeric (0–1); multiplicative noise amplitude for MNAR imputation
#' @param pmin Numeric (0–1); scaling factor applied in "min" MNAR strategy
#' @param mar_method Method for MAR imputation ("mean" or "median")
#' @param SDmean Numeric (0–1); multiplicative noise amplitude for MAR imputation
#' @param advanced Modifies MNAR/MAR classification threshold behavior
#' @param digits Rounding precision applied during imputation
#'
#' @returns A list containing 3 W4M-compatible tables with no missing values: Datamatrix sampleMetadata VariableMetadata
#'
#' @return Imputed intensity data matrix with no NA values
#' @return Original sample metadata with added MV percentage per sample
#' @return Original variable metadata with added missing value classification per group
#'
#' @author
#' Original code: Mathieu CLADIERE
#' Function and minor modifications: Anthony SLEIMAN
#' Maintainer: Mélanie PETERA
imputemsmv <- function(DataMatrix,
                       sampleMetadata,
                       variableMetadata,
                       colClass,
                       LimNan,
                       setSeed,
                       seed,
                       mnar_method = c("lod", "min"),
                       percLOD,
                       SDLOD,
                       pmin,
                       mar_method = c("mean", "median"),
                       SDmean,
                       mnar_computation = c("quantile", "mean", "median"),
                       advanced,
                       digits) {
    # Match user-selected methods
    mnar_method <- match.arg(mnar_method)
    mar_method <- match.arg(mar_method)
    mnar_computation <- match.arg(mnar_computation)

    # Optional reproductibility
    if (setSeed == TRUE) {
        set.seed(seed)
    } else {
        random_seed <- sample(1:10000, 1)
        cat("Seed number assigned value:", random_seed, ".\n")
    }

    cat(strrep("-", 60), "\nCounting features: \n")

    # Create feature index to preserve original order after processing
    numero <- seq_len(nrow(DataMatrix))
    if (length(numero) == 0) {
        stop("No features found in Data matrix. Aborting job.")
    }
    ordre <- data.frame(numero = numero, row.names = row.names(DataMatrix))

    # Extract unique group labels from sample metadata
    tUclass <- t(unique(sampleMetadata[colClass]))

    # Global missing valu statitics
    nMV <- sum(is.na(DataMatrix))
    pMV <- (nMV / (nrow(DataMatrix) * ncol(DataMatrix))) * 100
    cat("In raw data matrix there are: \n .", length(numero), "feature(s) \n .", length(tUclass), "group(s)")
    if (length(tUclass[1, ]) < 12) {
        cat("-", tUclass[1, ])
    }
    cat("\n .", (nrow(DataMatrix) * ncol(DataMatrix)), "values (missing included) \n .", nMV, "missing values:", round(pMV, digits = 2), "% of total values \n")
    cat(strrep("-", 60))
    cat("\n")

    # Compute missing values percentage per sample (used later in output metadata)
    SM_new <- t(rbind(apply(DataMatrix, 2, function(x) round(((sum(is.na(x)) / length(numero)) * 100), 2))))
    colnames(SM_new) <- c("MV_percentage")

    # Prepare vector containers for per-group data
    nom_sampleMedata_gp <- seq_along(tUclass)
    nom_DataMatrix_gp <- seq_along(tUclass)

    # Split data matrix and sample metadata by group (later, each group is imputed apart)
    for (i in seq_along(tUclass)) {
        nom_sampleMedata_gp[i] <- paste("SampleMetadata_", tUclass[i], sep = "")
        assign(nom_sampleMedata_gp[i], subset(sampleMetadata, get(colClass) == tUclass[i]))
        nom_DataMatrix_gp[i] <- paste("DataMatrix_", tUclass[i], sep = "")
        assign(nom_DataMatrix_gp[i], cbind(ordre, subset(DataMatrix, select = rownames(get(nom_sampleMedata_gp[i])))))
    }

    # Missing Value management ##############################.

    # Prepare vector containers for per-group data
    nom_DataMatrix_gpMV <- seq_along(nom_DataMatrix_gp)
    nom_DataMatrix_gpMARMNAR <- seq_along(nom_DataMatrix_gp)

    # Count MV for each group and calculate mean
    group <- c()
    replicate_count <- c()
    No_MV_count <- c()
    MAR_count <- c()
    MNAR_count <- c()
    for (i in seq_along(nom_DataMatrix_gp)) {
        # Count NA per feature within group
        NbNaN <- apply(get(nom_DataMatrix_gp[i]), 1, function(x) sum(is.na(x)))
        pNbNaN <- NbNaN / (ncol(get(nom_DataMatrix_gp[i])) - 1)

        # Store NA counts and their proportions
        nom_DataMatrix_gpMV[i] <- paste(nom_DataMatrix_gp[i], "_MV", sep = "")
        assign(nom_DataMatrix_gpMV[i], cbind(get(nom_DataMatrix_gp[i]), "NbNaN" = NbNaN, "pNbNaN" = pNbNaN))

        # Classify missingness type per feature
        MARMNAR <- seq_along(numero)
        rep_num <- ncol(get(nom_DataMatrix_gp[i])) - 1
        if (advanced == TRUE) {
            repTH <- 2
        } else {
            repTH <- 3
        }
        if (rep_num < repTH) {
            for (j in seq_along(numero)) {
                if (pNbNaN[j] == 0) {
                    MARMNAR[j] <- "No_MV"
                } else {
                    MARMNAR[j] <- "MNAR"
                }
            }
        } else {
            for (j in seq_along(numero)) {
                if (pNbNaN[j] == 0) {
                    MARMNAR[j] <- "No_MV"
                } else if (pNbNaN[j] <= LimNan) {
                    MARMNAR[j] <- "MAR"
                } else {
                    MARMNAR[j] <- "MNAR"
                }
            }
        }
        group[i] <- tUclass[i]
        replicate_count[i] <- rep_num
        No_MV_count[i] <- sum(MARMNAR == "No_MV")
        MAR_count[i] <- sum(MARMNAR == "MAR")
        MNAR_count[i] <- sum(MARMNAR == "MNAR")

        nom_DataMatrix_gpMARMNAR[i] <- paste(nom_DataMatrix_gp[i], "_MARMNAR", sep = "")
        assign(nom_DataMatrix_gpMARMNAR[i], MARMNAR)
    }

    informations <- data.frame(group, replicate_count, No_MV_count, MAR_count, MNAR_count)
    cat("\nPer group count summary:\n")
    if (nrow(informations) < 12) {
        cat("\n")
        print(informations)
    } else {
        cat("(main groups extract)\n\n")
        informations <- informations[order(informations[, 2], decreasing = TRUE), ]
        print(informations[1:12, ])
    }

    # Combine classification across groups into variable metadata as final output
    VariableData_MARMNAR <- cbind(as.data.frame(mget(nom_DataMatrix_gpMARMNAR[seq(nom_DataMatrix_gp)])))
    tUclassName <- paste0("mv_type_", tUclass)
    colnames(VariableData_MARMNAR) <- tUclassName

    # MAR imputation ########################################.

    # Prepare vector containers for per-group data
    nom_DataMatrix_gp_MAR <- seq_along(nom_DataMatrix_gp)
    nom_DataMatrix_gp_noMAR <- seq_along(nom_DataMatrix_gp)

    for (i in seq_along(nom_DataMatrix_gp)) {
        # Select MAR features
        nom_DataMatrix_gp_MAR[i] <- paste(nom_DataMatrix_gp[i], "_MAR", sep = "")
        assign(nom_DataMatrix_gp_MAR[i], subset(get(nom_DataMatrix_gpMV[i]), pNbNaN <= LimNan))
        Mat_noMAR_temp <- t(get(nom_DataMatrix_gp_MAR[i]))
        Mat_noMar_temp_2 <- Mat_noMAR_temp[!(row.names(Mat_noMAR_temp) %in% c("NbNaN", "pNbNaN")), ]

        # Impute according to chosen method (mean or median)
        if (mar_method == "mean") {
            for (k in seq_len(nrow(Mat_noMar_temp_2))) {
                for (j in seq_len(ncol(Mat_noMar_temp_2))) {
                    Mat_noMar_temp_2[k, j][is.na(Mat_noMar_temp_2[k, j])] <- round((mean(Mat_noMar_temp_2[seq(2, nrow(Mat_noMar_temp_2)), j], na.rm = TRUE) * round(runif(1, min = (1 - SDmean), max = (1 + SDmean)), digits = 2)), digits)
                }
            }
        } else if (mar_method == "median") {
            for (k in seq_len(nrow(Mat_noMar_temp_2))) {
                for (j in seq_len(ncol(Mat_noMar_temp_2))) {
                    Mat_noMar_temp_2[k, j][is.na(Mat_noMar_temp_2[k, j])] <- round((median(Mat_noMar_temp_2[seq(2, nrow(Mat_noMar_temp_2)), j], na.rm = TRUE) * round(runif(1, min = (1 - SDmean), max = (1 + SDmean)), digits = 2)), digits)
                }
            }
        } else {
            stop("Non valid method, choose one method between mean or median")
        }
        nom_DataMatrix_gp_noMAR[i] <- paste(nom_DataMatrix_gp[i], "_noMAR", sep = "")
        assign(nom_DataMatrix_gp_noMAR[i], t(Mat_noMar_temp_2))
    }

    # MNAR imputation ######################################.

    # Prepare vector containers for per-group data
    nom_DataMatrix_gp_MNAR <- seq_along(nom_DataMatrix_gp)
    nom_DataMatrix_gp_noMNAR <- seq_along(nom_DataMatrix_gp)

    # Impute according to chosen method (LOD or minimum)
    if (mnar_method == "lod") {
        # Compute global LOD from quantile
        Lod <- quantile(DataMatrix, na.rm = TRUE, prob = percLOD, names = FALSE)
        if (mnar_computation != "quantile") {
            low_lod <- na.omit(unlist(DataMatrix[DataMatrix <= Lod]))
            if (mnar_computation == "mean") {
                Lod <- mean(low_lod)
            } else if (mnar_computation == "median") {
                Lod <- median(low_lod)
            } else {
                stop("Non valid LOD computation method.\nAborting job.")
            }
        }

        cat("\n")
        cat(strrep("-", 60))
        cat("\n")
        cat("The data-driven limit of detection computed from data matrix is: ", Lod, "\nIn log(10), LOD is: ", log10(Lod))
        for (i in seq_along(nom_DataMatrix_gp)) {
            # Select MNAR features
            nom_DataMatrix_gp_MNAR[i] <- paste(nom_DataMatrix_gp[i], "_MNAR", sep = "")
            assign(nom_DataMatrix_gp_MNAR[i], subset(get(nom_DataMatrix_gpMV[i]), pNbNaN > LimNan))
            Mat_noMNAR_temp <- get(nom_DataMatrix_gp_MNAR[i])

            # Impute with LOD + noise
            if (nrow(get(nom_DataMatrix_gp_MNAR[i])) > 0) {
                for (j in seq_len(nrow(Mat_noMNAR_temp))) {
                    for (k in seq_len(ncol(Mat_noMNAR_temp))) {
                        Mat_noMNAR_temp[j, k][is.na(Mat_noMNAR_temp[j, k])] <- round((Lod * round(runif(1, min = (1 - SDLOD), max = (1 + SDLOD)), digits = 2)), digits)
                    }
                }
            }
            Mat_noMNAR_temp_2 <- Mat_noMNAR_temp[, !names(Mat_noMNAR_temp) %in% c("NbNaN", "pNbNaN")]
            nom_DataMatrix_gp_noMNAR[i] <- paste(nom_DataMatrix_gp[i], "_noMNAR", sep = "")
            assign(nom_DataMatrix_gp_noMNAR[i], Mat_noMNAR_temp_2)
        }
        # Minimum-based per feature MNAR imputation
    } else if (mnar_method == "min") {
        for (i in seq_along(nom_DataMatrix_gp)) {
            nom_DataMatrix_gp_MNAR[i] <- paste(nom_DataMatrix_gp[i], "_MNAR", sep = "")
            assign(nom_DataMatrix_gp_MNAR[i], subset(get(nom_DataMatrix_gpMV[i]), pNbNaN > LimNan))
            Mat_noMNAR_temp <- get(nom_DataMatrix_gp_MNAR[i])

            if (nrow(get(nom_DataMatrix_gp_MNAR[i])) > 0) {
                for (j in seq_len(nrow(Mat_noMNAR_temp))) {
                    minimal_value <- min(Mat_noMNAR_temp[j, ], na.rm = TRUE)
                    for (k in seq_len(ncol(Mat_noMNAR_temp))) {
                        Mat_noMNAR_temp[j, k][is.na(Mat_noMNAR_temp[j, k])] <- round((pmin * (minimal_value * round(runif(1, min = (1 - SDLOD), max = (1 + SDLOD)), digits = 2))), digits)
                    }
                }
            }
            Mat_noMNAR_temp_2 <- Mat_noMNAR_temp[, !names(Mat_noMNAR_temp) %in% c("NbNaN", "pNbNaN")]
            nom_DataMatrix_gp_noMNAR[i] <- paste(nom_DataMatrix_gp[i], "_noMNAR", sep = "")
            assign(nom_DataMatrix_gp_noMNAR[i], Mat_noMNAR_temp_2)
        }
    } else {
        stop("Non valid method, choose one method between LOD or Min")
    }

    # Reassemble full matrix ################################.

    nom_DataMatrix_gp_noMV <- seq_along(nom_DataMatrix_gp)
    for (i in seq_along(nom_DataMatrix_gp)) {
        # Merge mAR and MNAR
        Mat_NoMV_temp <- rbind(get(nom_DataMatrix_gp_noMAR[i]), get(nom_DataMatrix_gp_noMNAR[i]))
        # Restore original feature order
        Mat_noMV_temp2 <- Mat_NoMV_temp[order(Mat_NoMV_temp$numero), ]
        # Remove index column
        Mat_noMV_temp3 <- Mat_noMV_temp2[, -1]

        nom_DataMatrix_gp_noMV[i] <- paste(nom_DataMatrix_gp[i], "_noMV", sep = "")
        assign(nom_DataMatrix_gp_noMV[i], Mat_noMV_temp3)
    }

    # Combine all groups back into one matrix
    DataMatrix_NoMV <- cbind(as.data.frame(mget(nom_DataMatrix_gp_noMV[seq_along(nom_DataMatrix_gp)])))
    colnames(DataMatrix_NoMV) <- colnames(DataMatrix)

    # Create final variableMetaData file with NaN type per feature per group
    variableMetadata_MV <- cbind(variableMetadata, VariableData_MARMNAR)

    # Create final sampleMetaData file with NaN percentage per sample
    sampleMetadata_MV <- cbind(sampleMetadata, SM_new)

    return(list(Datamatrix = DataMatrix_NoMV, sampleMetadata = sampleMetadata_MV, VariableMetadata = variableMetadata_MV))
}
# ========================================================.


# ================= fix_head =============================.
#' Restores row names from first column of a data frame in which column names are in first column
#'
#' This function assumes that the first column of the input data frame contains row names.
#' It converts this column into row names then removes it.
#' This function is typically used when reading tabular files where row names
#' were imported as a regular column (e.g., using read.table with row.names = FALSE).
#'
#' @param df Dataframe with header = TRUE and row.names = FALSE
#'
#' @return Dataframe with restored row names and without the first column
#'
#' @author
#' Original code: Anthony SLEIMAN
#' Maintainer: Mélanie PETERA
fix_head <- function(df) {
    if (any(duplicated(df[, 1]))) {
        stop("First column contains duplicated row names.")
    }
    if (any(is.na(df[, 1]))) {
        stop("First column contains NA values, cannot set rownames.")
    }

    rownames(df) <- df[, 1]
    df2 <- df[, -1, drop = FALSE]
    df2
}
# ========================================================.


# ====================== plt_density =====================.
#' Comapre intensity distribution before and after imputation
#'
#' This function plots density of Log10-transformed intensity values from two data
#' matrices: one before imputation and one after imputation.
#'
#' Missing values are removed prior to density estimation
#'
#' @param df1 Data matrix before imputation to be plotted
#' @param df2 Data matrix after imputation to be plotted
#'
#' @return A base R plot showing density curves before (solid cyan line) and after (dashed magenta line) imputation
#'
#' @author
#' Original code: Anthony SLEIMAN
#' Maintainer: Mélanie PETERA
plt_density <- function(df1, df2) {
    # Choose plot colors
    color <- c("#009B9F", "#C75DAA")

    # Flatten matrix and remove missing values
    logDataVector1 <- log10(na.omit(unlist(df1)))

    # Restore structure if first column is row identifiers (which is the case in our W4M case)
    df2_fixed <- fix_head(df2)
    # Flatten second matrix and remove missing values
    logDataVector2 <- log10(na.omit(unlist(df2_fixed)))

    # Density of original data
    plot(density(logDataVector1),
        col = adjustcolor(color[1], alpha.f = 0.7),
        lty = 2,
        lwd = 2,
        main = "Density before vs after imputation",
        xlab = expression(log[10](Intensity))
    )

    # Density after imputation
    lines(density(logDataVector2),
        col = adjustcolor(color[2], alpha.f = 0.7), lwd = 2
    )

    legend(
        x = "topright",
        legend = c("Before", "After"),
        lty = c(2, 1),
        col = color,
        lwd = 2
    )
}
# ========================================================.


# ====================== plt_matrix ======================.
#' Visualize missing value pattern map in a data matrix
#' Useful for visualizing the structure of missingness across samples and features.
#' Detect batch effects or systematic missingness patterns.
#'
#' Converts a data matrix into a binary presence/absence map.
#' Missing values are encoded as 0 and observed intensities as 1.
#'
#' @param df Numeric data matrix
#'
#' @return Presence/absence map of data matix intensities. NA in black, intensities in white
#'
#' @author
#' Original code: Anthony SLEIMAN
#' Maintainer: Mélanie PETERA
plt_matrix <- function(df_original) {
    # Convert to binary presence/absence matrix (Na -> 0 ; everything else -> 1)
    df <- ifelse(is.na(df_original), 0, 1)
    # Transpose so samples become columns
    df <- t(df)

    # Chosen colors
    color <- c("black", "white") # black for missing values, white for non missing values

    # Plot binary matrix
    image(as.matrix(df),
        col = color,
        xaxt = "n", yaxt = "n"
    )
    axis(1, at = seq(from = 0, to = 1, length.out = length(row.names(df))), labels = FALSE)
    text(
        x = seq(from = 0, to = 1, length.out = length(row.names(df))),
        y = par("usr")[3] - 0.01,
        labels = row.names(df),
        srt = 40,
        adj = 1,
        xpd = TRUE,
        cex = 0.6
    )
    title(main = "Presence/absence map of data matrix intensities", col.main = "black")
    mtext("Samples", side = 3, line = 0)
    mtext("Features", side = 2, line = 0)
    par(xpd = FALSE)
    legend(
        x = "topright", xpd = TRUE, inset = c(0, -0.07),
        legend = c("MV", "No_MV"), bty = "n",
        fill = color,
        border = "black",
        cex = 0.6, horiz = TRUE
    )
}
# ========================================================.


# ====================== plt_group =======================.
#' Visualize missing value categories per sample group

#' Computes the proportion of missing value categories (No_MV, MA, MNAR) for
#' each biological group and displays them as stacked bar plot
#'
#' @param SM Sample metadata containing samples grouping information
#' @param VM Variable metadata containing missing value classification (output of imputemsmv)
#' @param ID Column name in SM defining sample groups
#'
#' @return A stacked bar plot showing percentage distribution of MV types per sample group
#'
#' @author
#' Original code: Anthony SLEIMAN
#' Maintainer: Mélanie PETERA
plt_group <- function(SM, VM, ID) {
    # Extract group names from metadata
    group_name <- as.character(unique(SM[, ID]))
    # Corresponding VM columns
    group_cols <- paste0("mv_type_", group_name)
    groups <- as.data.frame(VM[, group_cols])

    n_features <- nrow(groups)

    # Compute percentage per category per group
    info <- data.frame(type = c("NMV", "MAR", "MNAR"))
    for (i in seq_along(group_name)) {
        # Count categories
        No_MV <- sum(groups[, i] == "No_MV")
        MAR <- sum(groups[, i] == "MAR")
        MNAR <- sum(groups[, i] == "MNAR")
        # convert to percentage
        col_vals <- c(No_MV, MAR, MNAR) / n_features * 100
        # Append as a new column
        info <- cbind(info, col_vals)
    }

    # Fix column structure (using helper function fix_head from above)
    info <- fix_head(info)
    colnames(info) <- group_name
    info <- as.matrix(info)

    # Choose colors per category
    color <- c("white", "#6551CC", "#B22C2C") # white for No_MV, purple for MAR and maroon for MNAR

    # Stacked bar plot
    bp <- barplot(info, col = color, yaxt = "n", log = "")
    # Add percentage labels (centered in stacks)
    # Compute label positions
    y_pos <- apply(info, 2, cumsum) - info / 2
    # Flatten everything
    x_vals <- rep(bp, each = nrow(info))
    y_vals <- as.vector(y_pos)
    labels <- paste0(round(as.vector(info), 1), "%")
    # Keep only non-zero values
    keep <- as.vector(info) > 0
    text(
        x = x_vals[keep],
        y = y_vals[keep],
        labels = labels[keep],
        cex = 0.6,
        font = 2
    )
    title(main = "Missing values attribution type per sample group")
    par(xpd = FALSE)
    legend(
        x = "topright", xpd = TRUE, inset = c(0, -0.13),
        legend = c("No_MV", "MAR", "MNAR"),
        fill = color,
        border = "black",
        cex = 0.6, horiz = FALSE
    )
}
# ========================================================.


# ====================== plt_type ========================.
#' Visualize distribution of groups within each missing value type.
#'
#' Computes how each missing value category (No_MN, MAR, MNAR) is distributed across sample groups.
#' For each type, proportions are normalized across groups so that each category sums to 100%.
#' Allows identification of whether certain groups are enriched in specific missingness mechanisms.
#'
#' @param SM Sample metadata containing samples grouping information
#' @param VM Variable metadata containing missing value classification (output of imputemsmv)
#' @param ID Column name in SM defining sample groups
#'
#' @return A stacked barplot showing the proportion of groups within each category
#'
#' @author
#' Original code: Anthony SLEIMAN
#' Maintainer: Mélanie PETERA
plt_type <- function(SM, VM, ID) {
    # Extract group names from metadata
    group_name <- as.character(unique(SM[, ID]))
    # Corresponding VM columns
    group_cols <- paste0("mv_type_", group_name)
    groups <- as.data.frame(VM[, group_cols])

    # Count occurences per group for each MV category
    No_MV <- colSums(groups == "No_MV")
    MAR <- colSums(groups == "MAR")
    MNAR <- colSums(groups == "MNAR")

    # Normalize occurences per group for each MV category
    No_MV <- (No_MV / sum(No_MV)) * 100
    MAR <- (MAR / sum(MAR)) * 100
    MNAR <- (MNAR / sum(MNAR)) * 100

    # Combine into matrix (rows are groups, columns are MV category)
    info <- cbind(No_MV, MAR, MNAR)

    # Choose colors
    color <- hcl.colors(nrow(info), palette = "Tropic")

    # Stacked bar plot
    bp <- barplot(info, col = color, yaxt = "n", log = "")
    text(
        x = rep(bp, each = nrow(info)),
        y = apply(info, 2, cumsum) - info / 2,
        labels = paste0(round(info, 2), "%"), cex = 0.6, font = 2
    )
    title(main = "Distribution of sample groups across NA type", col.main = "black")
    par(xpd = FALSE)
    legend(
        x = "topright", xpd = TRUE, inset = c(-0.025, 0),
        legend = group_name, title = "Groups",
        fill = hcl.colors(length(groups), palette = "Tropic"),
        border = "black",
        cex = 0.6, horiz = FALSE
    )
}
# ========================================================.


# ====================== plt_pie =========================.
#' Plot global distribution of values' category (No_MV, MAR, MNAR)
#'
#' Compute the overall proportion of missing value categories across the entire dataset and visualize it as a pie chart
#'
#' @param SM Sample metadata containing samples grouping information
#' @param VM Variable metadata containing missing value classification (output of imputemsmv)
#' @param ID Column name in SM defining sample groups
#'
#' @return A pie chart showing global propotions of No_MV, MAR, MNAR
#'
#' @author
#' Original code: Anthony SLEIMAN
#' Maintainer: Mélanie PETERA
plt_pie <- function(SM, VM, ID) {
    # Extract group names from metadata
    group_name <- as.character(unique(SM[, ID]))
    # Corresponding VM columns
    group_cols <- paste0("mv_type_", group_name)
    groups <- as.data.frame(VM[, group_cols])
    # Total number of entries (features x groups)
    numero <- nrow(groups) * ncol(groups)

    # Count occurences of each missing value type
    No_MV <- sum(groups == "No_MV")
    MAR <- sum(groups == "MAR")
    MNAR <- sum(groups == "MNAR")

    # Convert counts to percentages
    pNo_MV <- (No_MV / numero) * 100
    pMAR <- (MAR / numero) * 100
    pMNAR <- (MNAR / numero) * 100

    # Combine into vector
    info <- c(pNo_MV, pMAR, pMNAR)

    # Choose color
    color <- c("white", "#6551CC", "#B22C2C") # white for No_MV, purple for MAR and maroon for MNAR
    label <- c("No_MV", "MAR", "MNAR")

    # Check if a category is absent
    info <- info[info > 0]
    color <- color[info > 0]
    label <- label[info > 0]

    # Plot pie chart
    pie(info,
        labels = paste0(round(info, 2), "%"),
        col = color
    )
    title(main = "Pie chart of values category in datamatrix", col.main = "black")
    par(xpd = FALSE)
    legend(
        x = "topright",
        legend = label,
        fill = color,
        border = "black",
        cex = 0.6, horiz = FALSE
    )
}
# ========================================================.


# ====================== plt_distribution ================.
#' Plot total number of missing values per intensity bin
#'
#' This function explores the relationship between feature intensity and missingness
#' For each feature (row), the minimum observed intensity is computed, then features are binned on this value (log10 scale)
#' The total number within each bin is then summed and displayed as a bar plot
#' Help assess MNAR behavior when the minimum-method is chosen
#'
#' @param df Numeric data matrix
#' @param bins Number of bins used to cut intesity values
#'
#' @return A bar plot showing total NA per intensity bin
#'
#' @author
#' Original code: Anthony SLEIMAN
#' Maintainer: Mélanie PETERA
plt_distribution <- function(df_original, bins) {
    # Check total missingness and remove empty features
    df <- df_original[rowSums(is.na(df_original)) != ncol(df_original), ]
    # Compute per-feature minimum intensity
    row_min <- apply(df, 1, function(x) min(x, na.rm = TRUE))
    # Log-transform
    row_min <- log10(row_min + 1e-9)
    # Count number of NA per feature
    row_na <- apply(df, 1, function(x) sum(is.na(x)))
    # Combine into matrix
    added_min <- cbind(row_min, row_na)

    # Bin features based on their minimum intensity
    bins <- cut(added_min[, 1], breaks = bins, include.lowest = TRUE)
    # Sum total NA per bin
    na_sum <- tapply(added_min[, 2], bins, sum)

    # Plot bar plot
    par(mgp = c(2, 0.5, 0))
    barplot(na_sum,
        las = 2,
        col = "#009B9F",
        border = NA,
        angle = 45,
        main = "Repartition of missing values according to the minimum of their corresponding ions",
        ylab = "Missing values count",
        cex.names = 0.7, cex.lab = 0.7, cex.axis = 0.7
    )
    mtext(expression(log[10](minimum_intensity_values)),
        side = 1,
        line = 3, cex = 0.7
    )
}
# ========================================================.
