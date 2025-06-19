suppressPackageStartupMessages(library(argparse))

check.err <- function(err.stock) {
    # err.stock = vector of results returned by check functions
    if (length(err.stock) != 0) {
        stop("\n- - - - - - - - -\n", err.stock, "\n- - - - - - - - -\n")
    }
}


# To check if the 3 standard tables match regarding identifiers
table_match <- function(dataMatrix, sampleMetadata, variableMetadata) {
    err.stock <- character()

    # Sample identifiers
    dm_samples <- colnames(dataMatrix)[-1]
    sm_samples <- as.character(sampleMetadata[[1]])
    # Variable identifiers
    dm_vars <- as.character(dataMatrix[[1]])
    vm_vars <- as.character(variableMetadata[[1]])

    # Check sample IDs
    missing_in_sm <- setdiff(dm_samples, sm_samples)
    missing_in_dm <- setdiff(sm_samples, dm_samples)
    if (length(missing_in_sm) > 0 || length(missing_in_dm) > 0) {
        err.stock <- c(
            err.stock,
            "\nData matrix and sample metadata do not match regarding sample identifiers."
        )
        if (length(missing_in_sm) > 0) {
            prefix <- if (length(missing_in_sm) < 4) {
                "\n    The "
            } else {
                "\n    For example, the "
            }
            err.stock <- c(
                err.stock,
                prefix,
                "following identifiers found in the data matrix\n",
                "    do not appear in the sample metadata file:\n",
                "    ",
                paste(head(missing_in_sm, 3), collapse = "\n    "),
                "\n"
            )
        }
        if (length(missing_in_dm) > 0) {
            prefix <- if (length(missing_in_dm) < 4) {
                "\n    The "
            } else {
                "\n    For example, the "
            }
            err.stock <- c(
                err.stock,
                prefix,
                "following identifiers found in the sample metadata file\n",
                "    do not appear in the data matrix:\n",
                "    ",
                paste(head(missing_in_dm, 3), collapse = "\n    "),
                "\n"
            )
        }
    }

    # Check variable IDs
    missing_in_vm <- setdiff(dm_vars, vm_vars)
    missing_in_dm_var <- setdiff(vm_vars, dm_vars)
    if (length(missing_in_vm) > 0 || length(missing_in_dm_var) > 0) {
        err.stock <- c(
            err.stock,
            "\nData matrix and variable metadata do not match regarding variable identifiers."
        )
        if (length(missing_in_vm) > 0) {
            prefix <- if (length(missing_in_vm) < 4) {
                "\n    The "
            } else {
                "\n    For example, the "
            }
            err.stock <- c(
                err.stock,
                prefix,
                "following identifiers found in the data matrix\n",
                "    do not appear in the variable metadata file:\n",
                "    ",
                paste(head(missing_in_vm, 3), collapse = "\n    "),
                "\n"
            )
        }
        if (length(missing_in_dm_var) > 0) {
            prefix <- if (length(missing_in_dm_var) < 4) {
                "\n    The "
            } else {
                "\n    For example, the "
            }
            err.stock <- c(
                err.stock,
                prefix,
                "following identifiers found in the variable metadata file\n",
                "    do not appear in the data matrix:\n",
                "    ",
                paste(head(missing_in_dm_var, 3), collapse = "\n    "),
                "\n"
            )
        }
    }

    if (length(err.stock) > 0) {
        err.stock <- c(err.stock, "\nPlease check your data.\n")
    }

    return(err.stock)
}


intens_check <- function(
    DM.name,
    SM.name,
    VM.name,
    method,
    chosen.stat,
    class.col,
    test.fold,
    class1,
    fold.frac,
    logarithm,
    VM.output,
    graphs.output) {
    # Read input tables
    DM <- read.table(DM.name, header = TRUE, sep = "\t", check.names = FALSE)
    SM <- read.table(SM.name, header = TRUE, sep = "\t", check.names = FALSE)
    VM <- read.table(VM.name, header = TRUE, sep = "\t", check.names = FALSE)

    # Table match check
    table.check <- table_match(DM, SM, VM)
    check.err(table.check)

    # Transpose and align DM
    rownames(DM) <- DM[, 1]
    DM <- data.frame(t(DM[, -1]))
    DM <- merge(x = cbind(1:nrow(SM), SM), y = DM, by.x = 2, by.y = 0)
    DM <- DM[order(DM[, 2]), ]
    rownames(DM) <- DM[, 1]
    DM <- DM[, -c(1:(ncol(SM) + 1))]

    stat.list <- strsplit(chosen.stat, ",")[[1]]

    # Class assignment
    if (method == "no_class") {
        c_class <- rep("global", nrow(DM))
        classnames <- "global"
        nb_class <- 1
        test.fold <- "No"
    } else {
        class.col <- colnames(SM)[as.numeric(class.col)]
        if (!(class.col %in% colnames(SM))) {
            stop("The column ", class.col, " is not in sample metadata")
        }
        c_class <- as.factor(SM[, class.col])
        nb_class <- nlevels(c_class)
        classnames <- levels(c_class)
        if (nb_class < 2 && test.fold == "Yes") {
            cat(
                "The column",
                class.col,
                "contains only one class, fold calculation could not be executed.\n"
            )
        }
        if (nb_class > (nrow(SM)) / 3 && method == "each_class") {
            cat("There are too many classes, consider reducing the number.\n")
        }
        if (method == "one_class") {
            if (!(class1 %in% classnames)) {
                stop("The class ", class1, " does not appear in the column ", class.col)
            }
            c_class <- factor(
                ifelse(c_class == class1, class1, "Other"),
                levels = c(class1, "Other")
            )
            nb_class <- nlevels(c_class)
            classnames <- levels(c_class)
        }
    }

    # Check numeric
    if (!is.numeric(as.matrix(DM))) {
        findchar <- function(myval) {
            ifelse(
                is.na(myval),
                "ok",
                ifelse(is.na(as.numeric(as.character(myval))), "char", "ok")
            )
        }
        chardiag <- suppressWarnings(apply(DM, 2, vapply, findchar, "character"))
        charlist <- which(chardiag == "char")
        err.stock <- paste(
            "\n- - - - - - - - -\nYour dataMatrix contains",
            length(charlist),
            "non-numeric value(s).\n"
        )
        stop(c(
            err.stock,
            "The dataMatrix file is supposed to contain only numeric values.\n- - - - - - - - -\n"
        ))
    }

    DM <- cbind(c_class, DM)
    stat.res <- t(DM[0, -1, drop = FALSE])
    names <- NULL

    mean.res <- sd.res <- med.res <- quart.res <- dec.res <- NA.res <- pct_NA.res <- fold.res <- NULL
    mean.names <- sd.names <- med.names <- quart.names <- dec.names <- NA.names <- pct_NA.names <- fold.names <- NULL
    graphs <- if (("NA" %in% stat.list) || (test.fold == "Yes")) 1 else 0
    data_bp <- data.frame()

    for (j in 1:nb_class) {
        idx <- which(DM$c_class == classnames[j])
        # Mean
        if ("mean" %in% stat.list) {
            mean.res <- cbind(mean.res, colMeans(DM[idx, -1], na.rm = TRUE))
            mean.names <- cbind(mean.names, paste("Mean", classnames[j], sep = "_"))
            if (j == nb_class) {
                stat.res <- cbind(stat.res, mean.res)
                names <- cbind(names, mean.names)
            }
        }
        # SD
        if ("sd" %in% stat.list) {
            sd.res <- cbind(sd.res, apply(DM[idx, -1], 2, sd, na.rm = TRUE))
            sd.names <- cbind(sd.names, paste("Sd", classnames[j], sep = "_"))
            if (j == nb_class) {
                stat.res <- cbind(stat.res, sd.res)
                names <- cbind(names, sd.names)
            }
        }
        # Median
        if (("median" %in% stat.list) && (!("quartile" %in% stat.list))) {
            med.res <- cbind(med.res, apply(DM[idx, -1], 2, median, na.rm = TRUE))
            med.names <- cbind(med.names, paste("Median", classnames[j], sep = "_"))
            if (j == nb_class) {
                stat.res <- cbind(stat.res, med.res)
                names <- cbind(names, med.names)
            }
        }
        # Quartiles
        if ("quartile" %in% stat.list) {
            quart.res <- cbind(
                quart.res,
                t(apply(DM[idx, -1], 2, quantile, na.rm = TRUE))
            )
            quart.names <- cbind(
                quart.names,
                paste("Min", classnames[j], sep = "_"),
                paste("Q1", classnames[j], sep = "_"),
                paste("Median", classnames[j], sep = "_"),
                paste("Q3", classnames[j], sep = "_"),
                paste("Max", classnames[j], sep = "_")
            )
            if (j == nb_class) {
                stat.res <- cbind(stat.res, quart.res)
                names <- cbind(names, quart.names)
            }
        }
        # Deciles
        if ("decile" %in% stat.list) {
            dec.res <- cbind(
                dec.res,
                t(apply(DM[idx, -1], 2, quantile, na.rm = TRUE, seq(0, 1, 0.1)))
            )
            dec.names <- cbind(
                dec.names,
                t(matrix(paste("D", seq(0, 10, 1), sep = "_", classnames[j])))
            )
            if (j == nb_class) {
                stat.res <- cbind(stat.res, dec.res)
                names <- cbind(names, dec.names)
            }
        }
        # Missing values
        if ("NA" %in% stat.list) {
            nb_NA <- apply(DM[idx, -1], 2, function(x) sum(is.na(x)))
            pct_NA <- round(nb_NA / nrow(DM[idx, -1]) * 100, 4)
            NA.res <- cbind(NA.res, nb_NA)
            pct_NA.res <- cbind(pct_NA.res, pct_NA)
            NA.names <- cbind(NA.names, paste("NA", classnames[j], sep = "_"))
            pct_NA.names <- cbind(
                pct_NA.names,
                paste("Pct_NA", classnames[j], sep = "_")
            )
            if (j == nb_class) {
                stat.res <- cbind(stat.res, NA.res, pct_NA.res)
                names <- cbind(names, NA.names, pct_NA.names)
            }
            # Barplot data
            bins <- cut(pct_NA, breaks = c(-1, 20, 40, 60, 80, 100), labels = FALSE)
            for (b in 1:5) data_bp[b, j] <- sum(bins == b)
            rownames(data_bp) <- c(
                "0%-20%",
                "20%-40%",
                "40%-60%",
                "60%-80%",
                "80%-100%"
            )
            if (j == nb_class) {
                if (sum(NA.res) == 0) cat("Data Matrix contains no NA.\n")
                colnames(data_bp) <- classnames
                data_bp <- as.matrix(data_bp)
            }
        }
        # Mean fold change
        if (test.fold == "Yes" && nb_class >= 2 && j != nb_class) {
            if (method == "each_class") fold.frac <- "Top"
            for (k in (j + 1):nb_class) {
                ratio1 <- if (fold.frac == "Bottom") classnames[k] else classnames[j]
                ratio2 <- if (fold.frac == "Bottom") classnames[j] else classnames[k]
                fold <- colMeans(DM[which(DM$c_class == ratio1), -1], na.rm = TRUE) /
                    colMeans(DM[which(DM$c_class == ratio2), -1], na.rm = TRUE)
                if (logarithm == "log2") {
                    fold <- log2(fold)
                } else if (
                    logarithm == "log10"
                ) {
                    fold <- log10(fold)
                }
                fold.res <- cbind(fold.res, fold)
                fname <- if (logarithm == "none") {
                    paste("fold", ratio1, "VS", ratio2, sep = "_")
                } else {
                    paste(logarithm, "fold", ratio1, "VS", ratio2, sep = "_")
                }
                fold.names <- cbind(fold.names, fname)
            }
            if (j == nb_class) {
                stat.res <- cbind(stat.res, fold.res)
                names <- cbind(names, fold.names)
            }
        }
    }

    # Ensure unique column names
    VM.names <- colnames(VM)
    for (i in seq_along(VM.names)) {
        for (j in seq_along(names)) {
            if (VM.names[i] == names[j]) names[j] <- paste(names[j], "2", sep = "_")
        }
    }
    colnames(stat.res) <- names

    # Output
    VM <- cbind(VM, stat.res)
    write.table(VM, VM.output, sep = "\t", quote = FALSE, row.names = FALSE)

    # Graphics
    pdf(graphs.output)
    if (graphs == 1) {
        if ("NA" %in% stat.list) {
            graph.colors <- c("green3", "palegreen3", "lightblue", "orangered", "red")
            par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
            bp <- barplot(
                data_bp,
                col = graph.colors,
                main = "Proportion of NA",
                xlab = "Classes",
                ylab = "Variables"
            )
            legend(
                "topright",
                fill = graph.colors,
                rownames(data_bp),
                inset = c(-0.3, 0)
            )
            stock <- 0
            for (i in 1:nrow(data_bp)) {
                text(
                    bp,
                    stock + data_bp[i, ] / 2,
                    data_bp[i, ],
                    col = "white",
                    cex = 0.7
                )
                stock <- stock + data_bp[i, ]
            }
        }
        if ((test.fold == "Yes") && (nb_class >= 2)) {
            clean_fold <- fold.res
            clean_fold[is.infinite(clean_fold)] <- NA
            for (j in 1:ncol(clean_fold)) {
                boxplot(clean_fold[, j], main = fold.names[j])
            }
        }
    } else {
        plot.new()
        legend("center", "You did not select any option with graphical output.")
    }
    dev.off()
}


parser <- ArgumentParser(description = "Intensity Check Tool")

parser$add_argument(
    "--dataMatrix_in",
    required = TRUE,
    help = "Input data matrix file"
)
parser$add_argument(
    "--sampleMetadata_in",
    required = TRUE,
    help = "Input sample metadata file"
)
parser$add_argument(
    "--variableMetadata_in",
    required = TRUE,
    help = "Input variable metadata file"
)
parser$add_argument("--method", required = TRUE, help = "Computation method")
parser$add_argument(
    "--chosen_stat",
    required = TRUE,
    help = "Statistics to compute"
)
parser$add_argument(
    "--class_col",
    default = NULL,
    help = "Class column (optional)"
)
parser$add_argument(
    "--test_fold",
    default = NULL,
    help = "Test fold (optional)"
)
parser$add_argument("--class1", default = NULL, help = "Class1 (optional)")
parser$add_argument(
    "--fold_frac",
    default = NULL,
    help = "Fold fraction (optional)"
)
parser$add_argument(
    "--logarithm",
    default = NULL,
    help = "Logarithm (optional)"
)
parser$add_argument(
    "--variableMetadata_out",
    required = TRUE,
    help = "Output variable metadata file"
)
parser$add_argument(
    "--graphs_out",
    required = TRUE,
    help = "Output graphs file"
)

args <- parser$parse_args()

print(args)

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


if (is.null(err_no_option)) {
    intens_check(
        args$dataMatrix_in,
        args$sampleMetadata_in,
        args$variableMetadata_in,
        args$method,
        args$chosen_stat,
        class_col,
        test_fold,
        class1,
        fold_frac,
        logarithm,
        args$variableMetadata_out,
        args$graphs_out
    )
}


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
