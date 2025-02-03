if (FALSE) {
    rm(list = ls())

    DM1 <- data.frame(
        data = c("5d_-kkcùf", "npèt", "5PY4(*3"),
        `j 785` = c(0.356426723610756, 0.801750949101246, 0.875199970420953),
        `y54j 68y4j6` = c(0.380152310071702, 0.535593104115636, 0.0825428101366147),
        `5-6 4` = c(0.0306944207412024, 0.258351312473067, 0.253659010703906),
        hrrrrh = c(0.334137638848017, 0.599475573145688, 0.507762246807195),
        `5h -` = c(0.298147485608469, 0.0763319665667417, 0.856444177031262)
    )

    DM2 <- data.frame(
        data = c("5d_-kkcùf", "npèt", "5PY4(*3"),
        `j 785` = c(0.356426723610756, 0.801750949101246, 0.875199970420953),
        `y54j 68y4j6` = c(0.380152310071702, 0.535593104115636, 0.0825428101366147),
        `5-6 4` = c(0.0306944207412024, 0.258351312473067, 0.253659010703906),
        hrrrrh = c(0.334137638848017, 0.599475573145688, 0.507762246807195),
        `5h -` = c(0.298147485608469, 0.0763319665667417, 0.856444177031262)
    )

    M1 <- data.frame(
        samplename = c("j 785", "y54j 68y4j6", "5-6 4", "hrrrrh", "5h -"),
        ABD = c(19, 24, 2, 3, "y"), E = c(9, "p0", 45, 24, 29),
        AAA = c("r", "bg", "il", "d", "b"),
        fp = c("pj", "z", "e", "r", "t"),
        uv = c("s", "d", "f", "s", "d")
    )

    M2 <- data.frame(
        samplename = c("j 785", "y54j 68y4j6", "5-6 4", "hrrrrh", "5h -"),
        ABD = c(19, 24, 2, 3, "y"), E = c(9, "ici", 45, 24, 29),
        AAA = c("r", "bg", "il", "d", "b"),
        fp = c("pj", "z", "e", "r", "t"),
        uv = c("s", "d", "f", "s", "d")
    )
    type <- "sample"
    concatenation <- "unique"
    tab1 <- "tab1"
    tab2 <- "tab2"
    choice_keep <- "oui"
    keep <- 0
    concat(DM1, M1, DM2, M2, type, tab1, tab2, concatenation, choice_keep, keep)
}

#################################################################################################################


concat <- function(DM1, M1, DM2, M2, type, tab1, tab2, concatenation, choice_keep, keep) {
    # DM1/DM2 = data.frame containing data Matrix
    # M1/M2 = data.frame containing sample Metadata or variable Metadata
    # type = "sample" or "variable" depending on Metadata content
    # tab1/tab2 = Suffix for Metadata 1/2
    # concatenation = type of concatenation
    # choice_keep = choice of keeping columns with the same no or keeping just one
    # keep = keep the column in M1 or M2
    # returns the  concatenated metadata and the two Data Matrix

    identifiers_1 <- colnames(M1)[1]
    identifiers_2 <- colnames(M2)[1]

    err.stock <- NULL

    # Concatenation------------------------------------------------------------------

    # If Metadatas is Sample_Metadata we transpose
    if (type == "sample") {
        rownames(DM1) <- DM1[, 1]
        corner_DM1 <- colnames(DM1)[1]
        DM1 <- DM1[, -1, drop = FALSE]
        DM1 <- t(DM1)
        DM1 <- data.frame(sample = row.names(DM1), DM1, check.names = FALSE)
        rownames(DM1) <- NULL

        rownames(DM2) <- DM2[, 1]
        corner_DM2 <- colnames(DM2)[1]
        DM2 <- DM2[, -1, drop = FALSE]
        DM2 <- t(DM2)
        DM2 <- data.frame(sample = row.names(DM2), DM2, check.names = FALSE)
        rownames(DM2) <- NULL
    }

    # Add order of sample and Sort by order


    M1$order1 <- seq(1, nrow(M1))
    M2$order2 <- seq(nrow(M1) + 1, nrow(M2) + nrow(M1))

    M1_bf <- M1[order(M1[, 1]), ]
    M2_bf <- M2[order(M2[, 1]), ]


    # Check the variables in common and extract them.


    same <- check_features(M1_bf, M2_bf)
    same <- same[-which(same == identifiers_1)]

    # Check that shared variables have the same values.
    # If not, they are renamed or deleted according to the parameters chosen by the user.
    result2 <- compare_same_columns(M1_bf, M2_bf, same, choice_keep, keep, tab1, tab2)
    M1 <- result2$M1
    M2 <- result2$M2

    # Unique--------------------------------------------------------------------------
    if (concatenation == "unique") {
        # Table match check
        # We verify that the individuals are all the same
        err.stock <- match2_bis(M1, M2, type)
        check_err(err.stock)
        M_merge <- merge(M1, M2, by = 1)
    }


    # Intersection--------------------------------------------------------------------

    if (concatenation == "intersection") {
        # select individuals in common
        sample_common <- intersect(M1[, 1], M2[, 1])

        # if the list of individuals in common is null, an error message is sent
        if (length(sample_common) == 0) {
            err.stock <- c(err.stock, "\nThere are no individuals in common \n")
            check_err(err.stock)
        }
        # if the list of individuals in common is less than 5, then a Warning message is sent
        if (length(sample_common) < 5) {
            cat("\nWarning: Less than 5 individuals in common\n")
        }
        M_merge <- merge(M1, M2, by = 1)
    }

    # Union --------------------------------------------------------------------------
    if (concatenation == "union") {
        # select common ids
        id_common <- intersect(M1[, 1], M2[, 1])

        if (is.null(id_common)) {
            cat("\nT Warning : there are no individuals in common\n")
        }

        M2_common <- M2[M2[, 1] %in% id_common, ]
        # Store rows with individuals belonging only to M2
        M2_specifique <- M2[!M2[, 1] %in% id_common, ]
        # Merge the two tables only with the samples not in common
        M_merge <- bind_rows(M1, M2_specifique)
        col_names <- colnames(M2_common)
        col_names <- col_names[-which(col_names == identifiers_2)]
        feature_common <- check_features(M_merge, M2_bf)
        # Check if M_merge and M2_bf have columns in common. If so, complete the table with the values not taken.
        if (!is.null(feature_common)) {
            identifiers_3 <- M2_specifique[, 1]
            # We select the value in M2_bf, the M2 table before undergoing any changes, then insert it in the M_merge table.
            for (feature in feature_common) {
                for (id in identifiers_3) {
                    index_row <- which(M2_bf[, 1] == id)
                    index_col <- which(colnames(M2_bf) == feature)
                    new_value <- M2_bf[index_row, index_col]
                    index_row <- which(M_merge[, 1] == id)
                    index_col <- which(colnames(M_merge) == feature)
                    M_merge[index_row, index_col] <- new_value
                }
            }
        }
        # Fill in the table with common values
        for (col in col_names) {
            for (id in id_common) {
                index_row <- which(M2_common[, 1] == id)
                index_col <- which(colnames(M2_common) == col)
                new_value <- M2_common[index_row, index_col]
                index_row <- which(M_merge[, 1] == id)
                index_col <- which(colnames(M_merge) == col)
                M_merge[index_row, index_col] <- new_value
            }
        }
    }
    M_merge_sort <- M_merge[order(M_merge$order1, M_merge$order2), ]
    M_merge_sort <- M_merge_sort[, -which(colnames(M_merge_sort) == "order1")]
    M_merge_sort <- M_merge_sort[, -which(colnames(M_merge_sort) == "order2")]
    # DataMatrix ---------------------------------------------------------------------

    colnames_1 <- colnames(DM1)
    colnames_2 <- colnames(DM2)
    # Unique -------------------------------------------------------------------------

    if (concatenation == "unique") {
        if (type == "sample") {
            rownames(DM1) <- DM1[, 1]
            DM1 <- DM1[, -1]
            DM1 <- t(DM1)
            DM1 <- data.frame(sample = row.names(DM1), DM1, check.names = FALSE)
            colnames(DM1)[1] <- corner_DM1
            rownames(DM1) <- NULL

            rownames(DM2) <- DM2[, 1]
            DM2 <- DM2[, -1, drop = FALSE]
            DM2 <- t(DM2)
            DM2 <- data.frame(sample = row.names(DM2), DM2, check.names = FALSE)
            colnames(DM2)[1] <- corner_DM2
            rownames(DM2) <- NULL
        }
        result <- list(M_merge_sort = M_merge_sort, DM1 = DM1, DM2 = DM2)
        return(result)
    }

    # Intersection--------------------------------------------------------------------

    if (concatenation == "intersection") {
        id_in_common <- intersect(DM1[, 1], DM2[, 1])

        DM1_filter <- subset(DM1, DM1[, 1] %in% id_in_common)
        DM2_filter <- subset(DM2, DM2[, 1] %in% id_in_common)

        if (type == "sample") {
            rownames(DM1_filter) <- DM1_filter[, 1]
            DM1_filter <- DM1_filter[, -1]
            DM1_filter <- t(DM1_filter)
            DM1_filter <- data.frame(sample = row.names(DM1_filter), DM1_filter, check.names = FALSE)
            colnames(DM1_filter)[1] <- corner_DM1
            rownames(DM1_filter) <- NULL

            rownames(DM2_filter) <- DM2_filter[, 1]
            DM2_filter <- DM2_filter[, -1, drop = FALSE]
            DM2_filter <- t(DM2_filter)
            DM2_filter <- data.frame(sample = row.names(DM2_filter), DM2_filter, check.names = FALSE)
            colnames(DM2_filter)[1] <- corner_DM2
            rownames(DM2_filter) <- NULL
        }
        result <- list(M_merge_sort = M_merge_sort, DM1 = DM1_filter, DM2 = DM2_filter)
        return(result)
    }

    # Union --------------------------------------------------------------------------

    if (concatenation == "union") {
        common_individuals <- intersect(DM1[, 1], DM2[, 1])
        common_columns <- intersect(colnames_1, colnames_2)
        # check whether there are individuals or variables in common
        if (is.null(common_individuals) || is.null(common_columns)) {
            comparison_result <- FALSE
            # If the individuals in common take the same values for all variables, then comparison_result=TRUE
        } else {
            DM1_common <- subset(DM1, DM1[, 1] %in% common_individuals)
            DM2_common <- subset(DM2, DM2[, 1] %in% common_individuals)
            DM1_common <- DM1_common[, common_columns, drop = FALSE]
            DM2_common <- DM2_common[, common_columns, drop = FALSE]

            for (col in common_columns) {
                comparison_result <- identical(DM1_common$col, DM2_common$col)
            }
        }

        if (comparison_result) {
            DM1$order1 <- seq(1, nrow(DM1))
            DM2$order2 <- seq(nrow(DM1) + 1, nrow(DM2) + nrow(DM1))
            DM1_sort <- DM1[order(DM1[, 1]), ]
            DM2_sort <- DM2[order(DM2[, 1]), ]
            id_in_common <- intersect(DM1[, 1], DM2[, 1])
            DM1_filter <- subset(DM1, DM1[, 1] %in% id_in_common)
            DM2_filter <- subset(DM2, DM2[, 1] %in% id_in_common)
            different_DM2 <- colnames_2[!colnames_2 %in% colnames_1]
            DM2_specifique <- DM2[!DM2[, 1] %in% id_in_common, ]
            # Merge the two tables only with the samples not in common
            DM1_merge <- bind_rows(DM1, DM2_specifique)


            # Deletion of columns present only in DM2
            DM1_merge <- DM1_merge[, !names(DM1_merge) %in% different_DM2]
            different_DM1 <- colnames_1[!colnames_1 %in% colnames_2]
            DM1_specifique <- DM1[!DM1[, 1] %in% id_in_common, ]
            # Merge the two tables only with the samples not in common
            DM2_merge <- bind_rows(DM2, DM1_specifique)
            # Deletion of columns present only in DM2
            DM2_merge <- DM2_merge[, !names(DM2_merge) %in% different_DM1]
            # DM2_merge


            DM1_merge_sort <- DM1_merge[order(DM1_merge$order1, DM1_merge$order2), ]
            DM1_merge_sort <- DM1_merge_sort[, -which(colnames(DM1_merge_sort) == "order1")]
            DM1_merge_sort <- DM1_merge_sort[, -which(colnames(DM1_merge_sort) == "order2")]

            DM2_merge_sort <- DM2_merge[order(DM2_merge$order1, DM2_merge$order2), ]
            DM2_merge_sort <- DM2_merge_sort[, -which(colnames(DM2_merge_sort) == "order1")]
            DM2_merge_sort <- DM2_merge_sort[, -which(colnames(DM2_merge_sort) == "order2")]


            if (type == "sample") {
                rownames(DM1_merge_sort) <- DM1_merge_sort[, 1]
                DM1_merge_sort <- DM1_merge_sort[, -1]
                DM1_merge_sort <- t(DM1_merge_sort)
                DM1_merge_sort <- data.frame(sample = row.names(DM1_merge_sort), DM1_merge_sort, check.names = FALSE)
                colnames(DM1_merge_sort)[1] <- corner_DM1
                rownames(DM1_merge_sort) <- NULL

                rownames(DM2_merge_sort) <- DM2_merge_sort[, 1]
                DM2_merge_sort <- DM2_merge_sort[, -1, drop = FALSE]
                DM2_merge_sort <- t(DM2_merge_sort)
                DM2_merge_sort <- data.frame(sample = row.names(DM2_merge_sort), DM2_merge_sort, check.names = FALSE)
                colnames(DM2_merge_sort)[1] <- corner_DM2
                rownames(DM2_merge_sort) <- NULL
            }

            result <- list(M_merge_sort = M_merge_sort, DM1 = DM1_merge_sort, DM2 = DM2_merge_sort)
            return(result)
        } else {
            # selects line ids that are in DM2 and not in DM1
            id_diff_1 <- setdiff(DM2[, 1], DM1[, 1])
            # we store them in a dataframe
            row_add_1 <- data.frame(id = id_diff_1)
            # renames columns with their names in DM1
            colnames(row_add_1)[1] <- colnames(DM1)[1]
            # Merge
            DM1_add <- bind_rows(DM1, row_add_1)
            id_diff_2 <- setdiff(DM1[, 1], DM2[, 1])
            row_add_2 <- data.frame(id = id_diff_2)
            colnames(row_add_2)[1] <- colnames(DM2)[1]
            DM2_add <- bind_rows(DM2, row_add_2)

            if (type == "sample") {
                rownames(DM1_add) <- DM1_add[, 1]
                DM1_add <- DM1_add[, -1]
                DM1_add <- t(DM1_add)
                DM1_add <- data.frame(sample = row.names(DM1_add), DM1_add, check.names = FALSE)
                colnames(DM1_add)[1] <- corner_DM1
                rownames(DM1_add) <- NULL
                rownames(DM2_add) <- DM2_add[, 1]
                DM2_add <- DM2_add[, -1, drop = FALSE]
                DM2_add <- t(DM2_add)
                DM2_add <- data.frame(sample = row.names(DM2_add), DM2_add, check.names = FALSE)
                colnames(DM2_add)[1] <- corner_DM2
                rownames(DM2_add) <- NULL
            }
            result <- list(M_merge_sort = M_merge_sort, DM1 = DM1_add, DM2 = DM2_add)
            return(result)
        }
    }
}
