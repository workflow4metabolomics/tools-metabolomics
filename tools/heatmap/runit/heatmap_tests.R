test_input_cut4 <- function() {
    testDirC <- "input"
    argLs <- list(
        dis_c = "1-cor",
        cut_sam_n = "4",
        cut_var_n = "3",
        cor_met_c = "spearman",
        agg_met_c = "ward",
        col_c = "blueOrangeRed",
        sca_l = "TRUE",
        cex_n = "0.8"
    )

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["samDF"]][13, "heat_clust"], 4)
}

test_exa1_cut3 <- function() {
    testDirC <- "exa1"
    argLs <- list(
        dis_c = "1-cor",
        cut_sam_n = "3",
        cut_var_n = "4",
        cor_met_c = "spearman",
        agg_met_c = "ward",
        col_c = "blueOrangeRed",
        sca_l = "TRUE",
        cex_n = "1"
    )

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["varDF"]]["V24", "heat_clust"], 4)
}

test_exa2_cut4 <- function() {
    testDirC <- "exa2"
    argLs <- list(
        dis_c = "1-cor",
        cut_sam_n = "1",
        cut_var_n = "1",
        cor_met_c = "spearman",
        agg_met_c = "ward",
        col_c = "blueOrangeRed",
        sca_l = "TRUE",
        cex_n = "1"
    )

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEquals(outLs[["varDF"]]["V31", "meta2"], "AM")
}
