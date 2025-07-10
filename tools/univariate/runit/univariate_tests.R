test_input_anova <- function() {

    testDirC <- "input"
    argLs <- list(facC = "qual",
                  tesC = "anova",
                  adjC = "fdr",
                  thrN = "0.05")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["varDF"]]["v6", "qual_anova_fdr"], 1.924156e-03, tolerance = 1e-6)

    checkEqualsNumeric(outLs[["varDF"]]["v4", "qual_anova_D.C_fdr"], 0.01102016, tolerance = 1e-6)

}

test_input_kruskal <- function() {

    testDirC <- "input"
    argLs <- list(facC = "qual",
                  tesC = "kruskal",
                  adjC = "fdr",
                  thrN = "0.05")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["varDF"]]["v4", "qual_kruskal_fdr"], 0.0008194662, tolerance = 1e-7)

    checkEqualsNumeric(outLs[["varDF"]]["v6", "qual_kruskal_D.A_fdr"], 0.002945952, tolerance = 1e-7)

}

test_example1_wilcoxDif <- function() {

    testDirC <- "example1"
    argLs <- list(facC = "jour",
                  tesC = "wilcoxon",
                  adjC = "fdr",
                  thrN = "0.05")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
    
    checkEqualsNumeric(outLs[["varDF"]]["MT3", "jour_wilcoxon_J3.J10_dif"], 0.216480042, tolerance = 1e-8)

}

test_example1_ttestFdr <- function() {

    testDirC <- "example1"
    argLs <- list(facC = "jour",
                  tesC = "ttest",
                  adjC = "fdr",
                  thrN = "0.05")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["varDF"]]["MT3", "jour_ttest_J3.J10_fdr"], 0.7605966, tolerance = 1e-6)

}
