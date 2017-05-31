test_input_kruskal <- function() {

    testDirC <- "input"
    argLs <- list(facC = "ageGroup",
                  tesC = "kruskal",
                  adjC = "fdr",
                  thrN = "0.05")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["varDF"]]["HMDB01032", "ageGroup_kruskal_senior.experienced_pva"], 0.1231246, tolerance = 1e-6)

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
