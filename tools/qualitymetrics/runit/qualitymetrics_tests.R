test_input_default <- function() {

    testDirC <- "input"
    argLs <- list(CV = "FALSE",
                  Compa = "TRUE",
                  seuil = 1,
                  poolAsPool1L = "TRUE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["varDF"]]["met_033", "blankMean_over_sampleMean"], 0.004417387, tolerance = 1e-6)
 
}

test_formatOrder <- function() {

    ## two first samples swapped in sampleMetadata

    testDirC <- "formatOrder"
    argLs <- list(CV = "FALSE",
                  Compa = "TRUE",
                  seuil = 1,
                  poolAsPool1L = "TRUE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["varDF"]]["met_033", "blankMean_over_sampleMean"], 0.004417387, tolerance = 1e-6)
 
}
