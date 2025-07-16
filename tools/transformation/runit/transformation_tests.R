test_input_log10 <- function() {

    testDirC <- "input"
    argLs <- list(method = "log10")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[['datMN']]['HU_021', 'HMDB03072'], 4.817089, tolerance = 1e-6)
 
}

test_input_sqrt <- function() {

    testDirC <- "input"
    argLs <- list(method = "sqrt")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[['datMN']]['HU_021', 'HMDB03072'], 256.177640, tolerance = 1e-6)
 
}

