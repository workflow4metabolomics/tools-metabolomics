test_input_pca <- function() {

    testDirC <- "input"
    argLs <- list(respC = "none",
                  predI = "NA",
                  orthoI = "0",
                  testL = "FALSE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_input_pcaGender <- function() {

    testDirC <- "input"
    argLs <- list(respC = "none",
                  predI = "NA",
                  orthoI = "0",
                  testL = "FALSE",
                  parMahalC = "gender")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_input_plsdaGender <- function() {

    testDirC <- "input"
    argLs <- list(respC = "gender",
                  predI = "NA",
                  orthoI = "0",
                  testL = "FALSE")
    
    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_input_oplsAge <- function() {

    testDirC <- "input"
    argLs <- list(respC = "age",
                  predI = "1",
                  orthoI = "1",
                  testL = "FALSE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_input_oplsdaGender <- function() {

    testDirC <- "input"
    argLs <- list(respC = "gender",
                  predI = "1",
                  orthoI = "1",
                  testL = "FALSE")
    
    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_sacurine_pca <- function() {

    testDirC <- "sacurine"
    argLs <- list(respC = "none",
                  predI = "NA",
                  orthoI = "0",
                  testL = "FALSE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_sacurine_pcaGender <- function() {

    testDirC <- "sacurine"
    argLs <- list(respC = "none",
                  predI = "NA",
                  orthoI = "0",
                  testL = "FALSE",
                  parMahalC = "gender")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_sacurine_plsAge <- function() {

    testDirC <- "sacurine"
    argLs <- list(respC = "age",
                  predI = "NA",
                  orthoI = "0",
                  testL = "FALSE")
    
    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_sacurine_plsdaGender <- function() {

    testDirC <- "sacurine"
    argLs <- list(respC = "gender",
                  predI = "NA",
                  orthoI = "0",
                  testL = "FALSE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_sacurineTest_pls <- function() {

    testDirC <- "sacurineTest"
    argLs <- list(respC = "age",
                  predI = "2",
                  orthoI = "0",
                  testL = "TRUE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["samDF"]][181, "age_PLS_predictions"], 40.82252, tolerance = 1e-5)
 
}

test_sacurineTest_opls <- function() {

    testDirC <- "sacurineTest"
    argLs <- list(respC = "age",
                  predI = "1",
                  orthoI = "2",
                  testL = "TRUE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["samDF"]][181, "age_OPLS_predictions"], 40.28963, tolerance = 1e-5)
 
}

test_sacurineTest_plsda <- function() {

    testDirC <- "sacurineTest"
    argLs <- list(respC = "gender",
                  predI = "2",
                  orthoI = "0",
                  testL = "TRUE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEquals(outLs[["samDF"]][181, "gender_PLSDA_predictions"], "F")
 
}

test_sacurineTest_oplsda <- function() {

    testDirC <- "sacurineTest"
    argLs <- list(respC = "gender",
                  predI = "1",
                  orthoI = "1",
                  testL = "TRUE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEquals(outLs[["samDF"]][181, "gender_OPLSDA_predictions"], "F")
 
}

test_sacurine_oplsAge <- function() {

    testDirC <- "sacurine"
    argLs <- list(respC = "age",
                  predI = "1",
                  orthoI = "1",
                  testL = "FALSE")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[["varDF"]][1, "age_OPLS_VIP_ortho"], 0.3514378, tolerance = 1e-7)
}

test_example1_plsda <- function() {

    testDirC <- "example1"
    argLs <- list(respC = "traitment",
                  predI = "3",
                  orthoI = "0",
                  testL = "FALSE")
    
    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}

test_example2_pca <- function() {

    testDirC <- "example2"
    argLs <- list(respC = "none",
                  predI = "NA",
                  orthoI = "0",
                  testL = "FALSE")
    
    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
}
