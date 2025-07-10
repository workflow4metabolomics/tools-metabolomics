test_sacurine_default <- function() {

    testDirC <- "sacurine"
    argLs <- list(respC = "gender",
                  methodC = "all",
                  bootI = "5",
                  pvalN = "0.05",
                  seedI = "123")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEquals(outLs[["varDF"]]["Oxoglutaric acid", "gender_tier_S"], "plsda|randomforest|svm")

}
