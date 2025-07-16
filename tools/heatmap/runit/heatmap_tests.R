test_input_cut4 <- function() {

    testDirC <- "input"
    argLs <- list(disC = "1-cor",
                  cutSamN = "4",
                  cutVarN = "3",
                  corMetC = "spearman",
                  aggMetC = "ward",
                  colC = "blueOrangeRed",
                  scaL = "TRUE",
                  cexN = "0.8")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)
 
    checkEqualsNumeric(outLs[['samDF']][13, 'heat_clust'], 4)

}

test_exa1_cut3 <- function() {

    testDirC <- "exa1"
    argLs <- list(disC = "1-cor",
                  cutSamN = "3",
                  cutVarN = "4",
                  corMetC = "spearman",
                  aggMetC = "ward",
                  colC = "blueOrangeRed",
                  scaL = "TRUE",
                  cexN = "1")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEqualsNumeric(outLs[['varDF']]['V24', 'heat_clust'], 4)

}

test_exa2_cut4 <- function() {

    testDirC <- "exa2"
    argLs <- list(disC = "1-cor",
                  cutSamN = "1",
                  cutVarN = "1",
                  corMetC = "spearman",
                  aggMetC = "ward",
                  colC = "blueOrangeRed",
                  scaL = "TRUE",
                  cexN = "1")

    argLs <- c(defaultArgF(testDirC), argLs)
    outLs <- wrapperCallF(argLs)

    checkEquals(outLs[['varDF']]['V31', 'meta2'], 'AM')

}

