univariateF <- function(datMN,
                        samDF,
                        varDF,
                        facC,
                        tesC = c("ttest", "wilcoxon", "anova", "kruskal", "pearson", "spearman")[1],
                        adjC = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")[7],
                        thrN = 0.05,
                        pdfC) {


    ## Option

    strAsFacL <- options()$stringsAsFactors
    options(stingsAsFactors = FALSE)
    options(warn = -1)

    ## Getting the response (either a factor or a numeric)

    if(mode(samDF[, facC]) == "character") {
        facFcVn <- factor(samDF[, facC])
        facLevVc <- levels(facFcVn)
    } else
        facFcVn <- samDF[, facC]

    cat("\nPerforming '", tesC, "'\n", sep="")

    varPfxC <- paste0(make.names(facC), "_", tesC, "_")

    
    if(tesC %in% c("ttest", "wilcoxon", "pearson", "spearman")) {

        
        switch(tesC,
               ttest = {
                   staF <- function(y) diff(tapply(y, facFcVn, function(x) mean(x, na.rm = TRUE)))
                   tesF <- function(y) t.test(y ~ facFcVn)[["p.value"]]
               },
               wilcoxon = {
                   staF <- function(y) diff(tapply(y, facFcVn, function(x) median(x, na.rm = TRUE)))
                   tesF <- function(y) wilcox.test(y ~ facFcVn)[["p.value"]]
               },
               pearson = {
                   staF <- function(y) cor(facFcVn, y, method = "pearson", use = "pairwise.complete.obs")
                   tesF <- function(y) cor.test(facFcVn, y, method = "pearson", use = "pairwise.complete.obs")[["p.value"]]
               },
               spearman = {
                   staF <- function(y) cor(facFcVn, y, method = "spearman", use = "pairwise.complete.obs")
                   tesF <- function(y) cor.test(facFcVn, y, method = "spearman", use = "pairwise.complete.obs")[["p.value"]]
               })

        staVn <- apply(datMN, 2, staF)

        adjVn <- p.adjust(apply(datMN,
                                2,
                                tesF),
                          method = adjC)

        sigVn <- as.numeric(adjVn < thrN)

        if(tesC %in% c("ttest", "wilcoxon"))
            varPfxC <- paste0(varPfxC, paste(rev(facLevVc), collapse = "."), "_")

        varDF[, paste0(varPfxC, ifelse(tesC %in% c("ttest", "wilcoxon"), "dif", "cor"))] <- staVn

        varDF[, paste0(varPfxC, adjC)] <- adjVn

        varDF[, paste0(varPfxC, "sig")] <- sigVn

        ## graphic

        pdf(pdfC, onefile = TRUE)

        varVi <- which(sigVn > 0)

        if(tesC %in% c("ttest", "wilcoxon")) {

            facVc <- as.character(facFcVn)
            names(facVc) <- rownames(samDF)

            for(varI in varVi) {
                
                varC <- rownames(varDF)[varI]
                
                boxF(facFcVn,
                     datMN[, varI],
                     paste0(varC, " (", adjC, " = ", signif(adjVn[varI], 2), ")"),
                     facVc)
                
            }

        } else { ## pearson or spearman

            for(varI in varVi) {

                varC <- rownames(varDF)[varI]

                mod <- lm(datMN[, varI] ~  facFcVn)

                plot(facFcVn, datMN[, varI],
                     xlab = facC,
                     ylab = "",
                     pch = 18,
                     main = paste0(varC, " (", adjC, " = ", signif(adjVn[varI], 2), ", R2 = ", signif(summary(mod)$r.squared, 2), ")"))
            
                abline(mod, col = "red")

                }

        }

        dev.off()

        
    } else if(tesC == "anova") {

        
        ## getting the names of the pairwise comparisons 'class1Vclass2'
        prwVc <- rownames(TukeyHSD(aov(datMN[, 1] ~ facFcVn))[["facFcVn"]])

        prwVc <- gsub("-", ".", prwVc, fixed = TRUE) ## 2016-08-05: '-' character in dataframe column names seems not to be converted to "." by write.table on ubuntu R-3.3.1

        ## omnibus and post-hoc tests
        
        aovMN <- t(apply(datMN, 2, function(varVn) {

            aovMod <- aov(varVn ~ facFcVn)
            pvaN <- summary(aovMod)[[1]][1, "Pr(>F)"]
            hsdMN <- TukeyHSD(aovMod)[["facFcVn"]]
            c(pvaN, c(hsdMN[, c("diff", "p adj")]))

        }))

        difVi <- 1:length(prwVc) + 1

        ## difference of the means for each pairwise comparison
        
        difMN <- aovMN[, difVi]
        colnames(difMN) <- paste0(varPfxC, prwVc, "_dif")

        ## correction for multiple testing
        
        aovMN <- aovMN[, -difVi, drop = FALSE]
        aovMN <- apply(aovMN, 2, function(pvaVn) p.adjust(pvaVn, method = adjC))

        ## significance coding (0 = not significant, 1 = significant)
        
        adjVn <- aovMN[, 1]
        sigVn <-  as.numeric(adjVn < thrN)

        aovMN <- aovMN[, -1, drop = FALSE]
        colnames(aovMN) <- paste0(varPfxC, prwVc, "_", adjC)

        aovSigMN <- aovMN < thrN
        mode(aovSigMN) <- "numeric"
        colnames(aovSigMN) <- paste0(varPfxC, prwVc, "_sig")

        ## final aggregated table

        resMN <- cbind(adjVn, sigVn, difMN, aovMN, aovSigMN)
        colnames(resMN)[1:2] <- paste0(varPfxC, c(adjC, "sig"))

        varDF <- cbind.data.frame(varDF, as.data.frame(resMN))

        ## graphic

        pdf(pdfC, onefile = TRUE)
        
        for(varI in 1:nrow(varDF)) {
            
            if(sum(aovSigMN[varI, ]) > 0) {
                
                varC <- rownames(varDF)[varI]

                boxplot(datMN[, varI] ~ facFcVn,
                        main = paste0(varC, " (", adjC, " = ", signif(adjVn[varI], 2), ")"))
                
                for(prwI in 1:length(prwVc)) {
                    
                    if(aovSigMN[varI, paste0(varPfxC, prwVc[prwI], "_sig")] == 1) {
                        
                        claVc <- unlist(strsplit(prwVc[prwI], ".", fixed = TRUE))
                        aovClaVl <- facFcVn %in% claVc
                        aovFc <- facFcVn[aovClaVl, drop = TRUE]
                        aovVc <- as.character(aovFc)
                        names(aovVc) <- rownames(samDF)[aovClaVl]
                        boxF(aovFc,
                             datMN[aovClaVl, varI],
                             paste0(varC, " (", adjC, " = ", signif(aovMN[varI, paste0(varPfxC, prwVc[prwI], "_", adjC)], 2), ")"),
                             aovVc)
                        
                    }
                       
                }
                
            }
            
        }

        dev.off()

        
    } else if(tesC == "kruskal") {
        

        ## getting the names of the pairwise comparisons 'class1.class2'
        
        nemMN <- posthoc.kruskal.nemenyi.test(datMN[, 1], facFcVn, "Tukey")[["p.value"]]
        nemVl <- c(lower.tri(nemMN, diag = TRUE))
        nemClaMC <- cbind(rownames(nemMN)[c(row(nemMN))][nemVl],
                          colnames(nemMN)[c(col(nemMN))][nemVl])
        nemNamVc <- paste0(nemClaMC[, 1], ".", nemClaMC[, 2])
        pfxNemVc <- paste0(varPfxC, nemNamVc)

        ## omnibus and post-hoc tests
        
        nemMN <- t(apply(datMN, 2, function(varVn) {

            pvaN <- kruskal.test(varVn ~ facFcVn)[["p.value"]]
            varNemMN <- posthoc.kruskal.nemenyi.test(varVn, facFcVn, "Tukey")[["p.value"]]
            c(pvaN, c(varNemMN))

        }))

        ## correction for multiple testing
        
        nemMN <- apply(nemMN, 2,
                       function(pvaVn) p.adjust(pvaVn, method = adjC))
        adjVn <- nemMN[, 1]
        sigVn <- as.numeric(adjVn < thrN)
        nemMN <- nemMN[, c(FALSE, nemVl)]
        colnames(nemMN) <- paste0(pfxNemVc, "_", adjC)

        ## significance coding (0 = not significant, 1 = significant)
        
        nemSigMN <- nemMN < thrN
        mode(nemSigMN) <- "numeric"
        colnames(nemSigMN) <- paste0(pfxNemVc, "_sig")

        ## difference of the medians for each pairwise comparison
        
        difMN <- sapply(1:nrow(nemClaMC), function(prwI) {
            prwVc <- nemClaMC[prwI, ]
            prwVi <- which(facFcVn %in% prwVc)
            prwFacFc <- factor(as.character(facFcVn)[prwVi], levels = prwVc)
            apply(datMN[prwVi, ], 2, function(varVn) -diff(as.numeric(tapply(varVn, prwFacFc, function(x) median(x, na.rm = TRUE)))))
        })
        colnames(difMN) <- gsub("_sig", "_dif", colnames(nemSigMN))

        ## final aggregated table
        
        resMN <- cbind(adjVn, sigVn, difMN, nemMN, nemSigMN)
        colnames(resMN)[1:2] <- paste0(varPfxC, c(adjC, "sig"))

        varDF <- cbind.data.frame(varDF, as.data.frame(resMN))

        ## graphic

        pdf(pdfC, onefile = TRUE)
        
        for(varI in 1:nrow(varDF)) {
            
            if(sum(nemSigMN[varI, ]) > 0) {
                
                varC <- rownames(varDF)[varI]

                boxplot(datMN[, varI] ~ facFcVn,
                        main = paste0(varC, " (", adjC, " = ", signif(adjVn[varI], 2), ")"))
                
                for(nemI in 1:length(nemNamVc)) {
                    
                    if(nemSigMN[varI, paste0(varPfxC, nemNamVc[nemI], "_sig")] == 1) {
                        
                        nemClaVc <- nemClaMC[nemI, ]
                        nemClaVl <- facFcVn %in% nemClaVc
                        nemFc <- facFcVn[nemClaVl, drop = TRUE]
                        nemVc <- as.character(nemFc)
                        names(nemVc) <- rownames(samDF)[nemClaVl]
                        boxF(nemFc,
                             datMN[nemClaVl, varI],
                             paste0(varC, " (", adjC, " = ", signif(nemMN[varI, paste0(varPfxC, nemNamVc[nemI], "_", adjC)], 2), ")"),
                             nemVc)
                        
                    }
                       
                }
                
            }
            
        }

        dev.off()
        
    }
    
    names(sigVn) <- rownames(varDF)
    sigSumN <- sum(sigVn, na.rm = TRUE)
    if(sigSumN) {
        cat("\nThe following ", sigSumN, " variable", ifelse(sigSumN > 1, "s", ""), " (", round(sigSumN / length(sigVn) * 100), "%) ", ifelse(sigSumN > 1, "were", "was"), " found significant at the ", thrN, " level:\n", sep = "")
        cat(paste(rownames(varDF)[sigVn > 0], collapse = "\n"), "\n", sep = "")
    } else
        cat("\nNo significant variable found at the selected ", thrN, " level\n", sep = "")

    options(stingsAsFactors = strAsFacL)

    return(varDF)

}


boxF <- function(xFc,
                 yVn,
                 maiC,
                 xVc) {
    
    boxLs <- boxplot(yVn ~  xFc,
                     main = maiC)
    
    outVn <- boxLs[["out"]]
    
    if(length(outVn)) {
        
        for(outI in 1:length(outVn)) {
            levI <- which(levels(xFc) == xVc[names(outVn)[outI]])
            text(levI,
                 outVn[outI],
                 labels = names(outVn)[outI],
                 pos = ifelse(levI == 2, 2, 4))
            }
        
    }
    
}
