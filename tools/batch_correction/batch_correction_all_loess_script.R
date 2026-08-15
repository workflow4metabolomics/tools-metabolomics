loessF <- function(datVn, qcaVi, preVi, spnN) {
    if (length(qcaVi) < 5) {
        return(predict(lm(datVn[qcaVi] ~ qcaVi),
            newdata = data.frame(qcaVi = preVi)
        ))
    } else {
        return(predict(
            loess(datVn[qcaVi] ~ qcaVi,
                control = loess.control(surface = "direct"),
                span = spnN
            ),
            newdata = data.frame(qcaVi = preVi)
        ))
    }

    ## Note:
    ##  the surface = 'direct' argument allows extrapolation
} ## loessF

plotBatchF <- function(datMN, samDF.arg, spnN.arg) {
    maiC <- switch(gsub("MN", "", deparse(substitute(datMN))),
        raw = "Raw",
        nrm = "Normalized"
    )

    colVc <- c(
        samp = "green4",
        biol = "green4",
        pool = "red",
        blan = "black",
        other = "yellow"
    )

    par(font = 2, font.axis = 2, font.lab = 2, lwd = 2, pch = 18)

    layout(matrix(c(1, 1, 2, 3), nrow = 2),
        widths = c(0.7, 0.3)
    )

    obsNamVc <- rownames(datMN)

    obsColVc <- sapply(
        substr(samDF.arg[, args$sample_type_col_name], 1, 4),
        function(typC) {
            ifelse(typC %in% names(colVc), colVc[typC], colVc["other"])
        }
    )

    ## Graphic 1: Sum of intensities for each sample

    par(mar = c(3.6, 3.6, 3.1, 0.6))

    batTab <- table(samDF.arg[, args$batch_col_name])

    sumVn <- rowSums(datMN, na.rm = TRUE)

    plot(sumVn,
        cex = 1.2,
        col = obsColVc,
        pch = 18,
        xaxs = "i",
        xlab = "",
        ylab = ""
    )

    mtext("Injection order",
        line = 2.2,
        side = 1
    )
    mtext("Sum of variable intensities",
        line = 2.2,
        side = 2
    )

    mtext(maiC, cex = 1.2, line = 1.5, side = 3)

    abline(
        v = cumsum(batTab) + 0.5,
        col = "red"
    )

    mtext(names(batTab),
        at = batTab / 2 + c(0, cumsum(batTab[-length(batTab)]))
    )

    obsColVuc <- obsColVc[sort(unique(names(obsColVc)))]

    text(rep(batTab[1], times = length(obsColVuc)),
        par("usr")[3] + (0.97 - length(obsColVuc) * 0.03 + 1:length(obsColVuc) * 0.03) * diff(par("usr")[3:4]),
        col = obsColVuc,
        font = 2,
        labels = names(obsColVuc),
        pos = 2
    )

    for (batC in names(batTab)) {
        batSeqVi <- which(samDF.arg[, args$batch_col_name] == batC)
        batPooVi <- intersect(
            batSeqVi,
            which(samDF.arg[, args$sample_type_col_name] == "pool")
        )
        batSamVi <- intersect(
            batSeqVi,
            which(samDF.arg[, args$sample_type_col_name] == "sample")
        )
        if (length(batPooVi)) {
            lines(batSeqVi,
                loessF(sumVn, batPooVi, batSeqVi, spnN = spnN.arg),
                col = colVc["pool"]
            )
        }
        if (length(batSamVi)) {
            lines(batSeqVi,
                loessF(sumVn, batSamVi, batSeqVi, spnN = spnN.arg),
                col = colVc["samp"]
            )
        }
    }

    ## Graphics 2 and 3 (right): PCA score plots of components 1-4

    radVn <- seq(0, 2 * pi, length.out = 100)
    epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16

    pcaMN <- datMN

    if (any(is.na(pcaMN))) {
        minN <- min(pcaMN, na.rm = TRUE)
        pcaMN[is.na(pcaMN)] <- minN
    }

    pcaLs <- opls(pcaMN, predI = 4, algoC = "svd")
    tMN <- getScoreMN(pcaLs)
    vRelVn <- pcaLs@modelDF[, "R2X"]

    n <- nrow(tMN)
    hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))

    hotFisN <- hotN * qf(0.95, 2, n - 2)

    pcsLs <- list(c(1, 2), c(3, 4))

    par(mar = c(3.6, 3.6, 0.6, 1.1))

    for (pcsN in 1:length(pcsLs)) {
        pcsVn <- pcsLs[[pcsN]]

        tcsMN <- tMN[, pcsVn]

        micMN <- solve(cov(tcsMN))

        n <- nrow(tMN)
        hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))

        hotFisN <- hotN * qf(0.95, 2, n - 2)

        hotVn <- apply(
            tcsMN,
            1,
            function(x) 1 - pf(1 / hotN * t(as.matrix(x)) %*% micMN %*% as.matrix(x), 2, n - 2)
        )

        obsHotVi <- which(hotVn < 0.05)

        xLabC <- paste("t",
            pcsVn[1],
            "(",
            round(vRelVn[pcsVn[1]] * 100),
            "%)",
            sep = ""
        )

        yLabC <- paste("t",
            pcsVn[2],
            "(",
            round(vRelVn[pcsVn[2]] * 100),
            "%)",
            sep = ""
        )

        xLimVn <- c(-1, 1) * max(sqrt(var(tcsMN[, 1]) * hotFisN), max(abs(tcsMN[, 1])))
        yLimVn <- c(-1, 1) * max(sqrt(var(tcsMN[, 2]) * hotFisN), max(abs(tcsMN[, 2])))

        plot(tcsMN,
            main = "",
            type = "n",
            xlab = "",
            ylab = "",
            xlim = xLimVn,
            ylim = yLimVn
        )

        mtext(xLabC,
            line = 2.2,
            side = 1
        )
        mtext(yLabC,
            line = 2.2,
            side = 2
        )

        par(lwd = 1)

        abline(
            v = axTicks(1),
            col = "grey"
        )

        abline(
            h = axTicks(2),
            col = "grey"
        )

        abline(v = 0)
        abline(h = 0)

        lines(
            sqrt(var(tcsMN[, 1]) * hotFisN) * cos(radVn),
            sqrt(var(tcsMN[, 2]) * hotFisN) * sin(radVn)
        )

        points(tcsMN,
            col = obsColVc,
            pch = 18
        )

        if (length(obsHotVi)) {
            text(tcsMN[obsHotVi, 1],
                tcsMN[obsHotVi, 2],
                col = obsColVc[obsHotVi],
                labels = obsNamVc[obsHotVi],
                pos = 3
            )
        }
    } ## for(pcsN in 1:length(pcsLs)) {

    return(invisible(list(
        sumVn = sumVn,
        tcsMN = tcsMN
    )))
} ## plotBatchF

shiftBatchCorrectF <- function(rawMN.arg,
                               samDF.arg,
                               refC.arg,
                               spnN.arg) {
    cat("\nReference observations are: ", refC.arg, "\n")

    ## computing median off all pools (or samples) for each variable

    refMeaVn <- apply(
        rawMN.arg[samDF.arg[, args$sample_type_col_name] == refC.arg, ],
        2,
        function(feaRefVn) mean(feaRefVn, na.rm = TRUE)
    )

    ## splitting data and sample metadata from each batch

    batRawLs <- split(as.data.frame(rawMN.arg),
        f = samDF.arg[, args$batch_col_name]
    )
    batRawLs <- lapply(batRawLs, function(inpDF) as.matrix(inpDF))

    batSamLs <- split(as.data.frame(samDF.arg),
        f = samDF.arg[, args$batch_col_name]
    )

    ## checking extrapolation: are there pools at the first and last observations of each batch

    if (refC.arg == args$sample_type_tags$pool) {
        pooExtML <- matrix(FALSE,
            nrow = 2, ncol = length(batRawLs),
            dimnames = list(c("first", "last"), names(batRawLs))
        )
        for (batC in names(batSamLs)) {
            batSamTypVc <- batSamLs[[batC]][, args$sample_type_col_name]
            pooExtML["first", batC] <- head(batSamTypVc, 1) == args$sample_type_tags$pool
            pooExtML["last", batC] <- tail(batSamTypVc, 1) == args$sample_type_tags$pool
        }
        if (!all(c(pooExtML))) {
            cat("\nWarning: Pools are missing at the first and/or last position of the following batches:\n")
            pooExtBatVi <- which(!apply(pooExtML, 2, all))
            for (i in 1:length(pooExtBatVi)) {
                cat(names(pooExtBatVi)[i], ": ",
                    paste(rownames(pooExtML)[!pooExtML[, pooExtBatVi[i]]], collapse = ", "), "\n",
                    sep = ""
                )
            }
            cat("Extrapolating loess fits for these batches may result in inaccurate modeling!\n")
        }
    }

    ## normalizing

    nrmMN <- NULL ## normalized data matrix to be computed

    cat("\nProcessing batch:")

    for (batC in names(batRawLs)) { ## processing each batch individually

        cat("\n", batC)

        batRawMN <- batRawLs[[batC]]
        batSamDF <- batSamLs[[batC]]

        batAllVi <- 1:nrow(batRawMN)

        batRefVi <- which(batSamDF[, args$sample_type_col_name] == refC.arg)

        if (length(batRefVi) < 5) {
            cat("\nWarning: less than 5 '", refC.arg, "'; linear regression will be performed instead of loess regression for this batch\n", sep = "")
        }

        ## prediction of the loess fit

        batLoeMN <- apply(
            batRawMN,
            2,
            function(rawVn) loessF(rawVn, batRefVi, batAllVi, spnN = spnN.arg)
        )

        ## normalization

        batLoeMN[batLoeMN <= 0] <- NA

        batNrmMN <- batRawMN / batLoeMN

        nrmMN <- rbind(
            nrmMN,
            batNrmMN
        )
    }

    cat("\n")

    nrmMN <- sweep(nrmMN, MARGIN = 2, STATS = refMeaVn, FUN = "*")

    return(nrmMN)
} ## shiftBatchCorrectF
