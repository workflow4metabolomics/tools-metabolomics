## Etienne Thevenot
## CEA, MetaboHUB
## W4M Core Development Team
## etienne.thevenot@cea.fr
## 2015-05-30


heatmapF <- function(proMN,
                     obsDF,
                     feaDF,
                     disC,    ## dissimilarity
                     cutSamN, ## number of sample clusters
                     cutVarN, ## number of variable clusters
                     fig.pdfC,
                     corMetC, ## correlation method
                     aggMetC, ## agglomeration method
                     colC,    ## color scale
                     scaL,
                     cexN) {

    ncaN <- 14 ## Sample and variable name truncature for display

    if(aggMetC == "ward") {

        rvsLs <- R.Version()
        aggMetC <- paste0(aggMetC,
                          ifelse(as.numeric(rvsLs[["major"]]) > 3 ||
                                 as.numeric(rvsLs[["major"]]) == 3 && as.numeric(rvsLs[["minor"]]) > 0.3,
                                 ".D",
                                 ""))
        
    }
    
    if(disC %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
        
        obsHcl <- hclust(dist(proMN, method = disC),
                         method = aggMetC)
        
        feaHcl <- hclust(dist(t(proMN), method = disC),
                         method = aggMetC)
        
    } else if(disC == "1-cor") {
        
        obsHcl <- hclust(as.dist(1-cor(t(proMN),
                                       method = corMetC,
                                       use = "pairwise.complete.obs")),
                         method = aggMetC)

        feaHcl <- hclust(as.dist(1-cor(proMN,
                                       method = corMetC,
                                       use = "pairwise.complete.obs")),
                         method = aggMetC)

    } else if(disC == "1-abs(cor)") {

        obsHcl <- hclust(as.dist(1-abs(cor(t(proMN),
                                           method = corMetC,
                                           use = "pairwise.complete.obs"))),
                         method = aggMetC)

        feaHcl <- hclust(as.dist(1-abs(cor(proMN,
                                           method = corMetC,
                                           use = "pairwise.complete.obs"))),
                         method = aggMetC)

    }

    heaMN <- proMN <- proMN[obsHcl[["order"]], feaHcl[["order"]]]

    if(scaL)
        heaMN <- scale(heaMN)

    heaMN <- heaMN[, rev(1:ncol(heaMN)), drop = FALSE]

    switch(colC,
           blueOrangeRed = {
               imaPalVn <- colorRampPalette(c("blue", "orange", "red"),
                                            space = "rgb")(5)[1:5]
           },
           redBlackGreen = {
               imaPalVn <- colorRampPalette(c("red", "black", "green"),
                                            space = "rgb")(5)[1:5]
           })


    ## figure
    ##-------

    pdf(fig.pdfC,
        width = 14,
        height = 14)

    layout(matrix(1:4, nrow = 2),
           widths = c(1, 9), heights = c(1, 9))

    ## Color scale

    scaN <- length(imaPalVn)

    par(mar = c(0.6, 0.6, 0.6, 4.1))

    ylimVn <- c(0, scaN)
    ybottomVn <- 0:(scaN - 1)
    ytopVn <- 1:scaN

    plot(x = 0,
         y = 0,
         bty = "n",
         font.axis = 2,
         font.lab = 2,
         type = "n",
         xlim = c(0, 1),
         ylim = ylimVn,
         xlab = "",
         ylab = "",
         xaxs = "i",
         yaxs = "i",
         xaxt = "n",
         yaxt = "n")

    rect(xleft = 0.8,
         ybottom = ybottomVn,
         xright = 1,
         ytop = ytopVn,
         col = imaPalVn,
         border = NA)

    prtVn <- pretty(range(heaMN, na.rm = TRUE))
    axis(at = scaN / diff(range(prtVn)) * (prtVn - min(prtVn)),
         font = 2,
         font.axis = 2,
         labels = prtVn,
         las = 1,
         lwd = 2,
         lwd.ticks = 2,
         side = 4,
         xpd = TRUE)

    arrows(par("usr")[2],
           par("usr")[4],
           par("usr")[2],
           par("usr")[3],
           code = 0,
           lwd = 2,
           xpd = TRUE)

    ## Feature dendrogram

    par(mar = c(7.1, 0.6, 0, 0.1),
        lwd = 2)

    plot(rev(as.dendrogram(feaHcl)), horiz = TRUE,
         leaflab = "none",
         main = "", xaxs = "i", yaxs = "i",
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")

    revFeaHcl <- list(merge = cbind(feaHcl[["merge"]][, 2], feaHcl[["merge"]][, 1]),
                      height = feaHcl[["height"]],
                      order = rev(feaHcl[["order"]]),
                      labels = feaHcl[["labels"]])

    if(cutVarN > 1) {
        cluFeaVn <- cutree(revFeaHcl, k = cutVarN)[revFeaHcl[["order"]]]
        cutFeaVn <- which(abs(diff(cluFeaVn)) > 0)
        cutFeaTxtVn <- c(cutFeaVn[1] / 2, cutFeaVn + diff(c(cutFeaVn, length(cluFeaVn))) / 2) + 0.5
        cutFeaLinVn <- cutFeaVn + 0.5
        text(par("usr")[1] + 0.2 * diff(par("usr")[1:2]),
             cutFeaTxtVn,
             labels = unique(cluFeaVn),
             cex = 2,
             font = 2,
             las = 2)
    }

    ## Observation dendrogram

    par(mar = c(0.1, 0, 0.6, 7.1),
        lwd = 2)

    plot(as.dendrogram(obsHcl), leaflab = "none",
         main = "", xaxs = "i", yaxs = "i",
         yaxt = "n", xlab = "", ylab = "")

    if(cutSamN > 1) {
        cluObsVn <- cutree(obsHcl, k = cutSamN)[obsHcl[["order"]]]
        cutObsVn <- which(abs(diff(cluObsVn)) > 0)
        cutObsTxtVn <- c(cutObsVn[1] / 2, cutObsVn + diff(c(cutObsVn, length(cluObsVn))) / 2) + 0.5
        cutObsLinVn <- cutObsVn + 0.5
        text(cutObsTxtVn,
             0.8 * par("usr")[4],
             labels =  unique(cluObsVn),
             cex = 2,
             font = 2)
    }

    ## Heatmap

    par(mar = c(7.1, 0, 0, 7.1))

    image(x = 1:nrow(heaMN),
          y = 1:ncol(heaMN),
          z = round(heaMN),
          col = imaPalVn,
          font.axis = 2,
          font.lab = 2,
          xaxt = "n",
          yaxt = "n",
          xlab = "",
          ylab = "")

    obsOrdVc <- obsHcl[["labels"]][obsHcl[["order"]]]
    obsOrdLenVn <- sapply(obsOrdVc, nchar)
    obsOrdVc <- substr(obsOrdVc, 1, ncaN)
    obsOrdVc <- paste0(obsOrdVc, ifelse(obsOrdLenVn > ncaN, ".", ""), " ")

    mtext(obsOrdVc,
          at = 1:nrow(heaMN),
          cex = cexN,
          las = 2,
          side = 1)

    feaOrdVc <- feaHcl[["labels"]][feaHcl[["order"]]]
    feaOrdLenVn <- sapply(feaOrdVc, nchar)
    feaOrdVc <- substr(feaOrdVc, 1, ncaN)
    feaOrdVc <- paste0(" ", feaOrdVc, ifelse(feaOrdLenVn > ncaN, ".", ""))

    mtext(feaOrdVc,
          at = ncol(heaMN):1,
          cex = cexN,
          las = 2,
          side = 4)

    if(cutVarN > 1)
        abline(h = cutFeaLinVn)
    if(cutSamN > 1)
        abline(v = cutObsLinVn)

    box()

    dev.off()

    ## Returning
    ##----------

    if(cutSamN > 1)
        obsDF[, "heat_clust"] <- cutree(obsHcl, k = cutSamN)
    obsDF <- obsDF[obsHcl[["order"]], , drop = FALSE]

    if(cutVarN > 1)
        feaDF[, "heat_clust"] <- cutree(feaHcl, k = cutVarN)
    feaDF <- feaDF[feaHcl[["order"]], , drop = FALSE]

    return(invisible(list(proMN = proMN,
                          obsDF = obsDF,
                          feaDF = feaDF)))


} ## end of heatmapF
