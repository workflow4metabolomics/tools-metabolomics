################################################################################################
# ANALYSES FOR QUALITY CONTROL                                                                 #
#                                                                                              #
# Author: Melanie PETERA                                                                       #
# User: Galaxy                                                                                 #
# Starting date: 04-09-2014                                                                    #
# V-1.0: Restriction of old filter script to CV filter                                         #
# V-1.1: Addition of data check                                                                #
# V-1.2: Substitution of deletion by addition of indicator variable                            #
# V-1.3: Handling special characters                                                           #
#                                                                                              #
#                                                                                              #
# Input files: dataMatrix ; sampleMetadata ; variableMetadata                                  #
# Output files: dataMatrix ; sampleMetadata ; variableMetadata                                 #
#                                                                                              #
################################################################################################

# Parameters (for dev)
if(FALSE){

    ion.file.in <- "test/ressources/inputs/ex_data_IONS.txt" #tab file
    meta.samp.file.in <- "test/ressources/inputs/ex_data_PROTOCOLE1.txt" #tab file
    meta.ion.file.in <- "test/ressources/inputs/ex_data_METAION.txt" #tab file

    ## ion.file.out <- "test/ressources/outputs/QCtest_ex_data_IONS.txt" #tab file
    meta.samp.file.out <- "test/ressources/outputs/QCtest_ex_data_PROTOCOLE1.txt" #tab file
    meta.ion.file.out <- "test/ressources/outputs/QCtest_ex_data_METAION.txt" #tab file

    CV <- TRUE ; if(CV){Compa<-TRUE;seuil<-1.25}else{Compa<-NULL;seuil<-NULL}

    poolAsPool1L <- FALSE

    if(FALSE) { ## Sacuri dataset

        ## 'example' input dir
        exaDirInpC <- "example/input"

        ion.file.in <- file.path(exaDirInpC, "dataMatrix.tsv")
        meta.samp.file.in <- file.path(exaDirInpC, "sampleMetadata.tsv")
        meta.ion.file.in <- file.path(exaDirInpC, "variableMetadata.tsv")

        poolAsPool1L <- FALSE

        ## 'example' output dir
        exaDirOutC <- gsub("input", "output", exaDirInpC)

        mata.samp.file.out <- file.path(exaDirOutC, "sampleMetadata.tsv")
        meta.ion.file_out <- file.path(exaDirOutC, "variableMetadata.tsv")
        fig.out <- file.path(exaDirOutC, "figure.pdf")
        log.out <- file.path(exaDirOutC, "information.txt")

        stopifnot(file.exists(exaDirOutC))

    }

}

QualityControl <- function(ion.file.in, meta.samp.file.in, meta.ion.file.in,
                           CV, Compa, seuil, poolAsPool1L,
                           ion.file.out, meta.samp.file.out, meta.ion.file.out, fig.out, log.out){
  # This function allows to analyse data to check its quality
  # It needs 3 datasets: the data matrix, the variables' metadata, the samples' metadata.
  # It generates 3 new datasets corresponding to the 3 inputs with additional columns.
  #
  # Parameters:
  # - xxx.in: input files' names
  # - xxx.out: output files' names
  # - CV: CV calculation yes/no
  # | > Compa: comparing pool and sample CVs (TRUE) or simple pool CV calculation (FALSE)
  # | > seuil: maximum ratio tolerated between pool and sample CVs or maximum pool CV


# Input -----------------------------------------------------------------------------------

ion.data <- read.table(ion.file.in,sep="\t",header=TRUE,check.names=FALSE, stringsAsFactors = FALSE)
meta.samp.data <- read.table(meta.samp.file.in,sep="\t",header=TRUE,check.names=FALSE, stringsAsFactors = FALSE)
meta.ion.data <- read.table(meta.ion.file.in,sep="\t",header=TRUE,check.names=FALSE, stringsAsFactors = FALSE)

# Error vector
err.stock <- "\n"

# Table match check
table.check <- match3(ion.data,meta.samp.data,meta.ion.data)
check.err(table.check)

# StockID
samp.id <- stockID(ion.data,meta.samp.data,"sample")
ion.data <- samp.id$dataMatrix
meta.samp.data <- samp.id$Metadata
samp.id <- samp.id$id.match


# Function 1: CV calculation --------------------------------------------------------------
# Allows to class ions according to the Coefficient of Variation (CV):
# Compa=TRUE:
# 	CV of pools and CV of samples are compared (ration between pools' one and samples' one)
# 	and confronted to a given ration.
# Compa=FALSE:
# 	only CV of pools are considered ; compared to a given threshold

if(CV){

  # Checking the sampleType variable
  if(is.null(meta.samp.data$sampleType)){
    err.stock <- c(err.stock,"\n-------",
                   "\nWarning : no 'sampleType' variable detected in sample meta-data !",
                   "\nCV can not be calculated.\n-------\n")
  }else{
    if(!("pool"%in%levels(factor(meta.samp.data$sampleType)))){
      err.stock <- c(err.stock,"\n-------",
                     "\nWarning : no 'pool' detected in 'sampleType' variable (sample meta-data) !",
                     "\nCV can not be calculated.\n-------\n")
    }else{
      if((!("sample"%in%levels(factor(meta.samp.data$sampleType))))&(Compa)){
        err.stock <- c(err.stock,"\n-------",
                       "\nWarning : no 'sample' detected in 'sampleType' variable (sample meta-data) !",
                       "\nCV can not be calculated.\n-------\n")
      }else{

  # Statement
  tmp.ion <- data.frame(CV.ind=rep(NA,nrow(ion.data)),CV.samp=rep(NA,nrow(ion.data)),
                        CV.pool=rep(NA,nrow(ion.data)),ion.data,stringsAsFactors=FALSE)
  # CV samples
  tmp.samp <- which(colnames(tmp.ion)%in%meta.samp.data[which(meta.samp.data$sampleType=="sample"),1])
  tmp.ion$CV.samp <- apply(tmp.ion[,tmp.samp],1,function(x)sd(x, na.rm = TRUE)) / rowMeans(tmp.ion[,tmp.samp], na.rm = TRUE)
  tmp.ion$CV.samp[which(apply(tmp.ion[,tmp.samp],1,function(x)sd(x, na.rm = TRUE))==0)] <- 0
  # CV pools
  tmp.samp <- which(colnames(tmp.ion)%in%meta.samp.data[which(meta.samp.data$sampleType=="pool"),1])
  tmp.ion$CV.pool <- apply(tmp.ion[,tmp.samp],1,function(x)sd(x, na.rm = TRUE)) / rowMeans(tmp.ion[,tmp.samp], na.rm = TRUE)
  tmp.ion$CV.pool[which(apply(tmp.ion[,tmp.samp],1,function(x)sd(x, na.rm = TRUE))==0)] <- 0
  # CV indicator
  if(Compa){tmp.ion$CV.ind <- ifelse((tmp.ion$CV.pool)/(tmp.ion$CV.samp)>seuil,0,1)
  }else{tmp.ion$CV.ind <- ifelse((tmp.ion$CV.pool)>seuil,0,1)}
  # Addition of new columns in meta.ion.data
  if(Compa){tmp.ion<-tmp.ion[,c(4,2,3,1,1)]}else{tmp.ion<-tmp.ion[,c(4,3,1,1)]}
  tmp.ion[,ncol(tmp.ion)] <- 1:nrow(tmp.ion)
  meta.ion.data <- merge(x=meta.ion.data,y=tmp.ion,by.x=1,by.y=1)
  meta.ion.data <- meta.ion.data[order(meta.ion.data[,ncol(meta.ion.data)]),][,-ncol(meta.ion.data)]
  rownames(meta.ion.data) <- NULL

  rm(tmp.ion,tmp.samp)

      }}}

} # end if(CV)

## complementary metrics (ET)

datMN <- t(as.matrix(ion.data[, -1]))
colnames(datMN) <- ion.data[, 1]
datMN <- datMN[, meta.ion.data[, 1]] ## in case meta.ion.data has been re-ordered during the CV = TRUE computations
quaLs <- qualityMetricsF(datMN,
                         meta.samp.data,
                         meta.ion.data,
                         poolAsPool1L,
                         fig.out,
                         log.out)
meta.samp.data <- quaLs[["samDF"]]
meta.ion.data <- quaLs[["varDF"]]


# Output ----------------------------------------------------------------------------------

# Getting back original identifiers
id.ori <- reproduceID(ion.data,meta.samp.data,"sample",samp.id)
ion.data <- id.ori$dataMatrix
meta.samp.data <- id.ori$Metadata


# Error checking
if(length(err.stock)>1){
  stop(err.stock)
}else{

## write.table(ion.data, ion.file.out, sep="\t", row.names=FALSE, quote=FALSE)
write.table(meta.samp.data, meta.samp.file.out, sep="\t", row.names=FALSE, quote=FALSE)
write.table(meta.ion.data, meta.ion.file.out, sep="\t", row.names=FALSE, quote=FALSE)

}


} # end of QualityControl function


# Typical function call
# QualityControl(ion.file.in, meta.samp.file.in, meta.ion.file.in,
#       CV, Compa, seuil,
#       ion.file.out, meta.samp.file.out, meta.ion.file.out)


qualityMetricsF <- function(datMN,
                            samDF,
                            varDF,
                            pooAsPo1L = TRUE,
                            fig.pdfC = NULL,
                            log.txtC = NULL) {

    optWrnN <- options()$warn
    options(warn = -1)


    ##------------------------------
    ## Functions
    ##------------------------------


    allDigF <- function (string) { ## from the Hmisc package (all.digits)
        k <- length(string)
        result <- logical(k)
        for (i in 1:k) {
            st <- string[i]
            ls <- nchar(st)
            ex <- substring(st, 1:ls, 1:ls)
            result[i] <- all(match(ex, c("0", "1", "2", "3", "4",
                                         "5", "6", "7", "8", "9", "."), nomatch = 0) > 0)
        }
        result
    }

    datPloF <- function() { ## ploting data matrix

        thrVn <- c(pvalue=0.001,
                   poolCv=0.3)

        ## Constants

        marLs <- list(dri = c(2.1, 2.6, 1.1, 1.1),
                      ima = c(1.1, 2.6, 4.1, 1.1),
                      msd = c(2.1, 2.6, 1.1, 0.6),
                      sam = c(3.1, 3.6, 1.1, 0.6),
                      pca = c(2.6, 3.6, 1.1, 0.6),
                      sca = c(1.1, 4.1, 4.1, 0.6),
                      tit = c(0.1, 0.6, 1.1, 0.6))
        palHeaVc <- rev(rainbow(ceiling(256 * 1.5))[1:256])

        ## Functions

        axiPreF <- function(valVn,
                            lenN) {

            if(NA %in% valVn) {
                warning("NA in valVn")
                valVn <- as.vector(na.omit(valVn))
            }

            if(lenN < length(valVn))
                stop("The length of in vector must be inferior to the length of the length parameter.")

            if(length(valVn) < lenN)
                valVn <- seq(from = min(valVn), to = max(valVn), length.out = lenN)

            preValVn <- pretty(valVn)

            preLabVn <- preAtVn <- c()

            for(n in 1:length(preValVn))
                if(min(valVn) < preValVn[n] && preValVn[n] < max(valVn)) {
                    preLabVn <- c(preLabVn, preValVn[n])
                    preAtVn <- c(preAtVn, which(abs(valVn - preValVn[n]) == min(abs(valVn - preValVn[n])))[1])
                }

            return(list(atVn = preAtVn,
                        labVn = preLabVn))

        }

        colF <- function(vecVn)
            sapply(vecVn,
                   function(outN) {
                       if(outN < ploRgeVn[1])
                           return(palHeaVc[1])
                       else if(outN > ploRgeVn[2])
                           return(palHeaVc[256])
                       else return(palHeaVc[round((outN - ploRgeVn[1]) / diff(ploRgeVn) * 256 + 1)])})

        obsColF <- function(typVc) {

            ## available color palette
            palVc <- palette()

            ## colors for common types are set aside
            palVc <- palVc[!(palVc %in% c("black", "red", "green3"))]

            ## filling in the types with dedicated colors
            samTypVc <- sort(unique(samDF[, "sampleType"]))
            samColVc <- character(length(samTypVc))
            if("blank" %in% samTypVc)
                samColVc[grepl("blank", samTypVc)] <- "black"
            if("pool" %in% samTypVc)
                samColVc[grepl("pool", samTypVc)] <- "red"
            if("sample" %in% samTypVc)
                samColVc[grepl("sample", samTypVc)] <- "green4"

            ## filling in the other types
            palColI <- 1
            palColMaxI <- length(palVc)

            while(any(samColVc == "")) {
                typToColI <- which(samColVc == "")[1]
                if(palColI <= palColMaxI)
                    samColVc[typToColI] <- palVc[palColI]
                else
                    samColVc[typToColI] <- "gray"
                palColI <- palColI + 1
            }

            names(samColVc) <- samTypVc

            samColVc[typVc]

        }

        par(font = 2,
            font.axis = 2,
            font.lab = 2,
            pch=18)

        layout(matrix(c(1, 3, 4, 5, 5,
                        1, 7, 7, 7, 6,
                        2, 7, 7, 7, 6),
                      byrow = TRUE,
                      nrow = 3),
               heights = c(1.8, 1.2, 2.5),
               widths = c(3.5, 1.8, 2.8, 1, 0.8))

        ## Colors
        ##-------

        if("sampleType" %in% colnames(samDF)) {
            obsColVc <- obsColF(samDF[, "sampleType"])
        } else
            obsColVc <- rep("black", nrow(samDF))

        ## PCA and Hotelling ellipse
        ##--------------------------

        vVn <- getPcaVarVn(ropLs)
        vRelVn <- vVn / ncol(datMN)

        par(mar = marLs[["pca"]])

        plot(ropScoreMN,
             type = "n",
             xlab = "",
             ylab = "",
             xlim = range(ropScoreMN[, 1]) * 1.1)
        mtext(paste("t1 (", round(vRelVn[1] * 100), "%)", sep = ""),
              cex = 0.7,
              line = 2,
              side = 1)
        mtext(paste("t2 (", round(vRelVn[2] * 100), "%)", sep = ""),
              cex = 0.7,
              las = 0,
              line = 2,
              side = 2)
        abline(h = 0, lty = "dashed")
        abline(v = 0, lty = "dashed")
        radVn <- seq(0, 2 * pi, length.out = 100)

        hotFisN <- hotN * qf(1 - thrVn["pvalue"], 2, n - 2)
        lines(sqrt(var(ropScoreMN[, 1]) * hotFisN) * cos(radVn),
              sqrt(var(ropScoreMN[, 2]) * hotFisN) * sin(radVn))

        text(ropScoreMN[, 1],
             ropScoreMN[, 2],
             cex = 0.7,
             col = obsColVc,
             labels = rownames(datMN))

        if("sampleType" %in% colnames(samDF)) {
            obsColVuc <- obsColVc[sort(unique(names(obsColVc)))]
            legOrdVc <- c("blank", paste0("pool", 8:1), "pool", "other", "sample")
            obsColVuc <- obsColVuc[legOrdVc[legOrdVc %in% names(obsColVuc)]]

            text(rep(par("usr")[1], times = length(obsColVuc)),
                 par("usr")[3] + (0.97 - length(obsColVuc) * 0.03 + 1:length(obsColVuc) * 0.03) * diff(par("usr")[3:4]),
                 col = obsColVuc,
                 font = 2,
                 labels = names(obsColVuc),
                 pos = 4)
        }

        ## Missing/low intensities and decile values
        ##------------------------------------------

        par(mar = marLs[["sam"]])

        plot(missZscoVn,
             deciZscoMaxVn,
             type = "n",
             xlab = "",
             ylab = "",
             xlim = c(min(missZscoVn),
                 max(missZscoVn) + 0.5))
        mtext("amount of missing values (z-score)",
              cex = 0.7,
              line = 2,
              side = 1)
        mtext("deciles (zscore)",
              cex = 0.7,
              las = 0,
              line = 2,
              side = 2)
        abline(h = qnorm(1 - thrVn["pvalue"] / 2) * c(-1, 1), lty = "dashed")
        abline(v = qnorm(1 - thrVn["pvalue"] / 2) * c(-1, 1), lty = "dashed")
        text(missZscoVn,
             deciZscoMaxVn,
             cex = 0.7,
             col = obsColVc,
             labels = rownames(datMN))

        ## tit: Title
        ##-----------

        par(mar = marLs[["tit"]])
        plot(0:1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
        text(1.5, 1, cex = 1.3, labels = "Quality Metrics")
        text(1, 0.85, adj=0, cex = 1.1, labels = paste0("NAs: ",
                                            round(length(which(is.na(c(datMN)))) / cumprod(dim(datMN))[2] * 100), "%"))
        text(1, 0.75, adj=0, cex = 1.1, labels = paste0("0 values: ",
                                            round(sum(abs(datMN) < epsN, na.rm=TRUE) / cumprod(dim(datMN))[2] * 100, 2), "%"))
        text(1, 0.65, adj=0, cex = 1.1, labels = paste0("min: ", signif(min(datMN, na.rm=TRUE), 2)))
        text(1, 0.55, adj=0, cex = 1.1, labels = paste0("median: ", signif(median(datMN, na.rm=TRUE), 2)))
        text(1, 0.45, adj=0, cex = 1.1, labels = paste0("mean: ", signif(mean(datMN, na.rm=TRUE), 2)))
        text(1, 0.35, adj=0, cex = 1.1, labels = paste0("max: ", signif(max(datMN, na.rm=TRUE), 2)))
        if("sampleType" %in% colnames(samDF) &&
           "pool" %in% samDF[, "sampleType"]) {
            poolCvNanVl <- is.nan(varDF[, "pool_CV"])
            text(1,
                 0.25,
                 adj=0, cex = 1.1,
                 labels = paste0("pool CV < ",
                     round(thrVn["poolCv"] * 100), "%: ",
                     round(sum(varDF[!poolCvNanVl, "pool_CV", drop = FALSE] < thrVn["poolCv"]) / nrow(varDF) * 100),
                     "%"))
            }

        text(1, 0.1, adj=0, labels = paste0("Thresholds used in plots:"))
        text(1, 0, adj=0, labels = paste0("  p-value = ", thrVn["pvalue"]))

        ## dri: Analytical drift
        ##----------------------

        par(mar = marLs[["dri"]])

        ## ordering

        driDatMN <- datMN
        driSamDF <- samDF

        driSamDF[, "ordIniVi"] <- 1:nrow(driDatMN)

        if("injectionOrder" %in% colnames(driSamDF)) {
            if("batch" %in% colnames(driSamDF))
                ordVi <- order(driSamDF[, "batch"],
                               driSamDF[, "injectionOrder"])
            else
                ordVi <- order(driSamDF[, "injectionOrder"])
        } else
            ordVi <- 1:nrow(driDatMN)

        driDatMN <- driDatMN[ordVi, ]
        driSamDF <- driSamDF[ordVi, ]

        driColVc <- rep("black", nrow(driDatMN))
        if("sampleType" %in% colnames(driSamDF))
            driColVc <- obsColF(driSamDF[, "sampleType"])

        plot(rowSums(driDatMN, na.rm=TRUE),
             col = driColVc,
             pch = 18,
             xlab = "",
             ylab = "")

        mtext("injection order",
              cex = 0.7,
              line = 2,
              side = 1)

        mtext("Sum of intens. for all variables",
              cex = 0.7,
              line = 2,
              side = 2)

        ## msd: Sd vs Mean plot
        ##---------------------

        par(mar = marLs[["msd"]])
        plot(apply(datMN, 2, function(y) mean(y, na.rm = TRUE)),
             apply(datMN, 2, function(y) sd(y, na.rm = TRUE)),
             pch = 18,
             xlab = "",
             ylab = "")
        mtext("mean",
              cex = 0.7,
              line = 2,
              side = 1)
        mtext("sd",
              cex = 0.7,
              line = 2,
              side = 2)

        ## sca-6: Color scale
        ##-------------------

        par(mar = marLs[["sca"]])

        ylimVn <- c(0, 256)
        ybottomVn <- 0:255
        ytopVn <- 1:256

        plot(x = 0,
             y = 0,
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

        rect(xleft = 0,
             ybottom = ybottomVn,
             xright = 1,
             ytop = ytopVn,
             col = palHeaVc,
             border = NA)

        eval(parse(text = paste("axis(at = axiPreF(c(ifelse(min(datMN, na.rm = TRUE) == -Inf, yes = 0, no = min(datMN, na.rm = TRUE)) , max(datMN, na.rm = TRUE)), 256)$atVn,
             font = 2,
             font.axis = 2,
             labels = axiPreF(c(ifelse(min(datMN, na.rm = TRUE) == -Inf, yes = 0, no = min(datMN, na.rm = TRUE)), max(datMN, na.rm = TRUE)), 256)$labVn,
             las = 1,
             lwd = 2,
             lwd.ticks = 2,
             side = 2,
             xpd = TRUE)", sep = "")))

        arrows(par("usr")[1],
               par("usr")[4],
               par("usr")[1],
               par("usr")[3],
               code = 0,
               lwd = 2,
               xpd = TRUE)

        ## ima: Image
        ##-----------

        par(mar = marLs[["ima"]])

        ploRgeVn <- range(datMN, na.rm = TRUE)

        imaMN <- t(datMN)[, rev(1:nrow(datMN)), drop = FALSE]

        image(x = 1:nrow(imaMN),
              y = 1:ncol(imaMN),
              z = imaMN,
              col = palHeaVc,
              font.axis = 2,
              font.lab = 2,
              xaxt = "n",
              yaxt = "n",
              xlab = "",
              ylab = "")

        if(length(rownames(datMN)) == 0) {
            rowNamVc <- rep("", times = nrow(datMN))
        } else
            rowNamVc <- rownames(datMN)

        if(length(colnames(datMN)) == 0) {
            colNamVc <- rep("", times = ncol(datMN))
        } else
            colNamVc <- colnames(datMN)

        xlaVc <- paste(paste(rep("[", 2),
                             c(1, nrow(imaMN)),
                             rep("] ", 2),
                             sep = ""),
                       rep("\n", times = 2),
                       c(colNamVc[1], tail(colNamVc, 1)),
                       sep = "")

        for(k in 1:2)
            axis(side = 3,
                 hadj = c(0, 1)[k],
                 at = c(1, nrow(imaMN))[k],
                 cex = 0.8,
                 font = 2,
                 labels = xlaVc[k],
                 line = -0.5,
                 tick = FALSE)


        ylaVc <- paste(paste(rep("[", times = 2),
                             c(ncol(imaMN), 1),
                             rep("]", times = 2),
                             sep = ""),
                       rep("\n", times = 2),
                       c(tail(rowNamVc, 1), rowNamVc[1]),
                       sep = "")

        for(k in 1:2)
            axis(side = 2,
                 at = c(1, ncol(imaMN))[k],
                 cex = 0.8,
                 font = 2,
                 hadj = c(0, 1)[k],
                 labels = ylaVc[k],
                 las = 0,
                 line = -0.5,
                 lty = "blank",
                 tick = FALSE)

        box(lwd = 2)


    }


    zScoreF <- function(x) {
        sdxN <- sd(x, na.rm = TRUE)
        if(sdxN < epsN)
            return(rep(0, length(x)))
        else
            return((x - mean(x, na.rm = TRUE)) / sdxN)
    }


    ## Option
    ##-------

    strAsFacL <- options()$stringsAsFactors
    options(stingsAsFactors = FALSE)

    ## Constants
    ##----------

    modNamC <- "Quality Metrics" ## module name
    
    epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16


    ##------------------------------
    ## Start
    ##------------------------------

    if(!is.null(log.txtC)) {
        
        sink(log.txtC)
        
        cat("\nStart of the '", modNamC, "' Galaxy module call: ",
            format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

    }

    
    ## Checking the numerical type of the dataMatrix
    ##----------------------------------------------


    if(mode(datMN) != "numeric") {
        sink()
        stop("dataMatrix is not of numeric type;\ncheck your tables with the Check Format module\n",
             call. = FALSE)
    }

    
    ## Re-ordering dataMatrix samples if need (internally only)
    ##---------------------------------------------------------

    
    if(!identical(rownames(datMN), samDF[, 1])) {

        cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   WARNING   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\nThe sample order is not identical in dataMatrix and sampleMetadata;\nRe-ordering of the dataMatrix samples will be performed internally in this module\nfor the computation of the metrics,\nwithout changing the orders in the sampleMetadata output;\n\nTo get a re-ordered dataMatrix as output, please use the Check Format module\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

        datMN <- datMN[samDF[, 1], , drop = FALSE]
        
    }

    
    ## Description
    ##------------

    cat("\n\nData description:\n\n", sep = "")
    cat("observations:", nrow(datMN), "\n")
    cat("variables:", ncol(datMN), "\n")
    cat("missing:", sum(is.na(datMN)), "\n")
    cat("0 values (%):",
        sum(abs(datMN) < epsN, na.rm=TRUE) / cumprod(dim(datMN))[2] * 100, "\n")
    cat("min:", min(datMN, na.rm=TRUE), "\n")
    cat("mean:", signif(mean(datMN, na.rm=TRUE), 2), "\n")
    cat("median:", signif(median(datMN, na.rm=TRUE), 2), "\n")
    cat("max:", signif(max(datMN, na.rm=TRUE), 2), "\n")

    if("sampleType" %in% colnames(samDF)) {
        cat("\nSample types:\n", sep = "")
        print(table(samDF[, "sampleType"]))
        cat("\n", sep="")
    }


    ##------------------------------
    ## Variable metrics
    ##------------------------------


    ## 'blank' observations

    if("sampleType" %in% colnames(samDF) && "blank" %in% samDF[, "sampleType"]) {

        cat("\nVariables: Blank mean, sd, and CV\n", sep="")

        blkVl <- samDF[, "sampleType"] == "blank"

        if(sum(blkVl) == 1)
            varDF[, "blank_mean"] <- datMN[blkVl, ]
        else
            varDF[, "blank_mean"] <- apply(datMN[blkVl, , drop=FALSE], 2, function(varVn) mean(varVn, na.rm=TRUE))

        if(sum(blkVl) == 1)
            varDF[, "blank_sd"] <- rep(0, nrow(varDF))
        else
            varDF[, "blank_sd"] <- apply(datMN[blkVl, , drop=FALSE], 2, function(varVn) sd(varVn, na.rm=TRUE))

        varDF[, "blank_CV"] <- varDF[, "blank_sd"] / varDF[, "blank_mean"]

    }


    ## 'sample' observations

    if("sampleType" %in% colnames(samDF) && "sample" %in% samDF[, "sampleType"]) {

        cat("\nVariables: Sample mean, sd, and CV\n", sep="")

        samVl <- samDF[, "sampleType"] == "sample"

        if(sum(samVl) == 1)
            varDF[, "sample_mean"] <- datMN[samVl, ]
        else
            varDF[, "sample_mean"] <- apply(datMN[samVl, , drop=FALSE], 2, function(varVn) mean(varVn, na.rm=TRUE))

        if(sum(samVl) == 1)
            varDF[, "sample_sd"] <- rep(0, nrow(varDF))
        else
            varDF[, "sample_sd"] <- apply(datMN[samVl, , drop=FALSE], 2, function(varVn) sd(varVn, na.rm=TRUE))

        varDF[, "sample_CV"] <- varDF[, "sample_sd"] / varDF[, "sample_mean"]

    }

    ## 'blank' mean / 'sample' mean ratio

    if(all(c("blank_mean", "sample_mean") %in% colnames(varDF))) {

        cat("\nVariables: Blank mean over sample mean\n", sep="")

        varDF[, "blankMean_over_sampleMean"] <- varDF[, "blank_mean"] / varDF[, "sample_mean"]

    }

    ## 'pool' observations

    if("sampleType" %in% colnames(samDF) && "pool" %in% samDF[, "sampleType"]) {

        cat("\nVariables: Pool mean, sd, and CV\n", sep="")

        pooVl <- samDF[, "sampleType"] == "pool"

        if(sum(pooVl) == 1)
            varDF[, "pool_mean"] <- datMN[pooVl, ]
        else
            varDF[, "pool_mean"] <- apply(datMN[pooVl, , drop=FALSE], 2, function(varVn) mean(varVn, na.rm=TRUE))

        if(sum(pooVl) == 1)
            varDF[, "pool_sd"] <- rep(0, nrow(varDF))
        else
            varDF[, "pool_sd"] <- apply(datMN[pooVl, , drop=FALSE], 2, function(varVn) sd(varVn, na.rm=TRUE))

        varDF[, "pool_CV"] <- varDF[, "pool_sd"] / varDF[, "pool_mean"]

    }

    ## 'pool' CV / 'sample' CV ratio

    if(all(c("pool_CV", "sample_CV") %in% colnames(varDF))) {

        cat("\nVariables: Pool CV over sample CV\n", sep="")

        varDF[, "poolCV_over_sampleCV"] <- varDF[, "pool_CV"] / varDF[, "sample_CV"]

    }


    ## 'pool' dilutions

    if("sampleType" %in% colnames(samDF) && any(grepl("pool.+", samDF[, "sampleType"]))) {

        pooVi <- grep("pool.*", samDF[, "sampleType"]) ## pool, pool2, pool4, poolInter, ...

        pooNamVc <- samDF[pooVi, "sampleType"]

        if(pooAsPo1L) {

            pooNamVc[pooNamVc == "pool"] <- "pool1" ## 'pool' -> 'pool1'

        } else {

            pooVl <- pooNamVc == "pool"
            pooVi <- pooVi[!pooVl]
            pooNamVc <- pooNamVc[!pooVl]

        }

        pooDilVc <- gsub("pool", "", pooNamVc)

        pooDilVl <- sapply(pooDilVc, allDigF)

        if(sum(pooDilVl)) {

            cat("\nVariables: Pool dilutions\n", sep="")

            pooNamVc <- pooNamVc[pooDilVl] ## for the plot

            pooVi <- pooVi[pooDilVl]

            dilVn <- 1 / as.numeric(pooDilVc[pooDilVl])

            varDF[, "poolDil_correl"] <- apply(datMN[pooVi, , drop=FALSE], 2,
                                               function(varVn) cor(dilVn, varVn))

            varDF[, "poolDil_pval"] <- apply(datMN[pooVi, , drop=FALSE], 2,
                                             function(varVn) cor.test(dilVn, varVn)[["p.value"]])

        }

    }


    ##------------------------------
    ## Sample metrics
    ##------------------------------


    ## Hotelling: p-value associated to the distance from the center in the first PCA score plane

    cat("\nObservations: Hotelling ellipse\n", sep="")

    ropLs <- opls(datMN, predI = 2, plotL = FALSE, printL = FALSE)

    ropScoreMN <- getScoreMN(ropLs)

    invCovScoMN <- solve(cov(ropScoreMN))

    n <- nrow(datMN)
    hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))

    hotPvaVn <- apply(ropScoreMN,
                      1,
                      function(x)
                      1 - pf(1 / hotN * t(as.matrix(x)) %*% invCovScoMN %*% as.matrix(x), 2, n - 2))

    samDF[, "hotelling_pval"] <- hotPvaVn

    ## p-value associated to number of missing values

    cat("\nObservations: Missing values\n", sep="")

    missZscoVn <- zScoreF(apply(datMN,
                                1,
                                function(rowVn) {
                                    sum(is.na(rowVn))
                                }))

    samDF[, "missing_pval"] <- sapply(missZscoVn, function(zscoN) 2 * (1 - pnorm(abs(zscoN))))

    ## p-value associated to the deciles of the profiles

    cat("\nObservations: Profile deciles\n", sep="")

    deciMN <- t(as.matrix(apply(datMN,
                                1,
                                function(x) quantile(x, 0.1 * 1:9, na.rm = TRUE))))

    deciZscoMN <- apply(deciMN, 2, zScoreF)

    deciZscoMaxVn <- apply(deciZscoMN, 1, function(rowVn) rowVn[which.max(abs(rowVn))])

    samDF[, "decile_pval"] <- sapply(deciZscoMaxVn, function(zscoN) 2 * (1 - pnorm(abs(zscoN))))


    ##------------------------------
    ## Figure
    ##------------------------------

    cat("\nPlotting\n")

    if(!is.null(fig.pdfC)) {
        pdf(fig.pdfC, width=11, height=7)
    } else
        dev.new(width=11, height=7)

    datPloF()

    if(!is.null(fig.pdfC))
        dev.off()


    ##------------------------------
    ## End
    ##------------------------------


    if(!is.null(log.txtC)) {

        cat("\nEnd of '", modNamC, "' Galaxy module call: ",
            as.character(Sys.time()), "\n", sep = "")

        cat("\n\n\n============================================================================")
        cat("\nAdditional information about the call:\n")
        cat("\n1) Parameters:\n")
        print(args)

        cat("\n2) Session Info:\n")
        sessioninfo <- sessionInfo()
        cat(sessioninfo$R.version$version.string,"\n")
        cat("Main packages:\n")
        for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
        cat("Other loaded packages:\n")
        for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")        

        cat("============================================================================\n")
        
        sink()
    }

    options(stingsAsFactors = strAsFacL)
    options(warn = optWrnN)

    return(list(samDF=samDF,
                varDF=varDF))


} ## qualityMetricsF
