library(xcms)
library(faahKO)

#BiocManager::install("faahKO")

parallel <- TRUE

cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
        recursive = TRUE)

cdfs <- cdfs[c(1, 2, 7, 8)]

pd <- data.frame(
        sample_name <- sub(basename(cdfs), pattern = ".CDF",
        replacement = "", fixed = TRUE),
        sample_group = c(rep("KO", 2), rep("WT", 2)),
        stringsAsFactors = FALSE)

cwp <- CentWaveParam()

if (parallel) {
    for (cdf in cdfs) {
        raw_data <- readMSData(files = cdf, pdata = new("NAnnotatedDataFrame", pd),
                    mode = "onDisk")
        xdata <- findChromPeaks(raw_data, param = cwp)
        if (!exists("xdata_merged"))
            xdata_merged <- xdata
        else
            xdata_merged <- c(xdata_merged, xdata)
    }
    xdata <- xdata_merged
    if (!is.null(args$sampleMetadata)) {
        sampleMetadataFile <- args$sampleMetadata
        sampleMetadata <- getDataFrameFromFile(sampleMetadataFile, header = F)
        xdata@phenoData@data$sample_group <- sampleMetadata$V2[match(xdata@phenoData@data$sample_name, sampleMetadata$V1)]

        if (any(is.na(pData(xdata)$sample_group))) {
            sample_missing <- pData(xdata)$sample_name[is.na(pData(xdata)$sample_group)]
            error_message <- paste("Those samples are missing in your sampleMetadata:", paste(sample_missing, collapse = " "))
            print(error_message)
            stop(error_message)
        }
    }
} else {
    raw_data <- readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd),
                mode = "onDisk")
    xdata <- findChromPeaks(raw_data, param = cwp)
}

pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
        bw = 5,
        minFraction = 0.3,
        minSamples = 1,
        binSize = 0.01,
        maxFeatures = 50)
xdata <- groupChromPeaks(xdata, param = pdp)

# WARNING: not align with the tests parameters
pgp <- PeakGroupsParam(minFraction = 0.85)
xdata <- adjustRtime(xdata, param = pgp)

pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
        minFraction = 0.4, bw = 20)
xdata <- groupChromPeaks(xdata, param = pdp)

xdata <- fillChromPeaks(xdata)

# This function retrieve a xset like object
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return(xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, "xcmsSet"))
        if (is.null(xset@phenoData$sample_group))
            sampclass(xset) <- "."
        else
            sampclass(xset) <- xset@phenoData$sample_group
        return(xset)
    }
}

xset <- getxcmsSetObject(xdata)


library(CAMERA)
xa <- annotate(
        xset,
        nSlaves = 1,
        sigma = 6,
        perfwhm = 0.6,
        ppm = 5,
        mzabs = 0.015,
        maxcharge = 3,
        maxiso = 4,
        minfrac = 0.5,
        quick = TRUE,
        intval = "into")

diffrep <- diffreport(
        object = xset,
        class1 = "KO",
        class2 = "WT",
        filebase = "KO-vs-WT",
        eicmax = 200,
        eicwidth = 200,
        sortpval = TRUE,
        value = "into",
        h = 480,
        w = 640,
        mzdec = 2,
        missing = 0)
