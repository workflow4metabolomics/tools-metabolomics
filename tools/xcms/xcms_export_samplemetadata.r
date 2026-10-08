#!/usr/bin/env Rscript
# Authors:
#   - ABiMS Team
#   - LABERCA - PARC project founding

# Import the different functions
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("xcms", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n")

args <- parseCommandArgs(evaluate = FALSE)

sampleMetadata <- NULL
for (image in args$images) {
    load(image)
    if (exists("raw_data")) xdata <- raw_data
    if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS")
    if (is.null(sampleMetadata)) {
        sampleMetadata <- xdata@sampleData
    } else {
        sampleMetadata <- rbind(sampleMetadata, xdata@sampleData)
    }
}
sampleMetadata <- sampleMetadata[, c("sample_name", "sample_group")]
print(sampleMetadata)
colnames(sampleMetadata) <- c("sample_name", "sample_group")
sampleMetadata$sample_name <- make.names(sampleMetadata$sample_name)

# Create a sampleMetada file
write.table(sampleMetadata, file = "sampleMetadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
