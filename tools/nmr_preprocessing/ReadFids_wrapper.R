## ------------------------------
## Libraries laoding
## ------------------------------
# library(batch)
library(optparse) # argument parsing


## Outputs
dataMatrix <- "dataMatrix.tsv"
sampleMetadata <- "sampleMetadata.tsv"
nomGraphe <- "graphOut.pdf"
logOut <- "logOut.txt"
sink(logOut, append = TRUE, split = TRUE)

# ------------------------------
# Command line interface (optparse)
# ------------------------------
option_list <- list(
    make_option(c("-f", "--fidzipfile"), type = "character", help = "Path to zipped Bruker FID file"),
    make_option(c("-t", "--title_line"), type = "character", help = "Title line for output"),
    make_option(c("-s", "--subdirectories"), action = "store_true", dest = "subdirectories", default = FALSE, help = "Whether to use subdirectories (boolean flag)"),
    make_option(c("-d", "--dirs_names"), action = "store_true", dest = "dirs_names", default = FALSE, help = "Whether to use dirs_names (boolean flag)")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)
print(args)


## ======================================================
## ======================================================
## Parameters Loading
## ======================================================
## ======================================================

## Inputs
# Path
## Bruker FIDs
fileType <- "Bruker"
zipfile <- args$fidzipfile
directory <- unzip(zipfile, list = FALSE)
path <- paste(getwd(), strsplit(directory[1], "/")[[1]][2], sep = "/")
print(path)

# other inputs from ReadFids
l <- args$title_line
subdirs <- args$subdirectories
dirs.names <- args$dirs_names


## ======================================================
## ======================================================
## Computation
## ======================================================
## ======================================================
if (length(warnings()) > 0) { # or !is.null(warnings())
    print("something happened")
}

## Starting
cat("\nStart of 'ReadFids' Galaxy module call: ", as.character(Sys.time()), "\n\n", sep = "")

outputs <- ReadFids(path = path, l = l, subdirs = subdirs, dirs.names = dirs.names)
data_matrix <- outputs[["Fid_data"]] # Data matrix
data_sample <- outputs[["Fid_info"]] # Sample metadata

pdf(nomGraphe, onefile = TRUE, width = 13, height = 13)
DrawSignal(data_matrix,
    subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
    xlab = "Frequency", num.stacked = 4,
    main = "Raw FID data", createWindow = FALSE
)
invisible(dev.off())

## ======================================================
## ======================================================
## Saving
## ======================================================
## ======================================================

# Data matrix
write.table(data_matrix, file = dataMatrix, quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# Sample metadata
write.table(data_sample, file = sampleMetadata, quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")

# Reproductibility
print(sessionInfo())
## Ending

cat("\nEnd of 'ReadFids' Galaxy module call: ", as.character(Sys.time()), sep = "")

sink()
