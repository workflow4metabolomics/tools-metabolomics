## ------------------------------
## Libraries laoding
## ------------------------------
library(optparse) # argument parsing

## Outputs
datamatrix <- "datamatrix.tsv"
samplemetadata <- "samplemetadata.tsv"
nom_graphe <- "graphOut.pdf"
log_out <- "log_out.txt"
sink(log_out, append = TRUE, split = TRUE)

# ------------------------------
# Command line interface (optparse)
# ------------------------------
option_list <- list(
  make_option(c("-f", "--fidzipfile"),
    type = "character",
    help = "Path to zipped Bruker FID file"
  ),
  make_option(c("-t", "--title_line"),
    type = "character",
    help = "Title line for output"
  ),
  make_option(c("-s", "--subdirectories"),
    action = "store_true",
    dest = "subdirectories", default = FALSE,
    help = "Whether to use subdirectories (boolean flag)"
  ),
  make_option(c("-d", "--dirs_names"),
    action = "store_true",
    dest = "dirs_names", default = FALSE,
    help = "Whether to use dirs_names (boolean flag)"
  )
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)
print(args)

## ======================================================
## Parameters Loading
## ======================================================

## Inputs
# Path
## Bruker FIDs
file_type <- "Bruker"
zipfile <- args$fidzipfile
directory <- unzip(zipfile, list = FALSE)
path <- paste(getwd(), strsplit(directory[1], "/")[[1]][2], sep = "/")
print(path)

# other inputs from ReadFids
l <- args$title_line
subdirs <- args$subdirectories
dirs_names <- args$dirs_names

## ======================================================
## Computation
## ======================================================
if (length(warnings()) > 0) { # or !is.null(warnings())
  print("something happened")
}

## Starting
cat("\nStart of 'ReadFids' Galaxy module call: ",
  as.character(Sys.time()), "\n\n",
  sep = ""
)

outputs <- read_fids(
  path = path, l = l, subdirs = subdirs,
  dirs_names = dirs_names
)
data_matrix <- outputs[["fid_data"]] # Data matrix
data_sample <- outputs[["fid_info"]] # Sample metadata

pdf(nom_graphe, onefile = TRUE, width = 13, height = 13)
draw_signal(data_matrix,
  subtype = "stacked",
  re_im_mod_arg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
  xlab = "Frequency", num_stacked = 4, main = "Raw FID data",
  create_window_drawsignal = FALSE
)
invisible(dev.off())

cat("Data matrix dimensions: ")
print(dim(data_matrix))
cat("Sample metadata dimensions: ")
print(dim(data_sample))

## ======================================================
## Saving
## ======================================================

# Data matrix
write.table(data_matrix,
  file = datamatrix, quote = FALSE, row.names = TRUE,
  sep = "\t", col.names = TRUE
)

# Sample metadata
write.table(data_sample,
  file = samplemetadata, quote = FALSE, row.names = TRUE,
  sep = "\t", col.names = TRUE
)

# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")

# Reproductibility
print(sessionInfo())
## Ending

cat("\nEnd of 'ReadFids' Galaxy module call: ",
    as.character(Sys.time()), sep = "")

sink()
