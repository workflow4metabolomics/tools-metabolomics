## Libraries laoding
## -----------------
library(ptw)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)

# In-house function for argument parsing
parse_args <- function() {
  args <- commandArgs()
  start <- which(args == "--args")[1] + 1
  if (is.na(start)) {
    return(list())
  }
  seq_by2 <- seq(start, length(args), by = 2)
  result <- as.list(args[seq_by2 + 1])
  names(result) <- args[seq_by2]
  return(result) # nolint: return_linter.
}

# R script call
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep = "/"))
}
# Import the different functions
source_local("NmrPreprocessing_script.R")
source_local("DrawFunctions.R")

## ------------------------------
## Script
## ------------------------------
run_example_l <- FALSE

if (!run_example_l) {
  arg_ls <- unlist(parse_args())
}
# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")
print(arg_ls)

## ------------------------------
## Constants
## ------------------------------
top_env_c <- environment()
flag_c <- "\n"

## Starting
cat("\nStart of 'Preprocessing' Galaxy module call: ",
    as.character(Sys.time()), "\n", sep = "")

## ======================================================
## ======================================================
## Parameters Loading
## ======================================================
## ======================================================

# graphical inputs
first_opc_graph <- arg_ls[["first_opc_graph"]]
ss_graph <- arg_ls[["ss_graph"]]
apod_graph <- arg_ls[["apod_graph"]]
ft_graph <- arg_ls[["ft_graph"]]
sr_graph <- arg_ls[["sr_graph"]]
zero_opc_graph <- arg_ls[["zero_opc_graph"]]
bc_graph <- arg_ls[["bc_graph"]]
final_graph <- arg_ls[["final_graph"]]

# 1rst order phase correction ------------------------
# Inputs
## Data matrix
fid_data0 <- read.table(arg_ls[["data_matrix_fid"]], header = TRUE,
                        check.names = FALSE, sep = "\t")
fid_data0 <- as.matrix(fid_data0)

## Samplemetadata
samplemetadata_fid <- read.table(arg_ls[["sample_metadata_fid"]],
                                 check.names = FALSE, header = TRUE, sep = "\t")
samplemetadata_fid <- as.matrix(samplemetadata_fid)

# water and solvent(s) correction ------------------------
# Inputs
lambda <- as.numeric(arg_ls[["lambda"]])

# apodization -----------------------------------------
# Inputs
phase <- 0
rect_ratio <- 1 / 2
gauss_lb <- 1
exp_lb <- 1
apodization <- arg_ls[["apodization_method"]]

if (apodization == "exp") {
  exp_lb <- as.numeric(arg_ls[["exp_lb"]])
} else if (apodization == "cos2") {
  phase <- as.numeric(arg_ls[["phase"]])
} else if (apodization == "hanning") {
  phase <- as.numeric(arg_ls[["phase"]])
} else if (apodization == "hamming") {
  phase <- as.numeric(arg_ls[["phase"]])
} else if (apodization == "blockexp") {
  rect_ratio <- as.numeric(arg_ls[["rect_ratio"]])
  exp_lb <- as.numeric(arg_ls[["exp_lb"]])
} else if (apodization == "blockcos2") {
  rect_ratio <- as.numeric(arg_ls[["rect_ratio"]])
} else if (apodization == "gauss") {
  rect_ratio <- as.numeric(arg_ls[["rect_ratio"]])
  gauss_lb <- as.numeric(arg_ls[["gauss_lb"]])
}

# Fourier transform ----------------------------------
# Inputs

# Zero Order Phase Correction -------------------------------
# Inputs
angle <- NULL
exclude_zopc <- NULL

zero_order_phase_method <- arg_ls[["zero_order_phase_method"]]

if (zero_order_phase_method == "manual") {
  angle <- arg_ls[["angle"]]
}

exclude_zone_zerophase <- arg_ls[["exclude_zone_zerophase.choice"]]
if (exclude_zone_zerophase == "YES") {
  exclude_zone_zerophase_list <- list()
  for (i in which(names(arg_ls) == "exclude_zone_zerophase_left")) {
    exclude_zone_zerophase_left <- as.numeric(arg_ls[[i]])
    exclude_zone_zerophase_right <- as.numeric(arg_ls[[i + 1]])
    exclude_zone_zerophase_list <- c(exclude_zone_zerophase_list,
                                     list(c(exclude_zone_zerophase_left,
                                            exclude_zone_zerophase_right)))
  }
  exclude_zopc <- exclude_zone_zerophase_list
}

# Internal referencering ----------------------------------
# Inputs
shift_threshold <- 2
ppm <- TRUE
shift_referencing_range_list <- NULL # fromto_rc
pct_near_value <- 0.02 # pc
rowindex_graph <- NULL
ppm_ref <- 0 # ppm.ref

shift_referencing_range <- arg_ls[["shift_referencing_range"]]
if (shift_referencing_range == "nearvalue") {
  pct_near_value <- as.numeric(arg_ls[["pctNearValue"]])
}

if (shift_referencing_range == "window") {
  shift_referencing_range_list <- list()
  for (i in which(names(arg_ls) == "shift_referencing_range_left")) {
    shift_referencing_range_left <- as.numeric(arg_ls[[i]])
    shift_referencing_range_right <- as.numeric(arg_ls[[i + 1]])
    shift_referencing_range_list <- c(shift_referencing_range_list,
                                      list(c(shift_referencing_range_left,
                                             shift_referencing_range_right)))
  }
}
shift_handling <- arg_ls[["shift_handling"]]

ppmvalue <- as.numeric(arg_ls[["ppmvalue"]])
# }

# Baseline Correction -------------------------------
# Inputs
lambda_bc <- as.numeric(arg_ls[["lambda_bc"]])
p_bc <- as.numeric(arg_ls[["p_bc"]])
epsilon <- as.numeric(arg_ls[["epsilon"]])

exclude_bc <- NULL

exclude_zone_bc <- arg_ls[["exclude_zone_bc.choice"]]
if (exclude_zone_bc == "YES") {
  exclude_zone_bc_list <- list()
  for (i in which(names(arg_ls) == "exclude_zone_bc_left")) {
    exclude_zone_bc_left <- as.numeric(arg_ls[[i]])
    exclude_zone_bc_right <- as.numeric(arg_ls[[i + 1]])
    exclude_zone_bc_list <- c(exclude_zone_bc_list,
                              list(c(exclude_zone_bc_left,
                                     exclude_zone_bc_right)))
  }
  exclude_bc <- exclude_zone_bc_list
}

# transformation of negative values -------------------------------
# Inputs
negative_to_zero <- arg_ls[["negative_to_zero"]]

# Outputs
nom_graphe <- arg_ls[["graphOut"]]
log <- arg_ls[["logOut"]]

## Checking arguments
## -------------------
error_stock <- "\n"
if (length(error_stock) > 1) {
  stop(error_stock)
}

## ======================================================
## Computation
## ======================================================
pdf(nom_graphe, onefile = TRUE, width = 13, height = 13)

# FirstOrderPhaseCorrection ---------------------------------
print("group_delay_correction")
fid_data <- group_delay_correction(fid_data0, fid_info = samplemetadata_fid,
                                   group_delay = NULL)

if (first_opc_graph == "YES") {
  title <- "FIDs after Group Delay Correction"
  DrawSignal(fid_data, subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
    xlab = "Frequency", num.stacked = 4, main = title,
    createWindow = FALSE
  )
}

# SolventSuppression ---------------------------------
print("solvent_suppression")
fid_data <- solvent_suppression(fid_data, lambda_ss = lambda,
                                ptw_ss = TRUE, plot_solvent = FALSE,
                                return_solvent = FALSE)

if (ss_graph == "YES") {
  title <- "FIDs after Solvent Suppression "
  DrawSignal(fid_data, subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}

print(rect_ratio)
# Apodization ---------------------------------
print("apodization")
fid_data <- Apodization(fid_data, fid_info = samplemetadata_fid, dt = NULL,
  type_apod = apodization, phase = phase, rect_ratio = rect_ratio,
  gauss_lb = gauss_lb, exp_lb = exp_lb, plot_window = FALSE,
  return_factor = FALSE
)

if (apod_graph == "YES") {
  title <- "FIDs after Apodization"
  DrawSignal(fid_data, subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
    xlab = "Frequency", num.stacked = 4,
    main = title, createWindow = FALSE
  )
}

# FourierTransform ---------------------------------
print("fourier_transform")
spectrum_data <- fourier_transform(fid_data, fid_info = samplemetadata_fid,
                                   reverse_axis = TRUE)

if (ft_graph == "YES") {
  title <- "Fourier transformed spectra"
  DrawSignal(spectrum_data, subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
    xlab = "Frequency", num.stacked = 4, main = title, createWindow = FALSE
  )
}

# ZeroOrderPhaseCorrection ---------------------------------
print("zero_order_phase_correction")
spectrum_data <- ZeroOrderPhaseCorrection(spectrum_data,
  type_zopc = zero_order_phase_method, plot_rms = NULL, return_angle = FALSE,
  createWindow = TRUE, angle = angle, plot_spectra = FALSE,
  ppm.zopc = TRUE, exclude.zopc = exclude_zopc
)


# InternalReferencing ---------------------------------
# if (shiftReferencing=="YES") {
spectrum_data <- internal_referencing(spectrum_data, samplemetadata_fid,
  method = "max", range = shift_referencing_range,
  ppm_value = ppmvalue, shift_handling = shift_handling, ppm_ir = TRUE,
  fromto_rc = shift_referencing_range_list, pc = pct_near_value
)

if (sr_graph == "YES") {
  title <- "Spectra after Shift Referencing"
  DrawSignal(spectrum_data, subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
    xlab = "Frequency", num.stacked = 4, main = title, createWindow = FALSE
  )
}

if (zero_opc_graph == "YES") {
  title <- "Spectra after Zero Order Phase Correction"
  DrawSignal(spectrum_data, subtype = "stacked",
    ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
    xlab = "Frequency", num.stacked = 4, main = title, createWindow = FALSE
  )
}


# BaselineCorrection ---------------------------------
print("baseline_correction")
spectrum_data <- BaselineCorrection(spectrum_data, ptw.bc = TRUE,
  lambda.bc = lambda_bc, p.bc = p_bc, eps = epsilon, ppm.bc = TRUE,
  exclude.bc = exclude_bc, returnBaseline = FALSE
)

if (bc_graph == "YES") {
  title <- "Spectra after Baseline Correction"
  DrawSignal(spectrum_data, subtype = "stacked",
             ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
             xlab = "Frequency", num.stacked = 4, main = title,
             createWindow = FALSE)
}


# negative_values_zeroing ---------------------------------
print("negative_values_zeroing")
if (negative_to_zero == "YES") {
  spectrum_data <- negative_values_zeroing(spectrum_data)
}
print(spectrum_data[1:5, 1:5])
if (final_graph == "YES") {
  title <- "Final preprocessed spectra"
  DrawSignal(spectrum_data, subtype = "stacked",
             ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = TRUE,
             xlab = "Frequency", num.stacked = 4, main = title,
             createWindow = FALSE)
}
invisible(dev.off())

# Variable metadata creation
data_variable <- matrix(NA, nrow = 1, ncol = dim(spectrum_data)[2],
                        dimnames = list("ID", NULL))
colnames(data_variable) <- colnames(spectrum_data)
data_variable[1, ] <- colnames(data_variable)

## ======================================================
## Saving
## ======================================================

# Data Matrix
write.table(round(t(Re(spectrum_data)), 6), file = arg_ls[["dataMatrix"]],
            quote = FALSE, row.names = TRUE, sep = "\t", col.names = TRUE)

# Variable metadata
write.table(data_variable, file = arg_ls[["variableMetadata"]], quote = FALSE,
            row.names = TRUE, sep = "\t", col.names = TRUE)

# input arguments
cat("\n INPUT and OUTPUT ARGUMENTS :\n")
arg_ls

## Ending
cat("\nVersion of R librairies")
print(sessionInfo())
cat("\nEnd of 'Preprocessing' Galaxy module call: ",
    as.character(Sys.time()), sep = "")

rm(list = ls())