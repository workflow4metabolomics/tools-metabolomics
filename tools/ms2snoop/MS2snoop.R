#'
#' read and process mspurity W4M files
#' create a summary of fragment for each precursor and a graphics of peseudo
#' spectra + correlation on which checking of fragment is based on
#' V3 try to identify and process multiple files for 1 precursor which may
#' occur if different collision energy are used
#' V4 elimination of correlation = NA. Correlation is done with precursor, if
#' precursor is not present correlation with most intense peak
#' author: Jean-Francois Martin
#' V5 is versionned, lintR-compliant, packaged, unit-tested, documented and
#' tested against data from other labs.
#' new maintainer: Lain Pavot - lain.pavot@inrae.fr
#'
#' @import optparse
#'


assign("MS2SNOOP_VERSION", "2.0.0")
lockBinding("MS2SNOOP_VERSION", globalenv())

assign("MISSING_PARAMETER_ERROR", 1)
lockBinding("MISSING_PARAMETER_ERROR", globalenv())

assign("BAD_PARAMETER_VALUE_ERROR", 2)
lockBinding("BAD_PARAMETER_VALUE_ERROR", globalenv())

assign("MISSING_INPUT_FILE_ERROR", 3)
lockBinding("MISSING_INPUT_FILE_ERROR", globalenv())

assign("NO_ANY_RESULT_ERROR", 255)
lockBinding("NO_ANY_RESULT_ERROR", globalenv())

assign("DEFAULT_PRECURSOR_PATH", "peaklist_precursors.tsv")
assign("DEFAULT_FRAGMENTS_PATH", "peaklist_fragments.tsv")
assign("DEFAULT_COMPOUNDS_PATH", "compounds_pos.txt")
assign("DEFAULT_OUTPUT_PATH", "compound_fragments_result.txt")
assign("DEFAULT_TOLMZ", 0.01)
assign("DEFAULT_TOLRT", 20)
assign("DEFAULT_MZDECIMAL", 3)
assign("DEFAULT_R_THRESHOLD", 0.85)
assign("DEFAULT_MINNUMBERSCAN", 8)
assign("DEFAULT_SEUIL_RA", 0.05)
lockBinding("DEFAULT_PRECURSOR_PATH", globalenv())
lockBinding("DEFAULT_FRAGMENTS_PATH", globalenv())
lockBinding("DEFAULT_COMPOUNDS_PATH", globalenv())
lockBinding("DEFAULT_OUTPUT_PATH", globalenv())
lockBinding("DEFAULT_TOLMZ", globalenv())
lockBinding("DEFAULT_TOLRT", globalenv())
lockBinding("DEFAULT_MZDECIMAL", globalenv())
lockBinding("DEFAULT_R_THRESHOLD", globalenv())
lockBinding("DEFAULT_MINNUMBERSCAN", globalenv())
lockBinding("DEFAULT_SEUIL_RA", globalenv())

assign("DEFAULT_EXTRACT_FRAGMENTS_R_THRESHOLD", 0.85)
assign("DEFAULT_EXTRACT_FRAGMENTS_SEUIL_RA", 0.1)
assign("DEFAULT_EXTRACT_FRAGMENTS_TOLMZ", 0.01)
assign("DEFAULT_EXTRACT_FRAGMENTS_TOLRT", 60)
lockBinding("DEFAULT_EXTRACT_FRAGMENTS_R_THRESHOLD", globalenv())
lockBinding("DEFAULT_EXTRACT_FRAGMENTS_SEUIL_RA", globalenv())
lockBinding("DEFAULT_EXTRACT_FRAGMENTS_TOLMZ", globalenv())
lockBinding("DEFAULT_EXTRACT_FRAGMENTS_TOLRT", globalenv())


########################################################################


get_formulas <- function(
  mzref,
  spectra,
  nominal_mz_list,
  processing_parameters,
  background = !TRUE
) {
  if (is.vector(mzref) && length(mzref) > 1) {
    return(lapply(
      mzref,
      function(mz) {
        return(get_formulas(
          mzref = mz,
          spectra = spectra,
          nominal_mz_list = nominal_mz_list,
          processing_parameters = processing_parameters,
          background = background
        ))
      }
    ))
  }
  input <- sprintf(
    "%s-%s.ms",
    gsub("[[:space:]]", "_", processing_parameters$c_name),
    mzref
  )
  output <- sprintf(
    "out/%s-%s.out",
    gsub("[[:space:]]", "_", processing_parameters$c_name),
    mzref
  )
  file_content <- paste(
    sprintf(">compound %s", processing_parameters$c_name),
    sprintf(">ionization %s", processing_parameters$ionization),
    sprintf(">parentmass %s", mzref),
    sprintf(">formula %s", processing_parameters$elemcomposition),
    sep = "\n"
  )
  displayed_file_content <- sprintf(
    "%s\n>collision\n%s",
    file_content,
    paste(
      sprintf(
        "%s %s",
        spectra[1:3, "mz"],
        spectra[1:3, "intensities"]
      ),
      collapse = "\n"
    )
  )
  if (nrow(spectra) > 3) {
    displayed_file_content <- sprintf(
      "%s\n[... %s more rows of mz and intensities ...]",
      displayed_file_content,
      nrow(spectra) - 3
    )
  }
  catf(
    ">> MS file created for %s with content:\n%s\n",
    processing_parameters$c_name,
    displayed_file_content
  )
  file_content <- sprintf(
    "%s\n\n>collision\n%s",
    file_content,
    paste(
      sprintf("%s %s", spectra$mz, spectra$intensities),
      collapse = "\n"
    )
  )
  cat(file_content, file = input, append = FALSE)
  command <- sprintf(
    paste(
      "sirius",
      "--noCite",
      "--noSummaries",
      "--loglevel=WARNING",
      "-i='%s'",
      "-o='%s'",
      "tree 2>&1 | grep -P '^(WARNING|SEVERE)'"
    ),
    input,
    output
  )
  verbose_catf(
    ">> Sirius is running %swith the command: %s\n",
    if (background) "in the background " else "",
    command
  )
  system(
    command,
    wait = !background,
    ignore.stdout = background,
    ignore.stderr = background
  )
  return(extract_sirius_results(output, spectra[, "mz"]))
  # return(extract_sirius_results(output, nominal_mz_list))
}

extract_sirius_results <- function(
  output,
  nominal_mz_list,
  mz_tolerance = 2,
  max_error = 5,
  min_error = 0.01
) {

  if (mz_tolerance < 0) {
    stop_with_status(
      "mz tolerance cannot be negative",
      BAD_PARAMETER_VALUE_ERROR
    )
  }

  if (mz_tolerance == 0) {
    stop_with_status(
      "mz tolerance cannot be zero",
      BAD_PARAMETER_VALUE_ERROR
    )
  }

  if (max_error < 0) {
    stop_with_status(
      "max mz error cannot be negative",
      BAD_PARAMETER_VALUE_ERROR
    )
  }

  if (max_error == 0) {
    stop_with_status(
      "max mz error cannot be zero",
      BAD_PARAMETER_VALUE_ERROR
    )
  }

  if (min_error < 0) {
    stop_with_status(
      "min mz error cannot be negative",
      BAD_PARAMETER_VALUE_ERROR
    )
  }

  if (min_error == 0) {
    stop_with_status(
      "min mz error cannot be zero",
      BAD_PARAMETER_VALUE_ERROR
    )
  }
  if (max_error < mz_tolerance) {
    stop_with_status(
      "max mz error should not be lesser than mz tolerance",
      BAD_PARAMETER_VALUE_ERROR
    )
  }
  if (max_error <= min_error) {
    stop_with_status(
      "max mz error should not be greater than min error",
      BAD_PARAMETER_VALUE_ERROR
    )
  }

  out_dir <- sprintf(
    "%s/%s",
    list.dirs(output, recursive = FALSE)[[1]],
    "spectra"
  )
  result_filename <- sprintf(
    "%s/%s",
    out_dir,
    list.files(out_dir)[[1]]
  )

  if (!is.null(result_filename)) {
    formula_df <- get_csv_or_tsv(result_filename)
  } else {
    return(rep(NA, length(nominal_mz_list)))
  }

  formulas <- data.frame(formula = c(), mz = c(), error = c())
  for (nominal_mz in nominal_mz_list) {
    curent_mz_tolerance <- mz_tolerance

    ## just a bad value, to enter the while at least once
    match <- data.frame()

    while (nrow(match) != 1) {
      match <- formula_df[
        which(abs(nominal_mz - formula_df$mz) <= curent_mz_tolerance),
      ]
      if (nrow(match) == 0) {
        ## if no match, increase the tolerance
        curent_mz_tolerance <- curent_mz_tolerance * 1.5
      } else {
        ## decrease the tolerance to filter more
        curent_mz_tolerance <- curent_mz_tolerance / 2.0
      }
      if (curent_mz_tolerance <= min_error) {
        warning(
          paste(
            "The mz tolerance to try to match fragment with their formula",
            "has lowered too low! This means that the sirius spectra output",
            "contains multiple fragments with the same mz, and we cannot deal",
            "with that."
          )
        )
        match <- data.frame(
          mz = mean(match$mz),
          formula = paste(match$formula, collapse = ";")
        )
      }
    }
    match$error <- round(nominal_mz - match$mz, 4)
    if (abs(match$error) > max_error) {
      catf(
        paste(
          "Could not find formula for fragment %s.",
          "Error is greater than %s (%s with %s)\n"
        ),
        nominal_mz, max_error, abs(match$error), match$formula
      )
      match$formula <- NA
      match$error <- NA
    } else {
      catf(
        "Fragment with mz=%s matches %s with an error of %s\n",
        nominal_mz, match$formula, match$error
      )
    }
    curent_row <- nrow(formulas) + 1
    formulas[curent_row, "formula"] <- match$formula
    formulas[curent_row, "mz"] <- match$mz
    formulas[curent_row, "error"] <- match$error
  }
  return(formulas)
}

#' @title plot_pseudo_spectra
#' @param x
#' @param fid
#' @param sum_int
#' @param vmz
#' @param cor_abs_int
#' @param refcol
#' @param c_name
#' @description plot_pseudo_spectra
#' function to compute sum of intensities among scans for all
#' m/z kept (cor > r_threshold & minimum number of scans)
#' and plot pseudo spectra
#' x dataframe scan X fragments with scans number in the 1st column and
#' ions in next with intensities
#' fid file id when several a precursor has been detected in several files
plot_pseudo_spectra <- function(
  x,
  fid,
  sum_int,
  vmz,
  cor_abs_int,
  refcol,
  meaned_mz,
  processing_parameters
) {
  ## du fait de la difference de nombre de colonne entre la dataframe qui
  ## inclue les scans en 1ere col, mzRef se decale de 1
  refcol <- refcol - 1
  ## compute relative intensities max=100%
  rel_int <- sum_int[-1]
  rel_int <- rel_int / max(rel_int)

  if (processing_parameters$do_pdf) {
    ## define max value on vertical axis (need to increase in order to plot the
    ## label of fragments)
    ymax <- max(rel_int) + 0.2 * max(rel_int)

    par(mfrow = c(2, 1))
    plot(vmz, rel_int, type = "h", ylim = c(0, ymax),
      main = processing_parameters$c_name
    )
    ## low correl coef. will be display in grey
    cor_low <- which(round(cor_abs_int, 2) < processing_parameters$r_threshold)

    lbmzcor <- sprintf("%s(r=%s)", vmz, round(cor_abs_int, 2))

    if (length(cor_low) > 0) {
      text(
        vmz[cor_low],
        rel_int[cor_low],
        lbmzcor[cor_low],
        cex = 0.5,
        col = "grey",
        srt = 90,
        adj = 0
      )
      if (length(vmz) - length(cor_low) > 1) {
        text(
          vmz[-c(refcol, cor_low)],
          rel_int[-c(refcol, cor_low)],
          lbmzcor[-c(refcol, cor_low)],
          cex = 0.6,
          col = 1,
          srt = 90,
          adj = 0
        )
      }
    } else {
      if (length(vmz) > 1) {
        text(
          vmz[-c(refcol)],
          rel_int[-c(refcol)],
          lbmzcor[-c(refcol)],
          cex = 0.6,
          col = 1,
          srt = 90,
          adj = 0
        )
      }
    }

    text(
      vmz[refcol],
      rel_int[refcol],
      lbmzcor[refcol],
      cex = 0.8,
      col = 2,
      srt = 90,
      adj = 0
    )
  }

  ## prepare result file
  cor_valid <- (round(cor_abs_int, 2) >= processing_parameters$r_threshold)

  do_sirius <- TRUE
  verbose_catf("Checking sirius parameters...\n")
  if (is.null(processing_parameters$ionization)) {
    do_sirius <- FALSE
    verbose_catf("[KO] No ionization passed in parameter.\n")
  } else {
    verbose_catf("[OK] Ionization=%s.\n", processing_parameters$ionization)
  }
  if (is.na(processing_parameters$elemcomposition)) {
    do_sirius <- FALSE
    verbose_catf("[KO] Elemental composition is NA.\n")
  } else if (length(processing_parameters$elemcomposition) < 1) {
    do_sirius <- FALSE
    verbose_catf("[KO] No elemental composition is provided.\n")
  } else if (processing_parameters$elemcomposition == "") {
    do_sirius <- FALSE
    verbose_catf("[KO] Elemental composition is an empty string.\n")
  } else {
    verbose_catf(
      "[OK] Elemental composition=%s.\n",
      processing_parameters$elemcomposition
    )
  }

  cp_res_length <- length(vmz)
  errors <- rep(NA, cp_res_length)
  formulas <- rep(NA, cp_res_length)
  if (do_sirius) {
    verbose_catf("Everything is ok, preparing for sirius.\n")
    formulas <- get_formulas(
      mzref = processing_parameters$mzref,
      spectra = data.frame(mz = meaned_mz, intensities = sum_int[-1]),
      nominal_mz_list = vmz,
      processing_parameters = processing_parameters
    )
    if (nrow(formulas) == 0) {
      catf("No formula found.\n")
    } else {
      errors <- formulas$error
      formulas <- formulas$formula
      catf(
        "Found %s formula for %s fragments\n",
        length(formulas[which(!(is.na(formulas)))]),
        cp_res_length
      )
    }
  } else {
    verbose_catf("Sirius cannot be run.\n")
  }
  cp_res <- data.frame(
    rep(processing_parameters$c_name, cp_res_length),
    rep(processing_parameters$inchikey, cp_res_length),
    rep(processing_parameters$elemcomposition, cp_res_length),
    formulas,
    errors,
    rep(fid, cp_res_length),
    vmz,
    cor_abs_int,
    sum_int[-1],
    rel_int,
    cor_valid
  )

  colnames(cp_res) <- c(
    "compoundName",
    "inchikey",
    "elemcomposition",
    "fragment",
    "fragment_error",
    "fileid",
    "fragments_mz",
    "CorWithPrecursor",
    "AbsoluteIntensity",
    "relativeIntensity",
    "corValid"
  )
  return(cp_res)
}

#'
#' @title extract_fragments
#'
#' @param precursors the precursor list from mspurity
#' @param fragments the fragments list from ms purity
# ' @param mzref
# ' @param rtref
# ' @param c_name
# ' @param inchikey
# ' @param elemcomposition
#' @param processing_parameters
#' @returns
#'
#' @description
#' function for extraction of fragments corresponding to precursors
#' detected by MSPurity
extract_fragments <- function( ## nolint cyclocomp_linter
  precursors,
  fragments,
  # mzref,
  # rtref,
  # c_name,
  # inchikey,
  # elemcomposition,
  processing_parameters
) {
  ## filter precursor in the precursors file based on mz and rt in the
  ## compound list
  catf("processing %s\n", processing_parameters$c_name)
  verbose_catf("===\n")
  param <- processing_parameters
  selected_precursors <- which(
    (abs(precursors$precurMtchMZ - param$mzref) <= param$tolmz)
    & (abs(precursors$precurMtchRT - param$rtref) <= param$tolrt)
  )
  rm(param)

  verbose_catf(
    "> %s precursors selected with mz=%s±%s and rt=%s±%s\n",
    length(selected_precursors),
    processing_parameters$mzref,
    processing_parameters$tolmz,
    processing_parameters$rtref,
    processing_parameters$tolrt
  )

  ## check if there is the precursor in the file

  if (length(selected_precursors) < 1) {
    cat("> non detected in precursor file\n")
    show_end_processing()
    return(NULL)
  }

  precursors <- precursors[selected_precursors, ]

  ## check if fragments corresponding to precursor are found in several
  ## files (collision energy)
  ## this lead to a processing for each fileid
  file_ids <- as.character(sort(unique(precursors$fileid)))
  if (length(file_ids) > 1) {
    catf("> several files detected for this compounds :\n")
  } else if (length(file_ids) < 1 || nrow(precursors) < 1) {
    return(data.frame())
  }

  res_comp <- data.frame()
  for (curent_file_id in file_ids) {
    curent_precursors <- precursors[precursors$fileid == curent_file_id, ]
    selected_fragments <- fragments[
      fragments$grpid %in% as.character(curent_precursors$grpid)
      & fragments$fileid == curent_file_id,
    ]
    filtered_fragments <- selected_fragments[
      selected_fragments$ra > processing_parameters$seuil_ra,
    ]
    res_comp_by_file <- process_file(
      curent_file_id = curent_file_id,
      precursor_mz = curent_precursors$mz,
      filtered_fragments = filtered_fragments,
      processing_parameters = processing_parameters
    )
    if (!is.null(res_comp_by_file)) {
      res_comp <- rbind(res_comp, res_comp_by_file)
    }
  }
  return(unique(res_comp))
}


process_file <- function(
  curent_file_id,
  precursor_mz,
  filtered_fragments,
  processing_parameters
) {
  mznominal <- round(x = filtered_fragments$mz, digits = 0)
  meaned_mz <- round(
    aggregate(
      data.frame(
        mz = filtered_fragments$mz,
        mznominal = mznominal
      ),
      list(mznominal),
      FUN = mean
    )$mz,
    digits = processing_parameters$mzdecimal
  )
  filtered_fragments <- data.frame(filtered_fragments, mznominal)

  ## creation of cross table row=scan col=mz X=ra

  vmz <- as.character(sort(unique(filtered_fragments$mznominal)))

  ds_abs_int <- create_ds_abs_int(vmz, filtered_fragments)

  if (global_debug) {
    print(ds_abs_int)
    # write.table(
    #   x = ds_abs_int,
    #   file = paste0(c_name, "ds_abs_int.txt"),
    #   row.names = FALSE,
    #   sep = "\t"
    # )
  }

  ## elimination of mz with less than min_number_scan scans (user defined
  ## parameter)
  xmz <- rep(NA, ncol(ds_abs_int) - 1)
  sum_int <- rep(NA, ncol(ds_abs_int))
  nbxmz <- 0
  nb_scan_check <- min(nrow(ds_abs_int), processing_parameters$min_number_scan)

  for (j in 2:ncol(ds_abs_int)) {
    sum_int[j] <- sum(ds_abs_int[j], na.rm = TRUE)
    if (sum(!is.na(ds_abs_int[[j]])) < nb_scan_check) {
      nbxmz <- nbxmz + 1
      xmz[nbxmz] <- j
    }
  }

  xmz <- xmz[-which(is.na(xmz))]
  if (length(xmz) > 0) {
    ds_abs_int <- ds_abs_int[, -c(xmz)]
    sum_int <- sum_int[-c(xmz)]
    ## liste des mz keeped decale de 1 avec ds_abs_int
    vmz <- as.numeric(vmz[-c(xmz - 1)])
    meaned_mz <- meaned_mz[-c(xmz - 1)]
  }

  ## mz of precursor in data precursor to check correlation with
  mz_prec <- paste0(
    "mz",
    round(mean(precursor_mz), processing_parameters$mzdecimal)
  )
  ## reference ion for correlation computing = precursor OR maximum
  ## intensity ion in precursor is not present
  refcol <- which(colnames(ds_abs_int) == mz_prec)
  if (length(refcol) == 0) {
    refcol <- which(sum_int == max(sum_int, na.rm = TRUE))
  }

  if (processing_parameters$do_pdf) {
    start_pdf(processing_parameters$c_name, curent_file_id)
  }

  ## Pearson correlations between absolute intensities computing
  cor_abs_int <- rep(NA, length(vmz))

  if (length(refcol) > 0) {
    for (i in 2:length(ds_abs_int)) {
      cor_abs_int[i - 1] <- stats::cor(
        x = ds_abs_int[[refcol]],
        y = ds_abs_int[[i]],
        use = "pairwise.complete.obs",
        method = "pearson"
      )
      debug_catf(
        "Correlation between %s and %s: %s\n",
        paste(ds_abs_int[[refcol]], collapse = ";"),
        paste(ds_abs_int[[i]], collapse = ";"),
        paste(cor_abs_int[i - 1], collapse = ";")
      )
      if (processing_parameters$do_pdf) {
        pdf_plot_ds_abs_int(
          processing_parameters$c_name,
          ds_abs_int,
          refcol,
          i,
          round(cor_abs_int[i - 1], 2)
        )
      }
    }
    ## plot pseudo spectra
    res_comp_by_file <- plot_pseudo_spectra(
      x = ds_abs_int,
      fid = curent_file_id,
      sum_int = sum_int,
      vmz = vmz,
      cor_abs_int = cor_abs_int,
      refcol = refcol,
      meaned_mz = meaned_mz,
      processing_parameters = processing_parameters
    )
    catf(
      "%s has been processed and %s fragments have been found.\n",
      processing_parameters$c_name,
      nrow(res_comp_by_file)
    )
  } else {
    res_comp_by_file <- NULL
    cat(">> non detected in fragments file \n")
  }
  show_end_processing()
  if (processing_parameters$do_pdf) {
    end_pdf()
  }
  return(res_comp_by_file)
}

create_ds_abs_int <- function(vmz, filtered_fragments) {
  verbose_catf(
    ">> fragments: %s\n",
    paste(vmz, collapse = " ")
  )
  ds_abs_int <- create_int_mz(vmz[1], filtered_fragments)
  for (mz in vmz[-1]) {
    int_mz <- create_int_mz(mz, filtered_fragments)
    ds_abs_int <- merge(
      x = ds_abs_int,
      y = int_mz,
      by.x = 1,
      by.y = 1,
      all.x = TRUE,
      all.y = TRUE
    )
  }
  return(ds_abs_int)
}

create_int_mz <- function(mz, filtered_fragments) {
  ## absolute intensity
  int_mz <- filtered_fragments[
    filtered_fragments$mznominal == mz,
    c("acquisitionNum", "i")
  ]
  colnames(int_mz)[2] <- paste0("mz", mz)
  ## average intensities of mass in duplicate scans
  comp_scans <- aggregate(x = int_mz, by = list(int_mz[[1]]), FUN = mean)
  return(comp_scans[, -1])
}

show_end_processing <- function() {
  verbose_catf("==========\n")
  cat("\n")
}

start_pdf <- function(c_name, curent_file_id) {
  pdf(
    file = sprintf("%s_processing_file%s.pdf", c_name, curent_file_id),
    width = 8,
    height = 11
  )
  par(mfrow = c(3, 2))
}

pdf_plot_ds_abs_int <- function(c_name, ds_abs_int, refcol, i, r_coef) {
  plot(
    ds_abs_int[[refcol]],
    ds_abs_int[[i]],
     xlab = colnames(ds_abs_int)[refcol],
     ylab = colnames(ds_abs_int)[i],
     main = sprintf(
      "%s corr coeff r=%s", c_name, r_coef
    )
  )
}
end_pdf <- function() {
  dev.off()
}

set_global <- function(var, value) {
  assign(var, value, envir = globalenv())
}

set_debug <- function() {
  set_global("global_debug", TRUE)
}

unset_debug <- function() {
  set_global("global_debug", FALSE)
}

set_verbose <- function() {
  set_global("global_verbose", TRUE)
}

unset_verbose <- function() {
  set_global("global_verbose", FALSE)
}

verbose_catf <- function(...) {
  if (global_verbose) {
    cat(sprintf(...))
  }
}


debug_catf <- function(...) {
  if (global_debug) {
    cat(sprintf(...))
  }
}

catf <- function(...) {
  cat(sprintf(...))
}

create_parser <- function() {
  parser <- optparse::OptionParser()
  parser <- optparse::add_option(
    parser,
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = paste(
      "[default %default]",
      "Print extra output"
    )
  )
  parser <- optparse::add_option(
    parser,
    c("-V", "--version"),
    action = "store_true",
    default = FALSE,
    help = "Prints version and exits"
  )
  parser <- optparse::add_option(
    parser,
    c("-d", "--debug"),
    action = "store_true",
    default = FALSE,
    help = paste(
      "[default %default]",
      "Print debug outputs"
    )
  )
  parser <- optparse::add_option(
    parser,
    c("-o", "--output"),
    type = "character",
    default = DEFAULT_OUTPUT_PATH,
    action = "store",
    help = "Path to the output file [default %default]"
  )
  parser <- optparse::add_option(
    parser,
    c("-p", "--precursors"),
    type = "character",
    default = DEFAULT_PRECURSOR_PATH,
    action = "store",
    help = "Path to the precursors file [default %default]"
  )
  parser <- optparse::add_option(
    parser,
    c("-f", "--fragments"),
    type = "character",
    default = DEFAULT_FRAGMENTS_PATH,
    action = "store",
    help = "Path to the fragments file [default %default]"
  )
  parser <- optparse::add_option(
    parser,
    c("-c", "--compounds"),
    type = "character",
    default = DEFAULT_COMPOUNDS_PATH,
    action = "store",
    help = "Path to the compounds file [default %default]"
  )
  parser <- optparse::add_option(
    parser,
    c("--tolmz"),
    type = "numeric",
    action = "store",
    default = DEFAULT_TOLMZ,
    metavar = "number",
    help = paste(
      "[default %default]",
      "Tolerance for MZ (in Dalton) to match the standard in the compounds"
    )
  )
  parser <- optparse::add_option(
    parser,
    c("--tolrt"),
    type = "integer",
    action = "store",
    default = DEFAULT_TOLRT,
    metavar = "number",
    help = paste(
      "[default %default]",
      "RT (in seconds) to match the standard in the compounds"
    )
  )
  parser <- optparse::add_option(
    parser,
    c("--seuil_ra"),
    type = "numeric",
    action = "store",
    default = DEFAULT_SEUIL_RA,
    metavar = "number",
    help = paste(
      "[default %default]",
      "relative intensity threshold"
    ),
  )
  parser <- optparse::add_option(
    parser,
    c("--mzdecimal"),
    type = "integer",
    default = DEFAULT_MZDECIMAL,
    action = "store",
    help = paste(
      "[default %default]",
      "Number of decimal to write for MZ"
    ),
    metavar = "number"
  )
  parser <- optparse::add_option(
    parser,
    c("--r_threshold"),
    type = "integer",
    default = DEFAULT_R_THRESHOLD,
    action = "store",
    help = paste(
      "[default %default]",
      "R-Pearson correlation threshold between precursor and fragment",
      "absolute intensity"
    ),
    metavar = "number"
  )
  parser <- optparse::add_option(
    parser,
    c("--min_number_scan"),
    type = "numeric",
    action = "store",
    default = DEFAULT_MINNUMBERSCAN,
    help = paste(
      "[default %default]",
      "Fragments are kept if there are found in a minimum number",
      "of min_number_scan scans"
    ),
    metavar = "number"
  )
  parser <- optparse::add_option(
    parser,
    c("--ionization"),
    type = "character",
    action = "store",
    default = "None",
    help = paste(
      "[default %default]",
      "Which ionization to use for sirius"
    ),
    metavar = "character"
  )
  return(parser)
}

stop_with_status <- function(msg, status) {
  sink(stderr())
  message(sprintf("Error: %s", msg))
  message(sprintf("Error code: %s", status))
  sink(NULL)
  base::quit(status = status)
}

check_args_validity <- function(args) { ## nolint cyclocomp_linter
  if (length(args$output) == 0 || nchar(args$output[1]) == 0) {
    stop_with_status(
      "Missing output parameters. Please set it with --output.",
      MISSING_PARAMETER_ERROR
    )
  }
  if (length(args$precursors) == 0 || nchar(args$precursors[1]) == 0) {
    stop_with_status(
      "Missing precursors parameters. Please set it with --precursors.",
      MISSING_PARAMETER_ERROR
    )
  }
  if (length(args$fragments) == 0 || nchar(args$fragments[1]) == 0) {
    stop_with_status(
      "Missing fragments parameters. Please set it with --fragments.",
      MISSING_PARAMETER_ERROR
    )
  }
  if (length(args$compounds) == 0 || nchar(args$compounds[1]) == 0) {
    stop_with_status(
      "Missing compounds parameters. Please set it with --compounds.",
      MISSING_PARAMETER_ERROR
    )
  }
  if (!file.exists(args$precursors)) {
    stop_with_status(
      sprintf(
        "Precursors file %s does not exist or cannot be accessed.",
        args$precursors
      ),
      MISSING_INPUT_FILE_ERROR
    )
  }
  if (!file.exists(args$fragments)) {
    stop_with_status(
      sprintf(
        "Fragments file %s does not exist or cannot be accessed.",
        args$fragments
      ),
      MISSING_INPUT_FILE_ERROR
    )
  }
  if (!file.exists(args$compounds)) {
    stop_with_status(
      sprintf(
        "Compounds file %s does not exist or cannot be accessed.",
        args$compounds
      ),
      MISSING_INPUT_FILE_ERROR
    )
  }
  if (in_galaxy_env()) {
    check_galaxy_args_validity(args)
  }
}

in_galaxy_env <- function() {
  sysvars <- Sys.getenv()
  sysvarnames <- names(sysvars)
  return(
    "_GALAXY_JOB_HOME_DIR" %in% sysvarnames
    || "_GALAXY_JOB_TMP_DIR" %in% sysvarnames
    || "GALAXY_MEMORY_MB" %in% sysvarnames
    || "GALAXY_MEMORY_MB_PER_SLOT" %in% sysvarnames
    || "GALAXY_SLOTS" %in% sysvarnames
  )
}

check_galaxy_args_validity <- function(args) {
  if (!file.exists(args$output)) {
    stop_with_status(
      sprintf(
        "Output file %s does not exist or cannot be accessed.",
        args$output
      ),
      MISSING_INPUT_FILE_ERROR
    )
  }
}

get_csv_or_tsv <- function(
  path,
  sep_stack = c("\t", ",", ";"),
  sep_names = c("tab", "comma", "semicolon"),
  header = TRUE,
  quote = "\""
) {
  sep <- determine_csv_or_tsv_sep(
    path = path,
    sep_stack = sep_stack,
    header = header,
    quote = quote
  )
  verbose_catf(
    "%s separator has been determined for %s.\n",
    sep_names[sep_stack == sep],
    path
  )
  return(read.table(
    file = path,
    sep = sep,
    header = header,
    quote = quote
  ))
}

determine_csv_or_tsv_sep <- function(
  path,
  sep_stack = c("\t", ",", ";"),
  header = TRUE,
  quote = "\""
) {
  count <- -1
  best_sep <- sep_stack[1]
  for (sep in sep_stack) {
    tryCatch({
      table <- read.table(
        file = path,
        sep = sep,
        header = header,
        quote = quote,
        nrows = 1
      )
      if (ncol(table) > count) {
        count <- ncol(table)
        best_sep <- sep
      }
    })
  }
  return(best_sep)
}

uniformize_columns <- function(df) {
  cols <- colnames(df)
  for (func in c(tolower)) {
    cols <- func(cols)
  }
  colnames(df) <- cols
  return(df)
}

handle_galaxy_param <- function(args) {
  for (param in names(args)) {
    if (is.character(args[[param]])) {
      args[[param]] <- gsub("__ob__", "[", args[[param]])
      args[[param]] <- gsub("__cb__", "]", args[[param]])
    }
  }
  return(args)
}

main <- function(args) {
  if (args$version) {
    catf("%s\n", MS2SNOOP_VERSION)
    base::quit(status = 0)
  }
  if (in_galaxy_env()) {
    print(sessionInfo())
    cat("\n\n")
  }
  check_args_validity(args)
  args <- handle_galaxy_param(args)
  if (args$ionization == "None") {
    args$ionization <- NULL
  }
  if (args$debug) {
    set_debug()
  }
  if (args$verbose) {
    set_verbose()
  }
  precursors <- get_csv_or_tsv(args$precursors)
  fragments <- get_csv_or_tsv(args$fragments)
  compounds <- get_csv_or_tsv(args$compounds)

  compounds <- uniformize_columns(compounds)
  mandatory_columns <- c(
    "compound_name",
    "mz",
    "rtsec",
    "inchikey"
  )
  presents <- mandatory_columns %in% colnames(compounds)
  if (!all(presents)) {
    stop_with_status(
      sprintf(
        "Some columns are missing: %s",
        paste(mandatory_columns[which(!presents)], collapse = ", ")
      ),
      BAD_PARAMETER_VALUE_ERROR
    )
  }

  res_all <- data.frame()
  processing_parameters <- list(
    min_number_scan = args$min_number_scan,
    mzdecimal = args$mzdecimal,
    r_threshold = args$r_threshold,
    seuil_ra = args$seuil_ra,
    tolmz = args$tolmz,
    tolrt = args$tolrt,
    ionization = args$ionization,
    do_pdf = FALSE
  )
  for (i in seq_len(nrow(compounds))) {
    processing_parameters$mzref <- compounds[["mz"]][i]
    processing_parameters$rtref <- compounds[["rtsec"]][i]
    processing_parameters$c_name <- compounds[["compound_name"]][i]
    processing_parameters$inchikey <- compounds[["inchikey"]][i]
    processing_parameters$elemcomposition <- compounds[["elemcomposition"]][i]
    res_cor <- extract_fragments(
      precursors = precursors,
      fragments = fragments,
      # mzref = compounds[["mz"]][i],
      # rtref = compounds[["rtsec"]][i],
      # c_name = compounds[["compound_name"]][i],
      # inchikey = compounds[["inchikey"]][i],
      # elemcomposition = compounds[["elemcomposition"]][i],
      processing_parameters = processing_parameters
    )
    if (!is.null(res_cor)) {
      res_all <- rbind(res_all, res_cor)
    }
  }

  if (nrow(res_all) == 0) {
    stop_with_status("No result at all!", NO_ANY_RESULT_ERROR)
  }

  write.table(
    x = res_all,
    file = args$output,
    sep = "\t",
    row.names = FALSE
  )
}

global_debug <- FALSE
global_verbose <- FALSE
args <- optparse::parse_args(create_parser())
main(args)

warnings()
