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


assign("MS2SNOOP_VERSION", "1.0.1")
lockBinding("MS2SNOOP_VERSION", globalenv())

assign("DEFAULT_PRECURSOR_PATH", "peaklist_precursors.tsv")
assign("DEFAULT_FRAGMENTS_PATH", "peaklist_fragments.tsv")
assign("DEFAULT_COMPOUNDS_PATH", "compounds_pos.txt")
assign("DEFAULT_OUTPUT_PATH", "compound_fragments_result.txt")
assign("DEFAULT_TOLMZ", 0.01)
assign("DEFAULT_TOLRT", 20)
assign("DEFAULT_MZDECIMAL", 0)
assign("DEFAULT_R_THRESHOLD", 0.85)
assign("DEFAULT_MINNUMBERSCAN", 8)
assign("DEFAULT_SEUIL_RA", 0.5)
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


debug <- FALSE


########################################################################

#' @title plot_pseudo_spectra
#' @param x
#' @param r_threshold
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
  r_threshold,
  fid,
  sum_int,
  vmz,
  cor_abs_int,
  refcol,
  c_name
) {
  ## du fait de la difference de nombre de colonne entre la dataframe qui
  ## inclue les scans en 1ere col, mzRef se decale de 1
  refcol <- refcol - 1
  ## compute relative intensities max=100%
  rel_int <- sum_int[-1]
  rel_int <- rel_int / max(rel_int)

  ## define max value on vertical axis (need to increase in order to plot the
  ## label of fragments)
  ymax <- max(rel_int) + 0.2 * max(rel_int)

  par(mfrow = c(2, 1))
  plot(vmz, rel_int, type = "h", ylim = c(0, ymax), main = c_name)
  ## low correl coef. will be display in grey
  cor_low <- which(round(cor_abs_int, 2) < r_threshold)

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

  ## prepare result file
  corValid <- (round(cor_abs_int, 2) >= r_threshold) ##nolint object_name_linter
  cp_res <- data.frame(
    rep(c_name, length(vmz)),
    rep(fid, length(vmz)),
    vmz,
    cor_abs_int,
    sum_int[-1],
    rel_int,
    corValid
  )

  colnames(cp_res) <- c(
    "compoundName",
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
#' @param mzref
#' @param rtref
#' @param c_name
#' @param r_threshold default = DEFAULT_EXTRACT_FRAGMENTS_R_THRESHOLD
#' @param seuil_ra default = DEFAULT_EXTRACT_FRAGMENTS_SEUIL_RA
#' @param tolmz default = DEFAULT_EXTRACT_FRAGMENTS_TOLMZ
#' @param tolrt default = DEFAULT_EXTRACT_FRAGMENTS_TOLRT
#' @returns
#'
#' @description
#' function for extraction of fragments corresponding to precursors
#' detected by MSPurity
extract_fragments <- function( ## nolint cyclocomp_linter
  precursors,
  fragments,
  mzref,
  rtref,
  c_name,
  min_number_scan,
  mzdecimal,
  r_threshold=DEFAULT_EXTRACT_FRAGMENTS_R_THRESHOLD,
  seuil_ra=DEFAULT_EXTRACT_FRAGMENTS_SEUIL_RA,
  tolmz=DEFAULT_EXTRACT_FRAGMENTS_TOLMZ,
  tolrt=DEFAULT_EXTRACT_FRAGMENTS_TOLRT
) {
  ## filter precursor in the precursors file based on mz and rt in the
  ## compound list
  cat("processing ", c_name, "\n")
  selected_precursors <- which(
    (abs(precursors$precurMtchMZ - mzref) <= tolmz)
    & (abs(precursors$precurMtchRT - rtref) <= tolrt)
  )

  ## check if there is the precursor in the file
  if (length(selected_precursors) > 0) {

    sprecini <- precursors[selected_precursors, ]

    ## check if fragments corresponding to precursor are found in several
    ## files (collision energy)
    ## this lead to a processing for each fileid
    mf <- levels(as.factor(sprecini$fileid))
    if (length(mf) > 1) {
      cat(" several files detected for this compounds :\n")
    }

    for (f in seq_along(mf)) {

      sprec <- sprecini[sprecini$fileid == mf[f], ]

      ## selection of fragment in the fragments file with the grpid common in
      ## both fragments and precursors
      selfrgt <- levels(as.factor(sprec$grpid))
      sfrgt <- fragments[
        fragments$grpid %in% selfrgt
        & fragments$fileid == mf[f],
      ]

      ## filter fragments on relative intensity seuil_ra = user defined
      ## parameter (MSpurity flags could be used here)
      sfrgtfil <- sfrgt[sfrgt$ra > seuil_ra, ]

      mznominal <- round(x = sfrgtfil$mz, mzdecimal)
      sfrgtfil <- data.frame(sfrgtfil, mznominal)

      ## creation of cross table row=scan col=mz X=ra
      vmz <- levels(as.factor(sfrgtfil$mznominal))

      cat(" fragments :", vmz)

      ## mz of precursor in data precursor to check correlation with
      mz_prec <- paste0("mz", round(mean(sprec$mz), mzdecimal))

      for (m in seq_along(vmz)) {

        ## absolute intensity
        cln <- c(
          which(colnames(sfrgtfil) == "acquisitionNum"),
          which(colnames(sfrgtfil) == "i")
        )
        int_mz <- sfrgtfil[sfrgtfil$mznominal == vmz[m], cln]
        colnames(int_mz)[2] <- paste0("mz", vmz[m])

        ## average intensities of mass in duplicate scans
        comp_scans <- aggregate(x = int_mz, by = list(int_mz[[1]]), FUN = mean)
        int_mz <- comp_scans[, -1]

        if (m == 1) {
          ds_abs_int <- int_mz
        } else {
          ds_abs_int <- merge(
            x = ds_abs_int,
            y = int_mz,
            by.x = 1,
            by.y = 1,
            all.x = TRUE,
            all.y = TRUE
          )
        }
      }
      if (debug) {
        write.table(
          x = ds_abs_int,
          file = paste0(c_name, "ds_abs_int.txt"),
          row.names = FALSE,
          sep = "\t"
        )
      }

      ## elimination of mz with less than min_number_scan scans (user defined
      ## parameter)
      xmz <- rep(NA, ncol(ds_abs_int) - 1)
      sum_int <- rep(NA, ncol(ds_abs_int))
      nbxmz <- 0
      nb_scan_check <- min(nrow(ds_abs_int), min_number_scan)

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
      }

      ## reference ion for correlation computing = precursor OR maximum
      ## intensity ion in precursor is not present
      refcol <- which(colnames(ds_abs_int) == mz_prec)
      if (length(refcol) == 0) {
        refcol <- which(sum_int == max(sum_int, na.rm = TRUE))
      }
      pdf(
        file = sprintf("%s_processing_file%s.pdf", c_name, mf[f]),
        width = 8,
        height = 11
      )
      par(mfrow = c(3, 2))

      ## Pearson correlations between absolute intensities computing
      cor_abs_int <- rep(NA, length(vmz))

      if (length(refcol) > 0) {
        for (i in 2:length(ds_abs_int)) {
          cor_abs_int[i - 1] <- cor(
            x = ds_abs_int[[refcol]],
            y = ds_abs_int[[i]],
            use = "pairwise.complete.obs",
            method = "pearson"
          )
          plot(
            ds_abs_int[[refcol]],
            ds_abs_int[[i]],
             xlab = colnames(ds_abs_int)[refcol],
             ylab = colnames(ds_abs_int)[i],
             main = sprintf(
              "%s corr coeff r=%s", c_name, round(cor_abs_int[i - 1], 2)
            )
          )
        }
        ## plot pseudo spectra
        res_comp_by_file <- plot_pseudo_spectra(
          x = ds_abs_int,
          r_threshold = r_threshold,
          fid = mf[f],
          sum_int = sum_int,
          vmz = vmz,
          cor_abs_int = cor_abs_int,
          refcol = refcol,
          c_name = c_name
        )
        if (f == 1) {
          res_comp <- res_comp_by_file
        }
      } else {
        res_comp_by_file <- NULL
        cat(" non detected in fragments file \n")
      }
      if (!is.null(res_comp_by_file)) {
        res_comp <- rbind(res_comp, res_comp_by_file)
      }
      cat("\n")
      dev.off()
    }
  } else {
    res_comp <- NULL
    cat(" non detected in precursor file \n")
  }
  return(res_comp)
}


create_parser <- function() {
  parser <- optparse::OptionParser()
  parser <- optparse::add_option(
    parser,
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Print extra output [default %default]"
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
    metavar = "number"
  )
  parser <- optparse::add_option(
    parser,
    c("--tolrt"),
    type = "integer",
    action = "store",
    default = DEFAULT_TOLRT,
    metavar = "number"
  )
  parser <- optparse::add_option(
    parser,
    c("--seuil_ra"),
    type = "numeric",
    action = "store",
    help = "relative intensity threshold",
    default = DEFAULT_SEUIL_RA,
    metavar = "number"
  )
  parser <- optparse::add_option(
    parser,
    c("--mzdecimal"),
    type = "integer",
    default = DEFAULT_MZDECIMAL,
    action = "store",
    help = "nb decimal for mz",
    metavar = "number"
  )
  parser <- optparse::add_option(
    parser,
    c("--r_threshold"),
    type = "integer",
    default = DEFAULT_R_THRESHOLD,
    action = "store",
    help = paste0(
      "r pearson correlation threshold between precursor and fragment ",
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
    help = paste0(
      "fragments are kept if there are found in a minimum number ",
      "of scans"
    ),
    metavar = "number"
  )
  return(parser)
}

main <- function(args) {
  if (args$version) {
    cat(sprintf("%s\n", MS2SNOOP_VERSION))
    base::quit(status = 0)
  }
  sessionInfo()
  ## FOLDER AND FILES
  ## MSpurity precursors file
  precursors <- read.table(
    file = args$precursors,
    header = TRUE,
    sep = "\t",
    quote = "\""
  )
  ## MSpurity fragments file
  fragments <- read.table(
    file = args$fragments,
    header = TRUE,
    sep = "\t",
    quote = "\""
  )
  ## list of compounds : col1=Name of molecule, col2=m/z, col3=retention time
  compounds <- read.table(
    file = args$compounds,
    sep = "\t",
    quote = "\"",
    header = TRUE
  )
  ## PARAMETERS
  ## tolerance for mz(dalton) rt(seconds) to match the standard in the compounds
  ## list with the precursor MSpurity file
  tolmz <- args$tolmz
  tolrt <- args$tolrt

  ##  relative intensity threshold
  seuil_ra <- args$seuil_ra
  ## nb decimal for mz
  mzdecimal <- args$mzdecimal
  ## r pearson correlation threshold between precursor and
  # #fragment absolute intensity
  r_threshold <- args$r_threshold
  ## fragments are kept if there are found in a minimum number of scans
  min_number_scan <- args$min_number_scan

  res_all <- NULL
  for (i in seq_len(nrow(compounds))) {
    ## loop execution for all compounds in the compounds file
    res_cor <- NULL
    res_cor <- extract_fragments(
      precursors = precursors,
      fragments = fragments,
      mzref = compounds[[2]][i],
      rtref = compounds[[3]][i],
      c_name = compounds[[1]][i],
      min_number_scan = min_number_scan,
      mzdecimal = mzdecimal,
      r_threshold = r_threshold,
      seuil_ra = seuil_ra,
      tolmz = tolmz,
      tolrt = tolrt
    )
    if (!is.null(res_cor)) {
      if (is.null(res_all)) {
        res_all <- res_cor
      } else {
        res_all <- rbind(res_all, res_cor)
      }
    }
  }

  if (is.null(res_all)) {
    stop("No result at all!")
  }
  write.table(
    x = res_all,
    file = args$output,
    sep = "\t",
    row.names = FALSE
  )
}

args <- optparse::parse_args(create_parser())
main(args)

warnings()
