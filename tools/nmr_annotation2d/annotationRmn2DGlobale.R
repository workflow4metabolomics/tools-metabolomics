#########################################################################
# ANNOTATION SPECTRE 2D MATRICE COMPLEXE BASEE SUR UNE (OU PLUSIEURS)   #
# SEQUENCE (s)                                                          #
# template : dataframe contenant la liste des couples de deplacements   #
#   chimiques de la matrice complexe a annoter                          #
# cosy : 1 si sequence a utiliser / 0 sinon                             #
# hmbc : 1 si sequence a utiliser / 0 sinon                             #
# hsqc : 1 si sequence a utiliser / 0  sinon                            #
# jres : 1 si sequence a utiliser / 0 sinon                             #
# tocsy : 1 si sequence a utiliser / 0 sinon                            #
# tolPpm1 : tolerance autorisee autour de la valeur1 du couple de       #
#   deplacements chimiques                                              #
# tolPpm2HJRes : tolerance autorisee autour de la valeur2 du couple de  #
#   deplacements chimiques si H dans dimension 2                        #
# tolPpm2C : tolerance autorisee autour de la valeur2 du couple de      #
#   deplacements chimiques si C dans dimension 2                        #
# seuil : valeur du score de presence en dela de laquelle les           #
#   metabolites annotes ne sont pas retenus                             #
# unicite : boolean pour ne retenir que les ...                         #
#########################################################################
## CALCUL MOYENNE SANS VALEUR(S) MANQUANTE(S)
mean_rm_na <- function(x) {
  mean(x, na.rm = TRUE)
}

annotationRmn2DGlobale <- function( ## nolint
  template,               ## nolint
  tolPpm1 = 0.01,         ## nolint
  tolPpm2HJRes = 0.002,   ## nolint
  tolPpm2C = 0.5,         ## nolint
  cosy = 1,               ## nolint
  hmbc = 1,               ## nolint
  hsqc = 1,               ## nolint
  jres = 1,               ## nolint
  tocsy = 1,              ## nolint
  seuil,                  ## nolint
  unicite = "NO"          ## nolint
) {
  ## varnames cleanup without modifying signature
  tol_ppm1 <- tolPpm1
  tol_ppm2_hj_res <- tolPpm2HJRes
  tol_ppm2_c <- tolPpm2C

  ## Initialisation
  options(max.print = 999999999)
  annotation_cosy <- data.frame()
  annotation_hbmc <- data.frame()
  annotation_hsqc <- data.frame()
  annotation_jres <- data.frame()
  annotation_tocsy <- data.frame()

  data_cosy <- "NA"
  data_hmbc <- "NA"
  data_hsqc <- "NA"
  data_jres <- "NA"
  data_tocsy <- "NA"

  ## Application seuil seulement si annotation avec 1 seule sequence
  seuil_pls_2d <- seuil

  if (cosy == 1) {
    matrice_cosy <- read.xlsx(
      template,
      sheet = "COSY",
      startRow = 2,
      colNames = TRUE,
      rowNames = FALSE,
      cols = 1:3,
      na.strings = "NA"
    )
    matrice_cosy <- matrice_cosy[matrice_cosy$peak.index != "x", ]
    annotation_cosy <- annotationRmn2D(
      matrice_cosy,
      BdDReference_COSY,
      "COSY",
      ppm1Tol = tolPpm1,
      ppm2Tol = tolPpm1,
      seuil = seuil_pls_2d,
      unicite = unicite
    )
    data_cosy <- data.frame(
      Metabolite = str_to_lower(annotation_cosy$liste_resultat$Metabolite),
      score.COSY = annotation_cosy$liste_resultat$score
    )
    data_cosy <- unique.data.frame(data_cosy)
  }

  if (hmbc == 1) {
    matrice_hmbc <- read.xlsx(
      template,
      sheet = "HMBC",
      startRow = 2,
      colNames = TRUE,
      rowNames = FALSE,
      cols = 1:3,
      na.strings = "NA"
    )
    matrice_hmbc <- matrice_hmbc[matrice_hmbc$peak.index != "x", ]
    annotation_hbmc <- annotationRmn2D(
      matrice_hmbc,
      BdDReference_HMBC,
      "HMBC",
      ppm1Tol = tolPpm1,
      ppm2Tol = tolPpm2C,
      seuil = seuil_pls_2d,
      unicite = unicite
    )
    data_hmbc <- data.frame(
      Metabolite = str_to_lower(annotation_hbmc$liste_resultat$Metabolite),
      score.HMBC = annotation_hbmc$liste_resultat$score
    )
    data_hmbc <- unique.data.frame(data_hmbc)
  }

  if (hsqc == 1) {
    matrice_hsqc <- read.xlsx(
      template,
      sheet = "HSQC",
      startRow = 2,
      colNames = TRUE,
      rowNames = FALSE,
      cols = 1:3,
      na.strings = "NA"
    )
    matrice_hsqc <- matrice_hsqc[matrice_hsqc$peak.index != "x", ]
    annotation_hsqc <- annotationRmn2D(
      matrice_hsqc,
      BdDReference_HSQC,
      "HSQC",
      ppm1Tol = tolPpm1,
      ppm2Tol = tolPpm2C,
      seuil = seuil_pls_2d,
      unicite = unicite
    )
    data_hsqc <- data.frame(
      Metabolite = str_to_lower(annotation_hsqc$liste_resultat$Metabolite),
      score.HSQC = annotation_hsqc$liste_resultat$score
    )
    data_hsqc <- unique.data.frame(data_hsqc)
  }

  if (jres == 1) {
    matrice_jres <- read.xlsx(
      template,
      sheet = "JRES",
      startRow = 2,
      colNames = TRUE,
      rowNames = FALSE,
      cols = 1:3,
      na.strings = "NA"
    )
    matrice_jres <- matrice_jres[matrice_jres$peak.index != "x", ]
    annotation_jres <- annotationRmn2D(
      matrice_jres,
      BdDReference_JRES,
      "JRES",
      ppm1Tol = tolPpm1,
      ppm2Tol = tolPpm2HJRes,
      seuil = seuil_pls_2d,
      unicite = unicite
    )
    data_jres <- data.frame(
      Metabolite = str_to_lower(annotation_jres$liste_resultat$Metabolite),
      score.JRES = annotation_jres$liste_resultat$score
    )
    data_jres <- unique.data.frame(data_jres)
  }

  if (tocsy == 1) {
    matrice_tocsy <- read.xlsx(
      template,
      sheet = "TOCSY",
      startRow = 2,
      colNames = TRUE,
      rowNames = FALSE,
      cols = 1:3,
      na.strings = "NA"
    )
    matrice_tocsy <- matrice_tocsy[matrice_tocsy$peak.index != "x", ]
    annotation_tocsy <- annotationRmn2D(
      matrice_tocsy,
      BdDReference_TOCSY,
      "TOCSY",
      ppm1Tol = tolPpm1,
      ppm2Tol = tolPpm1,
      seuil = seuil_pls_2d,
      unicite = unicite
    )
    data_tocsy <- data.frame(
      Metabolite = str_to_lower(annotation_tocsy$liste_resultat$Metabolite),
      score.TOCSY = annotation_tocsy$liste_resultat$score
    )
    data_tocsy <- unique.data.frame(data_tocsy)
  }

  ## CONCATENATION RESULTATS DIFFERENTES SEQUENCES
  data_2d <- list(data_cosy, data_hmbc, data_hsqc, data_jres, data_tocsy)
  which_sequence_nan <- which((data_2d != "NA"))
  data_2d <- data_2d[which_sequence_nan]
  sequences_combination <- data.frame(data_2d[1])
  seq_combi_nean_score <- sequences_combination

  ## Si une seule sequence et seuil sur score = filtre applique dans
  ## la fonction annotationRmn2D
  if (length(data_2d) >= 2) {
    ## CONCATENATION SCORE PAR SEQUENCE
    for (l in 2:length(data_2d)) {
      sequences_combination <- merge.data.frame(
        sequences_combination,
        data_2d[l],
        by = "Metabolite",
        all.x = TRUE,
        all.y = TRUE
      )
    }

    ## Replacement of NA values due to mis annotation
    for (m in seq_len(nrow(sequences_combination))) {
      if (cosy == 1)
        cosy_compound <- sort(names(BdDReference_COSY))
      if (hmbc == 1)
        hmbc_compound <- sort(names(BdDReference_HMBC))
      if (hsqc == 1)
        hsqc_compound <- sort(names(BdDReference_HSQC))
      if (jres == 1)
        jres_compound <- sort(names(BdDReference_JRES))
      if (tocsy == 1)
        tocsy_compound <- sort(names(BdDReference_TOCSY))

      if (!is.null(sequences_combination$score.COSY)) {
        if (is.na(sequences_combination[m, which(colnames(sequences_combination) == "score.COSY")])) {
          compound <- as.character(sequences_combination[m, 1])
          for (c in seq_len(length(cosy_compound)))
            if (str_to_lower(compound) == str_to_lower(cosy_compound[c]))
              sequences_combination[m, which(colnames(sequences_combination) == "score.COSY")] <- 0
        }
      }

      if (!is.null(sequences_combination$score.HMBC)) {
        if (is.na(sequences_combination[m, which(colnames(sequences_combination) == "score.HMBC")])) {
          compound <- as.character(sequences_combination[m, 1])
          for (c in seq_len(length(hmbc_compound)))
            if (str_to_lower(compound) == str_to_lower(hmbc_compound[c]))
              sequences_combination[m, which(colnames(sequences_combination) == "score.HMBC")] <- 0
        }
      }

      if (!is.null(sequences_combination$score.HSQC)) {
        if (is.na(sequences_combination[m, which(colnames(sequences_combination) == "score.HSQC")])) {
          compound <- as.character(sequences_combination[m, 1])
          for (c in seq_len(length(hsqc_compound)))
            if (str_to_lower(compound) == str_to_lower(hsqc_compound[c]))
               sequences_combination[m, which(colnames(sequences_combination) == "score.HSQC")] <- 0
        }
      }

      if (!is.null(sequences_combination$score.JRES)) {
        if (is.na(sequences_combination[m, which(colnames(sequences_combination) == "score.JRES")])) {
          compound <- as.character(sequences_combination[m, 1])
          for (c in seq_len(length(jres_compound)))
            if (str_to_lower(compound) == str_to_lower(jres_compound[c]))
              sequences_combination[m, which(colnames(sequences_combination) == "score.JRES")] <- 0
        }
      }

      if (!is.null(sequences_combination$score.TOCSY)) {
        if (is.na(sequences_combination[m, which(colnames(sequences_combination) == "score.TOCSY")])) {
          compound <- as.character(sequences_combination[m, 1])
          for (c in seq_len(length(tocsy_compound)))
            if (str_to_lower(compound) == str_to_lower(tocsy_compound[c]))
              sequences_combination[m, which(colnames(sequences_combination) == "score.TOCSY")] <- 0
        }
      }
    }

    ## SCORE MOYEN (sans prise en compte valeurs manquantes)
    mean_score <- round(
      apply(sequences_combination[, -1], 1, FUN = mean_rm_na),
      2
    )
    seq_combi_nean_score <- cbind.data.frame(
      sequences_combination,
      averageScore = mean_score
    )
  }
  return(list(
    COSY = annotation_cosy,
    HMBC = annotation_hbmc,
    HSQC = annotation_hsqc,
    JRES = annotation_jres,
    TOCSY = annotation_tocsy,
    ## SUPPRESSION METABOLITE AVEC SCORE MOYEN < SEUIL
    combination = seq_combi_nean_score[
      seq_combi_nean_score$averageScore > seuil,
    ]
  ))
}
