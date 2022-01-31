###################################################################################################
# ANNOTATION SPECTRE 2D MATRICE COMPLEXE BASEE SUR UNE (OU PLUSIEURS) SEQUENCE(s)                 #
# template : dataframe contenant la liste des couples de deplacements chimiques de la matrice complexe a annoter									                          #
# cosy : 1 si sequence a utiliser / 0 sinon                                                       #
# hmbc : 1 si sequence a utiliser / 0 sinon                                                       #
# hsqc : 1 si sequence a utiliser / 0  sinon                                                      #
# jres : 1 si sequence a utiliser / 0 sinon                                                       #
# tocsy : 1 si sequence a utiliser / 0 sinon                                                      #
# tolPpm1 : tolerance autorisee autour de la valeur1 du couple de deplacements chimiques          #
# tolPpm2HJRes : tolerance autorisee autour de la valeur2 du couple de deplacements chimiques si H dans dimension 2                       								          #
# tolPpm2C : tolerance autorisee autour de la valeur2 du couple de deplacements chimiques si C dans dimension 2                           									  #
# seuil : valeur du score de presence en dela de laquelle les metabolites annotes ne sont pas retenus #
# unicite : boolean pour ne retenir que les ...                                                   #
###################################################################################################
## CALCUL MOYENNE SANS VALEUR(S) MANQUANTE(S)
mean.rmNa <- function(x) {
  mean(x, na.rm = TRUE)
}

annotationRmn2DGlobale <- function(template, tolPpm1 = 0.01, tolPpm2HJRes = 0.002, tolPpm2C = 0.5, cosy = 1, hmbc = 1, hsqc = 1, jres = 1, tocsy = 1, seuil, unicite = "NO") {
  ## Initialisation
  options(max.print = 999999999)
  annotationCOSY <- data.frame()
  annotationHMBC <- data.frame()
  annotationHSQC <- data.frame()
  annotationJRES <- data.frame()
  annotationTOCSY <- data.frame()

  dataCOSY <- "NA"
  dataHMBC <- "NA"
  dataHSQC <- "NA"
  dataJRES <- "NA"
  dataTOCSY <- "NA"

  ## Application seuil seulement si annotation avec 1 seule sequence
  seuilPls2D <- seuil

  if (cosy == 1) {
    matrice.cosy <- read.xlsx(template, sheet = "COSY", startRow = 2, colNames = TRUE, rowNames = FALSE, cols = 1:3, na.strings = "NA")
    matrice.cosy <- matrice.cosy[matrice.cosy$peak.index != "x", ]
    annotationCOSY <- annotationRmn2D(matrice.cosy, BdDReference_COSY, "COSY", ppm1Tol = tolPpm1, ppm2Tol = tolPpm1, seuil = seuilPls2D, unicite = unicite)
    dataCOSY <- data.frame(Metabolite = str_to_lower(annotationCOSY$liste_resultat$Metabolite), score.COSY = annotationCOSY$liste_resultat$score)
    dataCOSY <- unique.data.frame(dataCOSY)
  }

  if (hmbc == 1) {
    matrice.hmbc <- read.xlsx(template, sheet = "HMBC", startRow = 2, colNames = TRUE, rowNames = FALSE, cols = 1:3, na.strings = "NA")
    matrice.hmbc <- matrice.hmbc[matrice.hmbc$peak.index != "x", ]
    annotationHMBC <- annotationRmn2D(matrice.hmbc, BdDReference_HMBC, "HMBC", ppm1Tol = tolPpm1, ppm2Tol = tolPpm2C, seuil = seuilPls2D, unicite = unicite)
    dataHMBC <- data.frame(Metabolite = str_to_lower(annotationHMBC$liste_resultat$Metabolite), score.HMBC = annotationHMBC$liste_resultat$score)
    dataHMBC <- unique.data.frame(dataHMBC)
  }

  if (hsqc == 1) {
    matrice.hsqc <- read.xlsx(template, sheet = "HSQC", startRow = 2, colNames = TRUE, rowNames = FALSE, cols = 1:3, na.strings = "NA")
    matrice.hsqc <- matrice.hsqc[matrice.hsqc$peak.index != "x", ]
    annotationHSQC <- annotationRmn2D(matrice.hsqc, BdDReference_HSQC, "HSQC", ppm1Tol = tolPpm1, ppm2Tol = tolPpm2C, seuil = seuilPls2D, unicite = unicite)
    dataHSQC <- data.frame(Metabolite = str_to_lower(annotationHSQC$liste_resultat$Metabolite), score.HSQC = annotationHSQC$liste_resultat$score)
    dataHSQC <- unique.data.frame(dataHSQC)
  }

  if (jres == 1) {
    matrice.jres <- read.xlsx(template, sheet = "JRES", startRow = 2, colNames = TRUE, rowNames = FALSE, cols = 1:3, na.strings = "NA")
    matrice.jres <- matrice.jres[matrice.jres$peak.index != "x", ]
    annotationJRES <- annotationRmn2D(matrice.jres, BdDReference_JRES, "JRES", ppm1Tol = tolPpm1, ppm2Tol = tolPpm2HJRes, seuil = seuilPls2D, unicite = unicite)
    dataJRES <- data.frame(Metabolite = str_to_lower(annotationJRES$liste_resultat$Metabolite), score.JRES = annotationJRES$liste_resultat$score)
    dataJRES <- unique.data.frame(dataJRES)
  }

  if (tocsy == 1) {
    matrice.tocsy <- read.xlsx(template, sheet = "TOCSY", startRow = 2, colNames = TRUE, rowNames = FALSE, cols = 1:3, na.strings = "NA")
    matrice.tocsy <- matrice.tocsy[matrice.tocsy$peak.index != "x", ]
    annotationTOCSY <- annotationRmn2D(matrice.tocsy, BdDReference_TOCSY, "TOCSY", ppm1Tol = tolPpm1, ppm2Tol = tolPpm1, seuil = seuilPls2D, unicite = unicite)
    dataTOCSY <- data.frame(Metabolite = str_to_lower(annotationTOCSY$liste_resultat$Metabolite), score.TOCSY = annotationTOCSY$liste_resultat$score)
    dataTOCSY <- unique.data.frame(dataTOCSY)
  }

  seqCombiMeanScoreSeuil <- data.frame()
  seqCombiMeanScoreSeuilFiltre <- data.frame()

  ## CONCATENATION RESULTATS DIFFERENTES SEQUENCES
  data2D <- list(dataCOSY, dataHMBC, dataHSQC, dataJRES, dataTOCSY)
  whichSequenceNaN <- which((data2D != "NA"))
  data2D <- data2D[whichSequenceNaN]
  sequencesCombination <- data.frame(data2D[1])
  seqCombiMeanScore <- sequencesCombination

    ## Si une seule sequence et seuil sur score = filtre applique dans la fonction annotationRmn2D
  if (length(data2D) >= 2) {
    ## CONCATENATION SCORE PAR SEQUENCE
    for (l in 2:length(data2D))
        sequencesCombination <- merge.data.frame(sequencesCombination, data2D[l], by = "Metabolite", all.x = TRUE, all.y = TRUE)

  ## Replacement of NA values due to mis annotation
  for (m in 1:nrow(sequencesCombination)){
    COSYcompound <- sort(names(BdDReference_COSY))
    HMBCcompound <- sort(names(BdDReference_HMBC))
    HSQCcompound <- sort(names(BdDReference_HSQC))
    JREScompound <- sort(names(BdDReference_JRES))
    TOCSYcompound <- sort(names(BdDReference_TOCSY))
    
    if (is.na(sequencesCombination[m, 2])){
      compound <- as.character(sequencesCombination[m, 1])
      for (c in 1:length(COSYcompound))
        if (str_to_lower(compound) == str_to_lower(COSYcompound[c]))
          sequencesCombination[m, 2] <- 0
    }
    
    if (is.na(sequencesCombination[m, 3])){
      compound <- as.character(sequencesCombination[m, 1])
      for (c in 1:length(HMBCcompound))
        if (str_to_lower(compound) == str_to_lower(HMBCcompound[c]))
          sequencesCombination[m, 3] <- 0
    }
    
    if (is.na(sequencesCombination[m, 4])){
      compound <- as.character(sequencesCombination[m, 1])
      for (c in 1:length(HSQCcompound))
        if (str_to_lower(compound) == str_to_lower(HSQCcompound[c]))
          sequencesCombination[m, 4] <- 0
    }
    
    if (is.na(sequencesCombination[m, 5])){
      compound <- as.character(sequencesCombination[m, 1])
      for (c in 1:length(JREScompound))
        if (str_to_lower(compound) == str_to_lower(JREScompound[c]))
          sequencesCombination[m, 5] <- 0
    }
    
    if (is.na(sequencesCombination[m, 6])){
      compound <- as.character(sequencesCombination[m, 1])
      for (c in 1:length(TOCSYcompound))
        if (str_to_lower(compound) == str_to_lower(TOCSYcompound[c]))
          sequencesCombination[m, 6] <- 0
    }
  }

    ## SCORE MOYEN (sans prise en compte valeurs manquantes)
    meanScore <- round(apply(sequencesCombination[, -1], 1, FUN = mean.rmNa), 2)
    seqCombiMeanScore <- cbind.data.frame(sequencesCombination, averageScore = meanScore)

        ## SUPPRESSION METABOLITE AVEC SCORE MOYEN < SEUIL
    seqCombiMeanScoreSeuilFiltre <- seqCombiMeanScore[seqCombiMeanScore$averageScore > seuil, ]
  }

  return(list(COSY = annotationCOSY, HMBC = annotationHMBC, HSQC = annotationHSQC, JRES = annotationJRES, TOCSY = annotationTOCSY, combination = seqCombiMeanScoreSeuilFiltre))
}
