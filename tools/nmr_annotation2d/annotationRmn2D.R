##########################################################################
# ANNOTATION SPECTRE 2D MATRICE COMPLEXE BASEE SUR UNE SEQUENCE RMN      #
# matriceComplexe : data.frame liste couples ppm de la matrice a annoter #
# BdDStandards : objet contenant la base de donnees des composes standards #
# nom_sequence : nom sequence 2D a utiliser pour annotation ("JRES", "COSY", "TOCSY", "HMBC", "HSQC") #
# ppm1Tol : tolerance ppm axe abscisses                                                               #
# ppm2Tol : tolerance ppm axe ordonnees                                                               #
# nb_ligne_template : preciser le nombre total de ligne de la feuille de calcul a annoter             #
#######################################################################################################
annotationRmn2D <- function(matriceComplexe, BdDStandards, nom_sequence, ppm1Tol = 0.01, ppm2Tol = 0.01,
                            seuil = 0, unicite = "NO") {
  ## Longueur de la peak-list de la matrice a annoter
  PeakListLength <- length(matriceComplexe[, 1])

  ## Nombre de metabolites inclus dans BdD de composes standards
  nbMetabolitesBdD <- length(BdDStandards)
  matrixAnnotation <- data.frame()
  allMetabolitesList <- data.frame()
  seuil_score <- seuil

  ## Boucle sur les metabolites inclus dans BdD
  for (i in seq_len(nbMetabolitesBdD)) {
    ## Infos metabolite en cours
    iMetabolite <- BdDStandards[[i]]
    ppm1M <- iMetabolite[, 1]
    ppm2M <- iMetabolite[, 2]
    nbPeakMetabolite <- length(ppm1M)
    MetaboliteName <- names(BdDStandards[i])

    ## Initialisation
    k <- 0
    presenceScore <- 0
    annotatedPpmRef <- data.frame()
    annotatedPpmList <- data.frame()
    annotatedPeakLength <- 0
    metabolites <- data.frame()
    metabolitesList <- data.frame()

    ## Boucle sur les couples de pics de la matrice a annoter
    for (p in seq_len(PeakListLength)) {
      ppmAnnotationF1 <- as.numeric(matriceComplexe[p, 3])
      ppmAnnotationF2 <- as.numeric(matriceComplexe[p, 2])
      e <- simpleMessage("end of file")
      tryCatch({
        if (!is.na(ppmAnnotationF1)) {
          matrixAnnotation <- unique.data.frame(rbind.data.frame(matrixAnnotation, matriceComplexe[p, ]))
        }
        # Recherche du couple de pics de la matrice la liste des couples du metabolite standard
        metaboliteIn <- (ppm1M >= (ppmAnnotationF2 - ppm1Tol) & ppm1M <= (ppmAnnotationF2 + ppm1Tol) &
                         ppm2M >= (ppmAnnotationF1 - ppm2Tol) & ppm2M <= (ppmAnnotationF1 + ppm2Tol))
        WhichMetaboliteIn <- which(metaboliteIn)
        # Si au moins un couple de la matrice a annoter dans liste couples metabolite standard
        if (length(WhichMetaboliteIn) > 0) {
          for (a in seq_len(length(WhichMetaboliteIn))) {
            annotatedPpmList <- data.frame(ppm1 = ppm1M[WhichMetaboliteIn[a]], ppm2 = ppm2M[WhichMetaboliteIn[a]], theoricalLength = nbPeakMetabolite)
            annotatedPpmRef <- rbind(annotatedPpmRef, annotatedPpmList)
          }
        }
      }, error = function(e) {
        cat("End of file \n")
      })
    }

    # Au - 1 couple de ppm de la matrice complexe annote
    if (nrow(annotatedPpmRef) >= 1) {
      ## Nombre couples annotes
      annotatedPeakLength <- nrow(annotatedPpmRef)

      ## Recherche doublons
      annotatedDoublons <- duplicated(annotatedPpmRef)
      if (sum(duplicated(annotatedPpmRef)) > 0) {
        annotatedPeakLength <- nrow(annotatedPpmRef) - sum(duplicated(annotatedPpmRef))
        annotatedPpmRef <- annotatedPpmRef[-duplicated(annotatedPpmRef), ]
      }
      presenceScore <- round(annotatedPeakLength / nbPeakMetabolite, 2)
    }

    ## Conservation metabolites dont score > seuil
    if (presenceScore > seuil_score) {
      metabolites <- data.frame(Metabolite = MetaboliteName, score = presenceScore)
      metabolitesList <- cbind.data.frame(annotatedPpmRef, metabolites)
      allMetabolitesList <- rbind.data.frame(allMetabolitesList, metabolitesList)
    }
  }

  # Initialisation
  commonPpm <- data.frame()
  commonPpmList <- data.frame()
  metaboliteAdd <- data.frame()
  metaboliteAddList <- data.frame()
  commonMetabolitesList <- data.frame()
  commonMetabolitesPpmList <- data.frame()
  commonMetabolitesPpmAllList1 <- data.frame()
  commonMetabolitesPpmAllList <- data.frame()
  listeTotale_2D_unicite <- allMetabolitesList[, 1:4]
  allMetabolitesList <- allMetabolitesList[, -3]
  metabolitesAllUnicite <- data.frame()

  ## Boucle sur tous couples annotes
  for (j in seq_len(length(allMetabolitesList$ppm1))) {
    ## Boucle sur metabolites dans BdD composes standards
    for (i in seq_len(nbMetabolitesBdD)) {
      ppmMetaboliteBdD <- BdDStandards[[i]]
      ppm1M <- ppmMetaboliteBdD[, 1]
      ppm2M <- ppmMetaboliteBdD[, 2]
      # Nombre de couples metabolite
      nbPeakMetabolite <- length(ppm1M)
      MetaboliteName <- names(BdDStandards[i])

      metabolitesInAll <- (ppm1M >= (allMetabolitesList[j, 1] - ppm1Tol) & ppm1M <= (allMetabolitesList[j, 1] + ppm1Tol) &
                           ppm2M >= (allMetabolitesList[j, 2] - ppm2Tol) & ppm2M <= (allMetabolitesList[j, 2] + ppm2Tol))
      WhichMetabolitesInAll <- which(metabolitesInAll)

      if (MetaboliteName != allMetabolitesList[j, 3] & length(WhichMetabolitesInAll) > 0) {
        metabolitesAllUnicite <- rbind.data.frame(metabolitesAllUnicite, listeTotale_2D_unicite[j, ])
        commonPpm <- data.frame(ppm1 = allMetabolitesList[j, 1], ppm2 = allMetabolitesList[j, 2])
        commonPpmList <- rbind.data.frame(commonPpmList, commonPpm)
        commonPpmList <- unique(commonPpmList)
        metaboliteAdd <- data.frame(nom_metabolite = MetaboliteName)
        metaboliteAddList <- rbind.data.frame(metaboliteAddList, metaboliteAdd)
        commonMetabolitesList <- rbind.data.frame(data.frame(nom_metabolite = allMetabolitesList[j, 3]), metaboliteAddList)
        commonMetabolitesPpmList <- cbind.data.frame(commonPpm, commonMetabolitesList)
        commonMetabolitesPpmAllList1 <- rbind.data.frame(commonMetabolitesPpmAllList1, commonMetabolitesPpmList)
        commonMetabolitesPpmAllList1 <- unique.data.frame(commonMetabolitesPpmAllList1)
      }
    }
    commonMetabolitesPpmAllList <- rbind.data.frame(commonMetabolitesPpmAllList, commonMetabolitesPpmAllList1)
    commonMetabolitesPpmAllList <- unique.data.frame(commonMetabolitesPpmAllList)

    #initialisation des data.frame
    commonPpm <- data.frame()
    metaboliteAdd <- data.frame()
    metaboliteAddList <- data.frame()
    metabolite_ref <- data.frame()
    commonMetabolitesList <- data.frame()
    commonMetabolitesPpmList <- data.frame()
    commonMetabolitesPpmAllList1 <- data.frame()
  }

  unicityAllList <- listeTotale_2D_unicite
  if (nrow(listeTotale_2D_unicite) != 0 & nrow(metabolitesAllUnicite) != 0)
    unicityAllList <- setdiff(listeTotale_2D_unicite, metabolitesAllUnicite)

  unicitynbCouplesRectif <- data.frame()
  for (g in seq_len(nrow(unicityAllList))) {
    metaboliteUnicity <- (unicityAllList$Metabolite == unicityAllList$Metabolite[g])
    WhichMetaboliteUnicity <- which(metaboliteUnicity)
    nb_occurence <- length(WhichMetaboliteUnicity)
    unicitynbCouplesRectif <- rbind.data.frame(unicitynbCouplesRectif, nb_occurence)
  }
  names(unicitynbCouplesRectif) <- "NbCouplesAnnotes"
  unicityAllList <- cbind.data.frame(unicityAllList, unicitynbCouplesRectif)

  unicityAllList <- cbind.data.frame(unicityAllList, score_unicite = unicityAllList$NbCouplesAnnotes / unicityAllList$theoricalLength)
  unicityAllList <- unicityAllList[, -3]
  unicityAllList <- unicityAllList[, -4]

  unicityAllList <- unicityAllList[unicityAllList$score_unicite > seuil_score, ]

  listeTotale_metabo <- data.frame()
  if (nrow(commonPpmList) != 0) {
    for (o in seq_len(length(commonPpmList[, 1]))) {
      tf6 <- (commonMetabolitesPpmAllList$ppm1 == commonPpmList[o, 1] & commonMetabolitesPpmAllList$ppm2 == commonPpmList[o, 2])
      w6 <- which(tf6)

      for (s in seq_len(length(w6))) {
        metaboliteAdd <- data.frame(nom_metabolite = commonMetabolitesPpmAllList[w6[s], 3])
        commonMetabolitesList <- paste(commonMetabolitesList, metaboliteAdd[1, ], sep = " ")
      }
      liste_metabo_ppm <- cbind.data.frame(ppm1 = commonPpmList[o, 1], ppm2 = commonPpmList[o, 2], commonMetabolitesList)
      listeTotale_metabo <- rbind.data.frame(listeTotale_metabo, liste_metabo_ppm)
      commonMetabolitesList <- data.frame()
    }
  }

  # Representation graphique
  if (nom_sequence == "HSQC" | nom_sequence == "HMBC") {
    atome <- "13C"
    indice_positif <- 1
    indice_negatif <- -10
  } else {
    atome <- "1H"
    indice_positif <- 0.5
    indice_negatif <- -0.5
  }

  matriceComplexe <- matrixAnnotation
  ppm1 <- as.numeric(matriceComplexe[, 2])
  ppm2 <- as.numeric(matriceComplexe[, 3])

  if (unicite == "NO") {
    listeTotale_2D_a_utiliser <- allMetabolitesList
    d1.ppm <- allMetabolitesList$ppm1
    d2.ppm <- allMetabolitesList$ppm2
  } else {
    listeTotale_2D_a_utiliser <- unicityAllList
    d1.ppm <- listeTotale_2D_a_utiliser$ppm1
    d2.ppm <- listeTotale_2D_a_utiliser$ppm2
  }

  if (nrow(listeTotale_2D_a_utiliser) > 0) {
    ## Taches de correlations
    # Matrice biologique + Annotations
    maxX <- max(round(max(as.numeric(matriceComplexe[, 2]))) + 0.5, round(max(as.numeric(matriceComplexe[, 2]))))
    maxY <- max(round(max(as.numeric(matriceComplexe[, 3]))) + indice_positif, round(max(as.numeric(matriceComplexe[, 3]))))
    probability.score <- as.factor(round(listeTotale_2D_a_utiliser[, 4], 2))
    lgr <- length(unique(probability.score))
    sp <- ggplot(matriceComplexe, aes(x = ppm1, y = ppm2))
    sp <- sp + geom_point(size = 2) + scale_x_reverse(breaks = seq(maxX, 0, -0.5)) +
      scale_y_reverse(breaks = seq(maxY, 0, indice_negatif)) +
      xlab("1H chemical shift (ppm)") + ylab(paste(atome, " chemical shift (ppm)")) + ggtitle(nom_sequence) +
      geom_text(data = listeTotale_2D_a_utiliser, aes(d1.ppm, d2.ppm, label = str_to_lower(substr(listeTotale_2D_a_utiliser[, 3], 1, 3)), col = probability.score),
                size = 4, hjust = 0, nudge_x = 0.02, vjust = 0, nudge_y = 0.2) + scale_colour_manual(values = viridis(lgr))
    print(sp)
  }

  # Liste des resultats (couples pmm / metabolite / score) + liste ppms metabolites communs
  if (unicite == "NO") {
    return(list(liste_resultat = allMetabolitesList, listing_ppm_commun = listeTotale_metabo))
  } else {
    return(list(liste_resultat_unicite = unicityAllList, listing_ppm_commun_affichage = listeTotale_metabo))
  }
}
