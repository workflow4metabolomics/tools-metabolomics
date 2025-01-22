#################################################################################################################
# SPECTRA NORMALIZATION FROM SPECTRAL DATA                                                    #
# User : Galaxy                                                                                                 #
# Original data : --                                                                                            #
# Starting date : 20-10-2014                                                                                    #
# Version 1 : 27-01-2015                                                                                        #
# Version 2 : 27-02-2015                                                                                        #
#                                                                                                               #
# Input files:                                                                                                  #
#   - Data matrix containing bucketed and integrated spectra to normalize                                       #
#   - Sample metadata matrix containing at least biological factor of interest                                  #
#   - Scaling method: Total intensity/Probabilistic Quotient Normalization                                      #
#   - Control group: name of control to compute median reference spectra                                        #
#   - Graph: normalization result representation                                                                #
#################################################################################################################
NmrNormalization <- function(dataMatrix, scalingMethod = c("None", "Total", "PQN", "BioFactor"), sampleMetadata = NULL,
                             bioFactor = NULL, ControlGroup = NULL, graph = c("None", "Overlay", "One_per_individual"),
                             nomFichier = NULL, savLog.txtC = NULL) {
    ## Option
    ## ---------------
    strAsFacL <- options()$stringsAsFactors
    options(stingsAsFactors = FALSE)
    options(warn = -1)


    ## Constants
    ## ---------------
    topEnvC <- environment()
    flgC <- "\n"

    ## Log file (in case of integration into Galaxy)
    ## ----------------------------------------------
    if (!is.null(savLog.txtC)) {
        sink(savLog.txtC, append = TRUE)
    }

    ## Functions definition
    ## ---------------------
    #################################################################################################################
    # Total intensity normalization
    # Input parameters
    #   - data : bucketed spectra (rows=buckets; columns=samples)
    #################################################################################################################
    NmrBrucker_total <- function(data) {
        # Total intensity normalization
        data.total <- apply(data, 2, sum)
        data.normalized <- data[, 1] / data.total[1]
        for (i in 2:ncol(data)) {
            data.normalized <- cbind(data.normalized, data[, i] / data.total[i])
        }
        colnames(data.normalized) <- colnames(data)
        rownames(data.normalized) <- rownames(data)
        return(data.normalized)
    }


    #################################################################################################################
    # Biological factor normalization
    # Input parameters
    #   - data : bucketed spectra (rows=buckets; columns=samples)
    #   - sampleMetadata : dataframe with biological factor of interest measured for each invidual
    #   - bioFactor : name of the column cotaining the biological factor of interest
    #################################################################################################################
    NmrBrucker_bioFact <- function(data, sampleMetadata, bioFactor) {
        # Total intensity normalization
        data.normalized <- data[, 1] / bioFactor[1]
        for (i in 2:ncol(data)) {
            data.normalized <- cbind(data.normalized, data[, i] / bioFactor[i])
        }
        colnames(data.normalized) <- colnames(data)
        rownames(data.normalized) <- rownames(data)
        return(data.normalized)
    }


    #################################################################################################################
    # Probabilistic quotient normalization (PQN)
    # Input parameters
    #   - data : bucketed spectra (rows=buckets; columns=samples)
    #   - sampleMetadata : dataframe with treatment group of inviduals
    #   - pqnFactor : number of the column cotaining the biological facor of interest
    #   - nomControl : name of the treatment group
    #################################################################################################################
    NmrBrucker_pqn <- function(data, sampleMetadata, pqnFactor, nomControl) {
        # Total intensity normalization
        data.total <- apply(data, 2, sum)
        data.normalized <- data[, 1] / data.total[1]
        for (i in 2:ncol(data)) {
            data.normalized <- cbind(data.normalized, data[, i] / data.total[i])
        }
        colnames(data.normalized) <- colnames(data)
        rownames(data.normalized) <- rownames(data)

        # Reference spectrum
        # Recuperation spectres individus controle
        control.spectra <- data.normalized[, sampleMetadata[, pqnFactor] == nomControl]
        spectrum.ref <- apply(control.spectra, 1, median)
        for (j in 1:length(spectrum.ref))
        {
            if (spectrum.ref[j] == 0) {
                spectrum.ref[j] <- mean(control.spectra[j, ])
            }
            if (spectrum.ref[j] == 0) {
                spectrum.ref[j] <- 10^(-24)
            }
        }

        # Ratio between normalized and reference spectra
        data.normalized.ref <- data.normalized / spectrum.ref

        # Median ratio
        data.normalized.ref.median <- apply(data.normalized.ref, 1, median)
        for (j in 1:length(data.normalized.ref.median)) {
            if (data.normalized.ref.median[j] == 0 | is.na(data.normalized.ref.median[j]) | data.normalized.ref.median == "NaN" | data.normalized.ref.median == "NA") {
                data.normalized.ref.median[j] <- mean(data.normalized.ref[j, ])
            }
        }

        # Normalization
        data.normalizedPQN <- data.normalized[, 1] / data.normalized.ref.median
        for (i in 2:ncol(data)) {
            data.normalizedPQN <- cbind(data.normalizedPQN, data.normalized[, i] / data.normalized.ref.median)
        }
        colnames(data.normalizedPQN) <- colnames(data)
        rownames(data.normalizedPQN) <- rownames(data)

        return(data.normalizedPQN)
    }


    ## Tests
    if (scalingMethod == "QuantitativeVariable") {
        if (mode(sampleMetadata[, bioFactor]) == "character") {
            bioFact <- factor(sampleMetadata[, bioFactor])
        } else {
            bioFact <- sampleMetadata[, bioFactor]
        }
    }

    ## Spectra scaling depending on the user choice
    if (scalingMethod == "None") {
        NormalizedBucketedSpectra <- dataMatrix
    } else if (scalingMethod == "Total") {
        NormalizedBucketedSpectra <- NmrBrucker_total(dataMatrix)
    } else if (scalingMethod == "PQN") {
        NormalizedBucketedSpectra <- NmrBrucker_pqn(dataMatrix, sampleMetadata, bioFactor, ControlGroup)
    } else if (scalingMethod == "QuantitativeVariable") {
        NormalizedBucketedSpectra <- NmrBrucker_bioFact(dataMatrix, sampleMetadata, bioFact)
    }

    ## OUTPUTS
    return(list(NormalizedBucketedSpectra))
}
