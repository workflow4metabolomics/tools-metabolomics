Preprocessing workflow for Flow Injection Analysis coupled to High-Resolution Mass Spectrometry data (FIA-HRMS)
===============================================================================================================

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure  

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/profia.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/profia).

### Description

**Version:** 3.1.0   
**Date:** 2018-01-08     
**Author:** Alexis Delabriere and Etienne A. Thevenot (CEA, LIST, MetaboHUB, W4M Core Development Team)   
**Email:** [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Delabriere A., Hohenester U., Colsch B., Junot C., Fenaille F. and Thevenot E.A. (2017). proFIA: A data preprocessing workflow for Flow Injection Analysis coupled to High-Resolution Mass Spectrometry. *Bioinformatics*, 33:3767-3775 (https://doi.org/10.1093/bioinformatics/btx458).   
**Licence:** CeCILL  
**Reference history:**      
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)

### Installation

* Configuration file: `profia_config.xml`
* Image files: 
  + `static/images/profia_workflowPositionImage.png`   
  + `static/images/profia_workingExampleImage.png`
* Wrapper file: `profia_wrapper.R`
* R packages
 + **batch** from CRAN  
  
    ```r
    install.packages("batch", dep=TRUE)
    install.packages("minpack.lm", dep=TRUE)
    install.packages("missForest", dep=TRUE)  
    install.packages("pracma", dep=TRUE)
    ```
  + **profia** from Bioconductor  
  
    ```r
    source("http://www.bioconductor.org/biocLite.R")
    biocLite("Biobase")
    biocLite("BiocParallel")
    biocLite("xcms")
    biocLite("plasFIA")
    biocLite("proFIA")
    ```  

### Tests

The code in the wrapper can be tested by running the `runit/profia_runtests.R` R file

You will need to install **RUnit** package in order to make it run:
```r
install.packages('RUnit', dependencies = TRUE)
```

### Working example  

### News

###### CHANGES IN VERSION 3.1.0  

NEW FEATURE  

randomForest method implemented for imputation of missing values  

###### CHANGES IN VERSION 3.0.6  

NEW FEATURES  

 * New (advanced) parameters available  

###### CHANGES IN VERSION 3.0.4  

MINOR MODIFICATION  

 * Details added in the documentation  

###### CHANGES IN VERSION 3.0.2  

NEW FEATURE  

 * Parallel processing  

###### CHANGES IN VERSION 3.0.0  

NEW FEATURE  

 * Creation of the tool  
