Metrics and graphics to assess the quality of the data
======================================================

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure  

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/qualitymetrics.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/qualitymetrics).

### Description

**Version:** 2.2.8  
**Date:** 2018-01-11  
**Author:** Marion Landi (INRA, PFEM), Mélanie Pétéra (INRA, PFEM), and Etienne A. Thévenot (CEA, LIST)  
**Email:** [melanie.petera(at)clermont.inra.fr](mailto:melanie.petera@clermont.inra.fr), [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:**  
**Licence:** CeCILL  
**Reference history:** [W4M00001b_sacurine-complete](http://galaxy.workflow4metabolomics.org/history/list_published)   
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)  

### Installation

* Configuration file:
    + `qualitymetrics_config.xml`  
* Image files: 
    + `static/images/QualityControl.png`    
    + `static/images/qualitymetrics_workingExampleImage.png`      
* Wrapper file:
    + `qualitymetrics_wrapper.R`  
* Script file:
    + `qualitymetrics_script.R`  
* R packages
  + **batch** from CRAN  
  
    ```r
    install.packages("batch", dep=TRUE)  
    ```

  + **ropls** from Bioconductor  
  
    ```r
    source("http://www.bioconductor.org/biocLite.R")  
    biocLite("ropls")      
    ```

### Tests

The code in the wrapper can be tested by running the `runit/qualitymetrics_runtests.R` R file

You will need to install **RUnit** package in order to make it run:
```r
install.packages('RUnit', dependencies = TRUE)
```

### News  

##### CHANGES IN VERSION 2.2.8  

MINOR MODIFICATION  

In the case of a distinct sample order between dataMatrix and sampleMetadata, the sample order from the dataMatrix is matched to sampleMetadata internally for the computations and graphics without modifying the order in the sampleMetadata output (a warning is generated in the information file); to get the re-ordered dataMatrix as output, please use the Check Format module  

##### CHANGES IN VERSION 2.2.6  

MINOR MODIFICATION  

 * Graphic (pool_CV inferior to 30% metric): pools with a NaN value are now counted as having a value superior to 30% (to avoid generating a final NA metric value)  

##### CHANGES IN VERSION 2.2.4  

INTERNAL MODIFICATION    

 * Additional running and installation tests added with planemo, conda, and travis  

##### CHANGES IN VERSION 2.2.3  

INTERNAL MODIFICATION    

 * Modifications of the 'qualitymetrics_script.R' file to handle the recent 'ropls' package versions (i.e. 1.3.15 and above) which use S4 classes  

 * Creating tests for the R code  
    
##### CHANGES IN VERSION 2.2.2

INTERNAL MODIFICATION  

 * Minor internal modification
