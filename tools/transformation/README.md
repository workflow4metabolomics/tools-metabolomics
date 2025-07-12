Transformation of the data matrix  
=================================

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure  

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/transformation.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/transformation).

### Description

**Version:** 2.2.0  
**Date:** 2016-08-04  
**Author:** Etienne A. Thevenot (CEA, LIST, MetaboHUB, W4M Core Development Team)   
**Email:** [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research*, **14**:3322-3335. [doi:10.1021/acs.jproteome.5b00354](http://dx.doi.org/10.1021/acs.jproteome.5b00354)  
**Licence:** CeCILL
**Reference history:** [W4M00001b_sacurine_complete](http://galaxy.workflow4metabolomics.org/history/list_published)     
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)

### Installation

* Configuration file: `transformation_config.xml`  
* Image files: `static/images/transformation_workflowPositionImage.png`   
* Wrapper file: `transformation_wrapper.R`  
* R packages  
    + **batch** from CRAN  
 
    ```r
    install.packages("batch", dep=TRUE)  
    ```
    
### Tests

The code in the wrapper can be tested by running the `runit/transformation_runtests.R` R file

You will need to install **RUnit** package in order to make it run:
```r
install.packages('RUnit', dependencies = TRUE)
```

### Working example

See the **W4M00001b_sacurine-complete** or **W4M00002_mtbls2** shared histories in the **Shared Data/Published Histories** menu (https://galaxy.workflow4metabolomics.org/history/list_published)  

### News

##### CHANGES IN VERSION 2.2.0  

NEW FEATURE  

 * Square root transformation now available  

##### CHANGES IN VERSION 2.0.2  

INTERNAL MODIFICATION  

 * Creating tests for R code  