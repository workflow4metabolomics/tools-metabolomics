Heatmap of the data matrix  
==========================

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure  

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/heatmap.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/heatmap).

### Description

**Version:** 2.2.2    
**Date:** 2016-08-04    
**Author:** Etienne A. Thevenot (CEA, LIST, MetaboHUB, W4M Core Development Team)   
**Email:** [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research*, **14**:3322-3335. [doi:10.1021/acs.jproteome.5b00354](http://dx.doi.org/10.1021/acs.jproteome.5b00354)  
**Licence:** CeCILL  
**Reference history:** [W4M00001a_sacurine-subset-statistics](http://galaxy.workflow4metabolomics.org/history/list_published)     
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)

### Installation

* Configuration file: `heatmap_config.xml`
* Image files:
    + `static/images/heatmap_workflowPositionImage.png`  
    + `static/images/heatmap_workingExampleImage.png`     
* Wrapper file: `heatmap_wrapper.R`  
* R packages  
    + **batch** from CRAN
  
    ```r
    install.packages("batch", dep=TRUE)  
    ```

### Tests

The code in the wrapper can be tested by running the `runit/heatmap_runtests.R` R file

### Working example  

See the **W4M00001a_sacurine-subset-statistics** shared history in the **Shared Data/Published Histories** menu (https://galaxy.workflow4metabolomics.org/history/list_published)  

### News

##### CHANGES IN VERSION 2.2.2  

INTERNAL MODIFICATION  

 * Minor internal modifications  

##### CHANGES IN VERSION 2.2.0  

NEW FEATURES  

 * Default method for the correlation coefficient is now 'pearson', instead of 'spearman' previously (the latter can still be selected in the advanced parameters)

 * The 1-abs(correlation) dissimilarity is now available (in addition to the default '1-correlation') in case the sign of correlations between samples and between variables does not matter, as well as the euclidean, maximum, manhattan, canberra, binary, and minkowski dissimilarities

 * A new red-green color scale is available

##### CHANGES IN VERSION 2.1.2  

INTERNAL MODIFICATION  

 * Creating additional files for planemo and travis validation  

##### CHANGES IN VERSION 2.0.2  

INTERNAL MODIFICATION  

 * Creating tests for R code  