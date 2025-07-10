Multivariate analysis by PCA, PLS(-DA), and OPLS(-DA)
=====================================================

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure  

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/multivariate.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/multivariate).

### Description

**Version:** 2.3.14  
**Date:** 2018-08-16  
**Author:** Etienne A. Thevenot (CEA, LIST, MetaboHUB, W4M Core Development Team)   
**Email:** [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research*, **14**:3322-3335. [doi:10.1021/acs.jproteome.5b00354](http://dx.doi.org/10.1021/acs.jproteome.5b00354)  
**Licence:** CeCILL  
**Reference history:** [W4M00001a_sacurine-subset-statistics](http://galaxy.workflow4metabolomics.org/history/list_published), [W4M00001b_sacurine_complete](http://galaxy.workflow4metabolomics.org/history/list_published)     
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)

### Installation

* Configuration file: `multivariate_config.xml` 
* Image files: 
  + `static/images/multivariate_workflowPositionImage.png` 
  + `static/images/multivariate_workingExampleImage.png` 
* Wrapper file: `multivariate_wrapper.R` 
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

The code in the wrapper can be tested by running the `runit/multivariate_runtests.R` R file

You will need to install **RUnit** package in order to make it run:
```r
install.packages('RUnit', dependencies = TRUE)
```

### Working example  

See the **W4M00001a_sacurine-subset-statistics**, **W4M00001b_sacurine-subset-complete**, **W4M00002_mtbls2**, **W4M00003_diaplasma** shared histories in the **Shared Data/Published Histories** menu (https://galaxy.workflow4metabolomics.org/history/list_published) 

### News  

###### CHANGES IN VERSION 2.3.14  

MINOR CORRECTION  

 * Score plot (x-score) advanced parameters (labels and colors) in case of test samples in the dataset  

###### CHANGES IN VERSION 2.3.12  

NEW FEATURE 

 * Some of the parameters were renamed according to the controlled vocabulary of the Ontology for Biomedical Investigations (OBI; http://obi-ontology.org/) and the statistics ontology (STATO; http://stato-ontology.org/), in collaboration with Philippe Rocca-Serra (University of Oxford)  

###### CHANGES IN VERSION 2.3.8   

MINOR MODIFICATION  

 * (O)PLS(-DA) coefficients display in case of multiple quantitative (or multiclass) response: now the column names of the coefficients for each response are correctly labelled in the variableMetadata file  
 
###### CHANGES IN VERSION 2.3.6  

INTERNAL MODIFICATION  

 * minor internal modifications  

###### CHANGES IN VERSION 2.3.4

INTERNAL MODIFICATION  

 * minor update of .shed.yml

###### CHANGES IN VERSION 2.3.2

INTERNAL MODIFICATION  

 * Modification of the `multivariate_wrapper.R` file to handle **ropls** package versions of 1.3.15 and above (i.e. after switching to S4 classes)

###### CHANGES IN VERSION 2.3.0  

NEW FEATURES  

 * **Predictions** now available (see the 'Samples to be tested' argument)
 * OPLS(-DA): **Predictive and Orthogonal VIP** are now computed (see the 'comments' section)
 * **Multiclass PLS-DA** implemented (see the 'comments' section)
