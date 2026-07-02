Univariate parametric and non-parametric hypothesis testing with correction for multiple testing  
================================================================================================  

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure    

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/univariate.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/univariate).

### Description

**Version:** 2.2.4   
**Date:** 2018-01-11  
**Author:** Marie Tremblay-Franco (INRA, MetaToul, MetaboHUB, W4M Core Development Team) and Etienne A. Thevenot (CEA, LIST, MetaboHUB, W4M Core Development Team)    
**Email:** [marie.tremblay-franco(at)toulouse.inra.fr](mailto:marie.tremblay-franco@toulouse.inra.fr); [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research*, **14**:3322-3335. [doi:10.1021/acs.jproteome.5b00354](http://dx.doi.org/10.1021/acs.jproteome.5b00354)  
**Reference history:** [W4M00001a_sacurine-subset-statistics](http://galaxy.workflow4metabolomics.org/history/list_published), [W4M00004_mtbls1](http://galaxy.workflow4metabolomics.org/history/list_published)  
**Licence:** CeCILL  
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)

### Installation

 * Configuration file: `univariate_config.xml`  
 * Image file: 
  + `static/images/univariate_workflowPositionImage.png`   
 * Wrapper file: `univariate_wrapper.R`  
 * Script file: `univariate_script.R`  
 * R packages  
  + **batch** from CRAN  
  
    ```r
    install.packages("batch", dep=TRUE)  
    ```
  + **PMCMR** from CRAN  
  
    ```r  
    install.packages("PMCMR", dep=TRUE)  
    ```    
 
### Tests

The code in the wrapper can be tested by running the `runit/univariate_runtests.R` R file

You will need to install **RUnit** package in order to make it run:
```r
install.packages('RUnit', dependencies = TRUE)
```

### Working example  

See the **W4M00001a_sacurine-subset-statistics**, **W4M00001b_sacurine-subset-complete**, **W4M00002_mtbls2**, **W4M00003_diaplasma** shared histories in the **Shared Data/Published Histories** menu (https://galaxy.workflow4metabolomics.org/history/list_published)  

### News

###### CHANGES IN VERSION 2.2.4  

MINOR MODIFICATION  

 * Internal minor modifications for building and testing  

###### CHANGES IN VERSION 2.2.0  

MAJOR MODIFICATION  

 * ANOVA and Kruskal-Wallis: The p-values of the post-hoc tests (i.e. from pairwise comparisons) are now further corrected for multiple testing over all variables (previously, only the p-value of the -first- omnibus test was corrected over all variables)  

MINOR MODIFICATION  

 * All values in the 'dif', adjusted p-value, and 'sig' columns are now displayed (previously, the values were set to NA when the p-value of the omnibus test was not significant)  

NEW FEATURE  

 * Graphic: a single pdf file containing the graphics of all significant tests is now produced as '_figure.pdf' output: boxplots (respectively scatterplots with the regression line in red and the R2 value) are displayed when the factor of interest is qualitative (respectively quantitative). The corrected p-value is indicated in the title of each plot  

###### CHANGES IN VERSION 2.1.4

NEW FEATURE  

 * Level names are now separated by '.' instead of '-' previously in the column names of the output variableMetadata table (e.g., 'jour_ttest_J3.J10_fdr' instead of 'jour_ttest_J3-J10_fdr' previously)  

INTERNAL MODIFICATION  

 * Minor internal changes  

###### CHANGES IN VERSION 2.1.2  

INTERNAL MODIFICATION  

 * Minor internal changes in .shed.yml for toolshed export

###### CHANGES IN VERSION 2.1.1  

INTERNAL MODIFICATION  

 * Internal handling of 'NA' p-values (e.g. when intensities are identical in all samples).
