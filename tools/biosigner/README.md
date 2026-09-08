Discovery of significant signatures from omics data
===================================================

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure  

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/biosigner.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/biosigner).

### Description

**Version:** 3.0.0  
**Date:** 2026-06-30     
**Author:** Philippe Rinaudo and Etienne A. Thevenot (CEA, LIST, MetaboHUB, W4M Core Development Team)   
**Email:** [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Rinaudo P., Boudah S., Junot C. and Thevenot E.A. (2016). *biosigner*: a new method for the discovery of significant molecular signatures from omics data. *Frontiers in Molecular Biosciences*, **3** (http://dx.doi.org/10.3389/fmolb.2016.00026).   
**Licence:** CeCILL
**Reference history:** [W4M00003_diaplasma](http://galaxy.workflow4metabolomics.org/history/list_published)      
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)

### Important Notice - Version 3.0.0 Breaking Changes

This version (3.0.0) includes significant updates to comply with modern Galaxy best practices and the latest bioconductor-biosigner package (v1.38.0):

1. **Parameter Renaming**: All parameters now use lowercase with underscores instead of camelCase. Existing workflows will need to be updated.

2. **Updated Function API**: The wrapper has been updated to use the correct function signatures for biosigner v1.38.0:
   - `biosign()` no longer accepts `printL`, `plotL`, or `.sinkC` arguments
   - `plot()` method uses `typeC` parameter ("tier" or "boxplot") instead of `tierMaxC`

See the tool's NEWS section below for the complete list of parameter name changes.

### Installation

* Configuration file: `biosigner_config.xml`
* Image files: 
  + `static/images/biosigner_workflowPositionImage.png`   
  + `static/images/biosigner_workingExampleImage.png`
* Wrapper file: `biosigner_wrapper.R`
* R packages (installed via Galaxy's conda dependencies):
  + **r-batch** (version 1.1_5+) from CRAN  
  + **bioconductor-biosigner** (version 1.38.0+) from Bioconductor  
  + **ropls** (version 1.4.0+) from Bioconductor

For manual installation outside Galaxy:
```r
# Install batch from CRAN
install.packages("batch", dep=TRUE)

# Install biosigner from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biosigner")
BiocManager::install("ropls")
```

### Tests

The code in the wrapper can be tested by running the `runit/biosigner_runtests.R` R file

You will need to install **RUnit** package in order to make it run:
```r
install.packages('RUnit', dependencies = TRUE)
```

### Working example  

See the **W4M00001a_sacurine-subset-statistics** or **W4M00003_diaplasma** shared histories in the **Shared Data/Published Histories** menu (https://galaxy.workflow4metabolomics.org/history/list_published)  

### News

###### CHANGES IN VERSION 3.0.0

BREAKING CHANGES

 * Updated Galaxy profile to 25.1 for modern Galaxy compatibility
 * Renamed all parameters to follow Galaxy best practices (lowercase with underscores):
   - `opcC` → `computational_mode`
   - `methodC` → `classification_method`
   - `bootI` → `num_bootstraps`
   - `tierC` → `selection_tier`
   - `pvalN` → `pvalue_threshold`
   - `seedI` → `random_seed`
   - `respC` → `response_column`
   - `dataMatrix_in` → `data_matrix`
   - `sampleMetadata_in` → `sample_metadata`
   - `variableMetadata_in` → `variable_metadata`
   - `variableMetadata_out` → `variable_metadata_out`
   - `figure_tier` → `tier_figure`
   - `figure_boxplot` → `boxplot_figure`
   - `information` → `log_information`
 * Added `detect_errors="exit_code"` to command element for better error handling
 * Updated wrapper to use correct biosigner v1.38.0 API:
   - Removed deprecated `printL`, `plotL`, `.sinkC` arguments from `biosign()` call
   - Updated `plot()` calls to use `typeC` parameter instead of `tierMaxC`

INTERNAL MODIFICATIONS

 * Tool now targets Galaxy profile 25.1
 * Updated R package requirements to biosigner v1.38.0 and r-batch v1.1_5
     
###### CHANGES IN VERSION 2.2.6  

INTERNAL MODIFICATION  

 * Minor internal modifications
 
###### CHANGES IN VERSION 2.2.4  

INTERNAL MODIFICATION  

 * Internal updates for planemo and travis validation  
 
###### CHANGES IN VERSION 2.2.2  

INTERNAL MODIFICATION  

 * Internal updates to biosigner package versions of 1.0.0 and above, and ropls versions of 1.4.0 and above (i.e. using S4 methods instead of S3)
    
###### CHANGES IN VERSION 2.2.1

NEW FEATURE

 * Creation of the tool: with S3 versions of biosigner (< 1.0.0) and ropls (< 1.4.0)
