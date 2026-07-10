# ImputeMSmv


Metadata
-----------

 * **@name**: ImputeMSmv
 * **@galaxyID**: imoutemsmv
 * **@version**: 1.0.0
 * **@authors**:
  - Original R code: Mathieu Cladière
  - Original Galaxy wrapper and R code adjustment: Anthony Sleiman
  - Maintainer: Melanie Petera (PFEM - INRAE - MetaboHUB)
 * **@init date**: 2026, April
 * **@main usage**: This tool is used for the classification and imputation of missing values in mass-spectrometry-based metabolomics data using group-aware MAR/MNAR strategies

 
Context
-----------

This Galaxy tool is provided as one of the [Workflow4Metabolomics](http://workflow4metabolomics.org) Galaxy instance data handling tools. 
W4M is an international infrastructure providing software tools to process, analyse and annotate metabolomics data. 

User interface is based on the Galaxy platform (homepage: https://galaxyproject.org/). It is an open, web-based platform for data intensive biomedical research. 
Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses.

The tool was created as part of a Master's level internship.
 
Configuration
-----------

### Requirement:
 * R software: version = 4.5.2 recommended
 * r-w4mrutils = 1.2.1

Technical description
-----------

Main files:

- ImputeMSmv_function.R: R function (core script)
- ImputeMSmv_wrapper.R: R script to link the main R function to inputs
- ImputeMSmv.xml: XML wrapper (interface for Galaxy)


Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10
                     


License
-----------

 * LGPL-3.0-or-later