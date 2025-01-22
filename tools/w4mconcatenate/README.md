# W4M concatenate


Metadata
-----------

 * **@name**: W4M concatenate
 * **@galaxyID**: W4Mconcatenate
 * **@version**: 1.0.0
 * **@authors**: Original code: Hanane Nourine (Intern - PFEM - INRAE) - Maintainer: Melanie Petera (PFEM - INRAE - MetaboHUB)
 * **@init date**: 2024, May
 * **@main usage**: This tool enables the concatenation of two matrices of Metadata and returns a matrix of Metadata and two DataMatrix

 
Context
-----------

This tool is provided as one of the [Workflow4Metabolomics](http://workflow4metabolomics.org) Galaxy instance data handling tools. 
W4M is an international infrastructure providing software tools to process, analyse and annotate metabolomics data. 

User interface is based on the Galaxy platform (homepage: https://galaxyproject.org/). It is an open, web-based platform for data intensive biomedical research. 
Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses.

The tool was created as part of a Master's level internship.
 
Configuration
-----------

### Requirement:
 * R software: version = 4.3.3 recommended
 * r-dplyr = 1.1.4
 * r-w4mrutils = 1.0.0

Technical description
-----------

Main files:

- concatenation.R: R function (core script)
- fonctions_auxiliaires.R : R auxiliary functions.
- concatenation_wrapper.R: R script to link the main R function to inputs
- concat.xml: XML wrapper (interface for Galaxy)


Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10
                     


License
-----------

 * Cea Cnrs Inria Logiciel Libre License, version 2.1 (CECILL-2.1)
