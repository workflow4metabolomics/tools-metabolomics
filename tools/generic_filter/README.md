Generic Filter
=======

Metadata
-----------

 * **@name**: Generic Filter
 * **@version**: 1.2.2
 * **@authors**: Marion Landi and Melanie Petera for first version - Maintainer: Melanie Petera (PFEM ; INRAE ; MetaboHUB)
 * **@init date**: 2014, december
 * **@main usage**: This tool allows to remove all samples and/or variables corresponding to specific values regarding designated factors or numerical variables. 

 
Context
-----------

This tool is provided as one of the [Workflow4Metabolomics](http://workflow4metabolomics.org) Galaxy instance utilities. W4M is a French infrastructure providing software tools to process, analyse and annotate metabolomics data. 

User interface is based on the Galaxy platform (homepage: https://galaxyproject.org/). It is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses.

 
Configuration
-----------

### Requirement:
 * R software: version > 3.0.0 recommended
 * Specific R library: 'W4MRUtils'


Technical description
-----------

Main files:

- filter_wrap.R: R script to format parameters and call the main R function
- generic_filter.xml: XML wrapper (interface for Galaxy)


Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10


License
-----------

 * Cea Cnrs Inria Logiciel Libre License, version 2.1 (CECILL-2.1)
