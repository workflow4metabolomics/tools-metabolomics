Intensity check
=======

Metadata
-----------

 * **@name**: Intensity check
 * **@version**: 2.0.0
 * **@authors**: Original code: Anthony Fernandes (PFEM - INRA) - Maintainer: Melanie Petera (PFEM - INRAE - MetaboHUB)
 * **@init date**: 2018, September
 * **@main usage**: This tool computes some statistical measures, number of missing values and mean fold change.

 
Context
-----------

This tool is provided as one of the [Workflow4Metabolomics](http://workflow4metabolomics.org) Galaxy instance statistical tools. 
W4M is a French infrastructure providing software tools to process, analyse and annotate metabolomics data. 

User interface is based on the Galaxy platform (homepage: https://galaxyproject.org/). It is an open, web-based platform for data intensive biomedical research. 
Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses.

 
Configuration
-----------

### Requirement:
 * R software: version > 3.0.0 recommended
 * Specific R libraries: 'argparse' (for the Galaxy wrapper only)


Technical description
-----------

Main files:

- Script_intensity_check.R: R function (core script)
- wrapper_intensity_check.R: R script to link the main R function to inputs
- xml_intensity_check.xml: XML wrapper (interface for Galaxy)


Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10


License
-----------

 * Cea Cnrs Inria Logiciel Libre License, version 2.1 (CECILL-2.1)
