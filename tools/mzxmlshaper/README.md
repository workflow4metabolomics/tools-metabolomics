# mzXML Shaper


Metadata
-----------

 * **@name**: mz(X)ML Shaper
 * **@galaxyID**: mzxmlshaper
 * **@version**: 1.0.0
 * **@authors**: Original code: Quentin Ruin (Contractual engineer - PFEM - INRAE - MetaboHUB) - Maintainer: Melanie Petera (PFEM - INRAE - MetaboHUB)
 * **@init date**: 2024, September
 * **@main usage**: This tool enables the conversion netCDF, mzML or mzXML files into W4M's XCMS mz(X)ML supported file formats

 
Context
-----------

The tool was created to cope with unsupported file formats that may not be read by W4M's Galaxy XCMS.

It can be used for any purpose necessitating standardized mzML or mzXML files, be it visualization in a third-party software, local workflows or W4M XCMS workflows.  

User interface is based on the Galaxy platform (homepage: https://galaxyproject.org/). It is an open, web-based platform for data intensive biomedical research. 
Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses.

 
Configuration
-----------

### Requirement:
 * R software: version = 4.3.3 recommended
 * bioconductor-msdata = 0.42.0
 * bioconductor-mzr = 2.36.0
 * r-w4mrutils = 1.0.0

Technical description
-----------

Main files:

- mzXMLShaper.R: R function (core script)
- mzXMLShaper.xml: XML wrapper (interface for Galaxy)
- 111-1_POS01.CDF (test file): Riker Metabolome Database (http://metabobank.riken.jp/metabo/db/plantMetabolomics/http:/metadb.riken.jp/db/plantMetabolomics/0.1/File/RPMM0054_111-1)
- BlancFin_POS_RA1_1_6869.mzML (test file): Internal blank sample, PFEM - INRAE (https://pfem.isc.inrae.fr/)
- example.mzXML (test file): PRIDE Toolsuite (https://github.com/PRIDE-Toolsuite/inspector-example-files/blob/master/peak-files/example.mzXML.gz)


Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10
                     


License
-----------

 * GPL-3.0-or-later
