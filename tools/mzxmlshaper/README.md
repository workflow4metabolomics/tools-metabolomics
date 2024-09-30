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
- ko16.CDF: netCDF test file from faahKO Bioconductor package (https://www.bioconductor.org/packages/release/data/experiment/html/faahKO.html)
- ko16.mzml: ko.CDF file converted to mzML using mzR R package
- ko16.mzXml: ko.CDF file converted to mzXML using mzR R package

Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10
                     


License
-----------

 * GPL-3.0-or-later
