# mz(X)MLMetadataCollector


Metadata
-----------

 * **@name**: mz(X)MLMetadataCollector
 * **@galaxyID**: mzxmlmetadatacollector
 * **@version**: 1.0.0
 * **@authors**: Original code: Quentin Ruin (Contractual engineer - PFEM - INRAE - MetaboHUB) - Maintainer: Melanie Petera (PFEM - INRAE - MetaboHUB)
 * **@init date**: 2024, November
 * **@main usage**: This tool enables the gathering of metadata encapsulated in the headers of mzML and mzXML files 

 
Context
-----------

The tool was created to efficiently and quickly collect essential metadata of mzML and mzXML files.

It is based on the parsing of the XML of the files, whose entries vary with the constructors and the conversion method used. The metadata header is streamed and parsing stops at the end of the first spectrum, so the encoded payloads of the remaining scans are never held in memory. Summing the spectral and chromatographic points needs a second pass over the whole file, in which every element is discarded as soon as it is complete, so the memory used stays independent of the number of scans. A missing metadata entry is reported as "Not found", while an unreadable file (missing input, malformed XML, unsupported extension, or root element contradicting the file extension) makes the tool fail.

User interface is based on the Galaxy platform (homepage: https://galaxyproject.org/). It is an open, web-based platform for data intensive biomedical research. 
Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses.

 
Configuration
-----------

### Requirement:
 * Python: version = 3.10
 * numpy: version =  2.1.0

Technical description
-----------

Main files:

- mzXMLMetadataCollector.py: Python function (core script)
- mzXMLMetadataCollector.xml: XML wrapper (interface for Galaxy)
- metadata_multiple.tabular: test file, expected resultats for the "multiple files or collection" parameter test
- metadata_single.tabular: test file, expected resultats for the "single file" parameter test
- metadata_mzxml.tabular: test file, expected resultats for the "single file" parameter test on an mzXML input

Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10
                     


License
-----------

 * GPL-3.0-or-later
