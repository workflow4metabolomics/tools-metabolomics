# mzXML Shaper


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

It is based on the recognition of XML-like tags in the headers of the files, that vary with the constructors and the conversion method used.

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
- BlancFin_POS_RA1_1_6869.mzML: test file (truncated to be lighter)
- BlancFin_POS_RA1_1_6869.mzXML: test file (truncated to be lighter)
- Blank1_POS_RA1_1_6727.mzXML: test file (truncated to be lighter)
- Bruker_inhouse_mzml.mzML: test file (truncated to be lighter)
- Impact_DA53_mzml.mzML: test file (truncated to be lighter)
- OE240_DS_QC_ID_01.mzML: test file (truncated to be lighter)
- OE240_DS_QC_ID_012.mzML: test file (truncated to be lighter)
- Pos_Xenobio_DDA_30_35CE_incl_90_Spk_3_A.mzML: test file (truncated to be lighter)
- Pos_Xenobio_DDA_30_35CE_incl_90_Spk_3_A2.mzML: test file (truncated to be lighter)
- QC_20autoMSMS_10_1_2969.mzXML: test file (truncated to be lighter)
- QC_autoMSMS_10_1_2969.mzML: test file (truncated to be lighter)

Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10
                     


License
-----------

 * GPL-3.0-or-later
