====================================
GC-MS data processing version 2.1.1
====================================

1) Date : 2019-05-20

2) Author and maintainer :

   Yann GUITTON

   LABERCA - PFCA Metabolomics Platform

   E-mail: yann.guitton@oniris-nantes.fr / yann.guitton@gmail.com

3) Funding :

   Developed within IDEALG project (http://www.idealg.ueb.eu/versionAnglaise/)

   Developed within IFB-MetaboHub W4M project

4) Usage restrictions :

   Use of this tool is restricted to the service conditions of the MetaboHUB-IFB infrastructures.
   For any question regarding the use of these services, please contact: yann.guitton@oniris-nantes.fr

5) Installation :

   4 files are required for installation :

   - 'README.rst'
         Instructions for installation
   
   - 'metaMS_plot.xml'
         Configuration file; to be put into the './galaxy-dist/tools/' directory 
         and 2.png files for illustration

   - 'metaMS_plot.r'
         Wrapper code written in R aimed at launching the runGC function from the metaMS package given the arguments entered by the user through the Galaxy interface
   
   - 'metaMS R package '
         The 'metaMS' package requires dependencies and can be installed with source("http://bioconductor.org/biocLite.R") or `biocLite("metaMS")`
 
         This code is for installation of the Galaxy module on the Workflow4metabolomics.org MetaboHUB-IFB platform only and must not be distributed without the author agreement

   
Changelog/News
--------------
**Version 1.0 - 20/05/2019**

- NEW : new tool extract from previous metaMS_runGC tool. EICs have been corrected.

Test Status
-----------

Planemo test using conda: passed


====================================
GC-MS data processing version 3.0.0
====================================

1) Date : 2019-05-20

2) Author and maintainer :

   Yann GUITTON

   LABERCA - PFCA Metabolomics Platform

   E-mail: yann.guitton@oniris-nantes.fr / yann.guitton@gmail.com

3) Funding :

   Developed within IDEALG project (http://www.idealg.ueb.eu/versionAnglaise/)

   Developed within IFB-MetaboHub W4M project

4) Usage restrictions :

   Use of this tool is restricted to the service conditions of the MetaboHUB-IFB infrastructures.
   For any question regarding the use of these services, please contact: yann.guitton@oniris-nantes.fr

5) Installation :

   4 files are required for installation :

   - 'README.rst'
         Instructions for installation
   
   - 'metaMS_runGC.xml'
         Configuration file; to be put into the './galaxy-dist/tools/' directory 
         and 2.png files for illustration

   - 'metaMS_runGC.r'
         Wrapper code written in R aimed at launching the runGC function from the metaMS package given the arguments entered by the user through the Galaxy interface
   
   - 'metaMS R package '
         The 'metaMS' package requires dependencies and can be installed with source("http://bioconductor.org/biocLite.R") or `biocLite("metaMS")`
 
         This code is for installation of the Galaxy module on the Workflow4metabolomics.org MetaboHUB-IFB platform only and must not be distributed without the author agreement

   
Changelog/News
--------------
**Version 3.0 - 20/05/2019**

- Divided tool into two new. One will be metaMS_plot and this one stay the same without plotting chromatograms and without zip files

**Version 2.0 - 15/05/2017**

- RI option, better usage of ri filter and easier file reading for RI database

**Version 1.1 - 11/07/2016**

- TEST: refactoring to pass planemo test using conda dependencies

**Version 1.0 - 01/06/2015**

- NEW: new tool

Test Status
-----------

Planemo test using conda: passed
