
Changelog/News
--------------

.. _News: https://bioconductor.org/packages/release/bioc/news/xcms/NEWS

**Version 3.12.0+galaxy* - 03/03/2020**

- UPGRADE: upgrade the xcms version from 3.6.1 to 3.12.0 (see XCMS News_)

**Version 3.6.1+galaxy1 - 13/02/2020**

- NEW: if a retention time ajustment have already been applied to your data. The function applyAdjustedRtime will replace raw retention times with adjusted retention times and so alloww to cumulate the ajustments.

**Version 3.6.1.0 - 03/09/2019**

- UPGRADE: upgrade the xcms version from 3.4.4 to 3.6.1 (see XCMS News_)

**Version 3.4.4.1 - 30/04/2019**

- BUGFIX: remove the pre-compute of the chromatograms which was memory consuming. Now, only xcms plot chromatogram will generate the Chromatograms.

**Version 3.4.4.0 - 08/02/2019**

- UPGRADE: upgrade the xcms version from 3.0.0 to 3.4.4 (see XCMS News_)

**Version 3.0.0.0 - 08/03/2018**

- UPGRADE: upgrade the xcms version from 1.46.0 to 3.0.0. So refactoring of a lot of underlying codes and methods. Some parameters may have been renamed.

- NEW: a bunch of new options: Obiwarp.(centerSample, response, distFun, gapInit, gapExtend, factorDiag, factorGap, localAlignment, initPenalty)

- IMPROVEMENT: the advanced options are now in sections. It will allow you to access to all the parameters and to know their default values.

- CHANGE: removing of the TIC and BPC plots. You can new use the dedicated tool "xcms plot chromatogram"


**Version 2.1.1 - 29/11/2017**

- BUGFIX: To avoid issues with accented letter in the parentFile tag of the mzXML files, we changed a hidden mechanim to LC_ALL=C


**Version 2.1.0 - 03/02/2017**

- IMPROVEMENT: xcms.retcor can deal with merged individual data


**Version 2.0.8 - 22/12/2016**

- BUGFIX: when having only one group (i.e. one folder of raw data) the BPC and TIC pdf files do not contain any graph


**Version 2.0.7 - 06/07/2016**

- UPGRADE: upgrate the xcms version from 1.44.0 to 1.46.0


**Version 2.0.6 - 04/04/2016**

- TEST: refactoring to pass planemo test using conda dependencies


**Version 2.0.5 - 10/02/2016**

- BUGFIX: better management of errors. Datasets remained green although the process failed

- BUGFIX: some pdf remained empty even when the process succeed

- UPDATE: refactoring of internal management of inputs/outputs

- UPDATE: refactoring to feed the new report tool


**Version 2.0.2 - 02/06/2015**

- IMPROVEMENT: new datatype/dataset formats (rdata.xcms.raw, rdata.xcms.group, rdata.xcms.retcor ...) will facilitate the sequence of tools and so avoid incompatibility errors.

- IMPROVEMENT: parameter labels have changed to facilitate their reading.
