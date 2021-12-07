
Changelog/News
--------------

.. _News: https://bioconductor.org/packages/release/bioc/news/xcms/NEWS

**Version 3.12.0+galaxy* - 03/03/2020**

- UPGRADE: upgrade the xcms version from 3.6.1 to 3.12.0 (see XCMS News_)

**Version 3.6.1.0 - 03/09/2019**

- UPGRADE: upgrade the xcms version from 3.4.4 to 3.6.1 (see XCMS News_)

**Version 3.4.4.0 - 08/02/2019**

- UPGRADE: upgrade the xcms version from 3.0.0 to 3.4.4 (see XCMS News_)

- BUGFIX: issue with Inf values in the exported DataMatrix: https://github.com/sneumann/xcms/issues/323#issuecomment-433044378

**Version 3.0.0.2 - 09/11/2018**

- BUGFIX: issue when the vector at peakidx is too long and is written in a new line during the export of the peaklist


**Version 3.0.0.1 - 09/10/2018**

- IMPROVEMENT: the export of the PeakList is now mandatory


**Version 3.0.0.0 - 08/03/2018**

- UPGRADE: upgrade the xcms version from 1.46.0 to 3.0.0. So refactoring of a lot of underlying codes and methods. Some parameters may have been renamed.

- UPDATE: since xcms 3.0.0, the selection of a method is no more needed (chrom or MSW). xcms will detect from the data the peak picking method used in findChromPeaks

- UPDATE: since xcms 3.0.0, new parameters are available: expandMz, expandRt and ppm


**Version 2.1.1 - 29/11/2017**

- BUGFIX: To avoid issues with accented letter in the parentFile tag of the mzXML files, we changed a hidden mechanim to LC_ALL=C


**Version 2.1.0 - 07/02/2017**

- IMPROVEMENT: change the management of the peaklist ids. The main ids remain the same as xcms generated. The export setiings now only add custom names in the variableMetadata tab (namecustom)

- IMPROVEMENT: xcms.fillpeaks can deal with merged individual data


**Version 2.0.8 - 22/12/2016**

- IMPROVEMENT: Add an option to export the peak list at this step without having to wait for CAMERA.annotate


**Version 2.0.7 - 06/07/2016**

- UPGRADE: upgrate the xcms version from 1.44.0 to 1.46.0


**Version 2.0.6 - 04/04/2016**

- TEST: refactoring to pass planemo test using conda dependencies


**Version 2.0.5 - 10/02/2016**

- BUGFIX: better management of errors. Datasets remained green although the process failed

- UPDATE: refactoring of internal management of inputs/outputs

- UPDATE: refactoring to feed the new report tool


**Version 2.0.2 - 02/06/2015**

- IMPROVEMENT: new datatype/dataset formats (rdata.xcms.raw, rdata.xcms.group, rdata.xcms.retcor ...) will facilitate the sequence of tools and so avoid incompatibility errors.

- IMPROVEMENT: parameter labels have changed to facilitate their reading.
