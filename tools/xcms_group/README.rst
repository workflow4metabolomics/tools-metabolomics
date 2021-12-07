
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

**Version 3.0.0.1 - 09/11/2018**

- BUGFIX: issue when the vector at peakidx is too long and is written in a new line during the export of the peaklist

**Version 3.0.0.0 - 08/03/2018**

- UPGRADE: upgrade the xcms version from 1.46.0 to 3.0.0. So refactoring of a lot of underlying codes and methods. Some parameters may have been renamed.

- NEW: a bunch of new options: PeakDensity.minSamples), MzClust.minSamples)

- NEW: a new density plot

- IMPROVEMENT: the advanced options are now in sections. It will allow you to access to all the parameters and to know their default values.


**Version 2.1.1 - 29/11/2017**

- BUGFIX: To avoid issues with accented letter in the parentFile tag of the mzXML files, we changed a hidden mechanim to LC_ALL=C


**Version 2.1.0 - 07/02/2017**

- IMPROVEMENT: Add an option to export the peak list at this step without have to wait camara.annotate

- IMPROVEMENT: xcms.group can deal with merged individual data from "xcms.xcmsSet Merger"

- BUGFIX: the default value of "density" -> "Maximum number of groups to identify in a single m/z slice" which was of 5 have been changed to fix with the XMCS default values to 50


**Version 2.0.6 - 06/07/2016**

- UPGRADE: upgrate the xcms version from 1.44.0 to 1.46.0


**Version 2.0.5 - 04/04/2016**

- TEST: refactoring to pass planemo test using conda dependencies


**Version 2.0.4 - 10/02/2016**

- BUGFIX: better management of errors. Datasets remained green although the process failed

- UPDATE: refactoring of internal management of inputs/outputs

- UPDATE: refactoring to feed the new report tool


**Version 2.0.2 - 02/06/2015**

- IMPROVEMENT: new datatype/dataset formats (rdata.xcms.raw, rdata.xcms.group, rdata.xcms.retcor ...) will facilitate the sequence of tools and so avoid incompatibility errors.

- IMPROVEMENT: parameter labels have changed to facilitate their reading.
