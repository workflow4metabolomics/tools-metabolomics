
Changelog/News
--------------

**Version 2.2.6+camera1.48.0 - 11/06/2021**

- UPGRADE: upgrade the CAMERA version from 1.46.0 to 1.48.0 (see CAMERA News_)

**Version 2.2.6+camera1.46.0 - 12/04/2021**

- UPGRADE: upgrade the CAMERA version from 1.42.0 to 1.46.0 (see CAMERA News_)

**Version 2.2.6+camera1.42.0+galaxy1 - 09/03/2020**

- BUGFIX: Fix the zip export of the pictures (eic and boxplot)

**Version 2.2.5+camera1.42.0 - 13/02/2020**

- UPGRADE: upgrade the CAMERA version from 1.38.1 to 1.42.0 (see CAMERA News_)

**Version 2.2.5 - 09/04/2019**

- NEW: zip export are back for pictures (eic and boxplot) and diffreport tables

- UPGRADE: upgrade the CAMERA version from 1.34.0 to 1.38.1 (see CAMERA News_)

- UPGRADE: refactoring of internal code

**Version 2.2.4 - 09/10/2018**

- NES: CAMERA.annotate no longer export a DataMatrix. fillChromPeaks does the job

**Version 2.2.3 - 30/04/2018**

- NEW: support the new xcms 3.0.0 wrapper

**Version 2.2.2 - 01/03/2018**

- UPGRADE: upgrate the CAMERA version from 1.26.0 to 1.32.0

**Version 2.2.1 - 29/11/2017**

- BUGFIX: To avoid issues with accented letter in the parentFile tag of the mzXML files, we changed a hidden mechanim to LC_ALL=C

**Version 2.2.0 - 03/02/2017**

- BUGFIX: the diffreport ids didn't convert the rt in minute as the other export

- UPDATE: the settings (digits, convertion in minutes) of the identifiants will no longer modify the native one. Because we want to be conservative and because it can be dangerous for the data integrity during a futur merge of the table, we decide to put those customization in a new column namecustom within the variableMetadata.

- IMPROVEMENT: add some sections within the form to separate the different parts of the process

- IMPROVEMENT: add the possibility to use defined ruleset

- IMPROVEMENT: add the possibility to set the MZ digit within the identifiants

- IMPROVEMENT: CAMERA.annotate is now compatible with merged individual data from xcms.xcmsSet


**Version 2.1.5 - 21/04/2016**

- UPGRADE: upgrate the CAMERA version from 1.22.0 to 1.26.0


**Version 2.1.4 - 18/04/2016**

- TEST: refactoring to pass planemo test using conda dependencies

**Version 2.1.3 - 10/02/2016**

- BUGFIX: better management of errors. Datasets remained green although the process failed

- BUGFIX: the conversion into minutes of the retention time was applied to the diffreport outputs (several conditions)

- IMPROVEMENT: when there are several conditions, the tool will generate individual datasets (tsv, pdf) instead of a zip file. The usual png (eic, boxplot) will from now be integrated in two pdf.

- UPDATE: refactoring of internal management of inputs/outputs


**VERSION 2.1.0 - 09/10/2015**

- BUGFIX: There was a bug with the CAMERA.annotate (generating a bad dataMatrix (intensities which don't match with the metabolites))


**VERSION 2.1.0 - 07/06/2015**

- IMPROVEMENT: new datatype/dataset formats (rdata.camera.positive, rdata.camera.negative, rdata.camera.quick ...) will facilitate the sequence of tools and so avoid incompatibility errors.

- IMPROVEMENT: parameter labels have changed to facilitate their reading.

- UPDATE: merged with annotateDiffreport. Some parameters are dedicated to experiences with several conditions
