Job Dynamic Destination Mapping
-------------------------------

**Why**

xcmsSet wrapper allow both individual file and zip file which can contain several samples.
Thus, it can be interesting to adjust the number of thread according to the input type.
For example: 1 thread for a single mzXML or NetCDF file and 8 threads for a zip file.

**What**

The [Dynamic Destination Mapping](https://galaxyproject.org/admin/config/jobs/#dynamic-destination-mapping) allow Galaxy to choose a destination at runtime based on factors such as the job inputs, user submitting the job, cluster status, etc...

**How**

1. The file [`destinations_input_type.py`](https://raw.githubusercontent.com/workflow4metabolomics/xcms/master/galaxy/xcms_xcmsset/destinations_input_type.py) (shiped with this tool) must be place in `lib/galaxy/jobs/rules/`
2. The `job_conf.xml` must be inspired by the [`job_conf.xml.sample`](https://github.com/workflow4metabolomics/xcms/blob/master/galaxy/xcms_xcmsset/job_conf.xml.sample) shiped with this tool
3. The final destination names must match between the `job_conf.xml` and the `destinations_input_type.py`
4. Restart Galaxy


Changelog/News
--------------

.. _News: https://bioconductor.org/packages/release/bioc/news/xcms/NEWS

**Version 3.12.0+galaxy* - 03/03/2020**

- UPGRADE: upgrade the xcms version from 3.6.1 to 3.12.0 (see XCMS News_)

**Version 3.6.1.0 - 03/09/2019**

- UPGRADE: upgrade the xcms version from 3.4.4 to 3.6.1 (see XCMS News_)

**Version 3.4.4.1 - 30/04/2019**

- BUGFIX: remove the pre-compute of the chromatograms which was memory consuming. Now, only xcms plot chromatogram will generate the Chromatograms.

**Version 3.4.4.0 - 08/02/2019**

- UPGRADE: upgrade the xcms version from 3.0.0 to 3.4.4 (see XCMS News_)

**Version 3.0.0.0 - 08/03/2018**

- UPGRADE: upgrade the xcms version from 1.46.0 to 3.0.0. So refactoring of a lot of underlying codes and methods. Some parameters may have been renamed.

- NEW: a bunch of new options: Spectra Filters (previously scanrange), CentWave.(mzCenterFun, fitgauss, verboseColumns), MatchedFilter.(sigma, impute, baseValue, max), MSW.(verboseColumns), ...

- UPDATE: since xcms 3.0.0, some options are no more available: scanrange (replace by filters), profmethod, MatchedFilter.step, MatchedFilter.sigma, MSW.winSize.noise, MSW.SNR.method

- IMPROVEMENT: the advanced options are now in sections. It will allow you to access to all the parameters and to know their default values.

- IMPROVEMENT: the tool "should" be now more flexible in term of file naming: it "should" accept space and comma. But don't be too imaginative :)

- CHANGE: removing of the TIC and BPC plots. You can new use the dedicated tool "xcms plot chromatogram"


**Version 2.1.1 - 29/11/2017**

- BUGFIX: To avoid issues with accented letter in the parentFile tag of the mzXML files, we changed a hidden mechanim to LC_ALL=C


**Version 2.1.0 - 22/02/2017**

- NEW: The W4M tools will be able now to take as input a single file. It will allow to submit in parallel several files and merge them afterward using "xcms.xcmsSet Merger" before "xcms.group".

- BUGFIX: the default value of "matchedFilter" -> "Step size to use for profile generation" which was of 0.01 have been changed to fix with the XMCS default values to 0.1


**Version 2.0.11 - 22/12/2016**

- BUGFIX: propose scanrange for all methods


**Version 2.0.10 - 22/12/2016**

- BUGFIX: when having only one group (i.e. one folder of raw data) the BPC and TIC pdf files do not contain any graph


**Version 2.0.9 - 06/07/2016**

- UPGRADE: upgrade the xcms version from 1.44.0 to 1.46.0


**Version 2.0.8 - 06/04/2016**

- TEST: refactoring to pass planemo test using conda dependencies


**Version 2.0.7 - 10/02/2016**

- BUGFIX: better management of errors. Datasets remained green although the process failed

- BUGFIX/IMPROVEMENT: New checking steps around the imported data in order to raise explicte error message before or after launch XCMS: checking of bad characters in the filenames, checking of the XML integrity and checking of duplicates which can appear in the sample names during the XCMS process because of bad characters

- BUGFIX/IMPROVEMENT: New step to check and delete bad characters in the XML: accented characters in the storage path of the mass spectrometer

- UPDATE: refactoring of internal management of inputs/outputs

- TEST: refactoring to feed the new report tool


**Version 2.0.2 - 18/01/2016**

- BUGFIX: Some zip files were tag as "corrupt" by R. We have changed the extraction mode to deal with thoses cases.


**Version 2.0.2 - 09/10/2015**

- BUGFIX: Some users reported a bug in xcms.xcmsSet. The preprocessing stops itself and doesn't import the whole dataset contained in the zip file without warning. But meanwhile, please check your samplemetadata dataset and the number of rows.


**Version 2.0.2 - 02/06/2015**

- NEW: The W4M workflows will now take as input a zip file to ease the transfer and to improve dataset exchange between tools and users. (See How_to_upload). The previous "Library directory name" is still available but we invite user to switch on the new zip system as soon as possible.

- IMPROVEMENT: new datatype/dataset formats (rdata.xcms.raw, rdata.xcms.group, rdata.xcms.retcor ...) will facilitate the sequence of tools and so avoid incompatibility errors.

- IMPROVEMENT: parameter labels have changed to facilitate their reading.
