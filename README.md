<<<<<<< HEAD
# Galaxy tools for metabolomics

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-xcms/README.html) ![Galaxy Tool Linting and Tests for push and PR](https://github.com/workflow4metabolomics/tools-metabolomics/workflows/Galaxy%20Tool%20Linting%20and%20Tests%20for%20push%20and%20PR/badge.svg)

## Purpose
This repository aims to gather tools and contributors from the metabolomics world.

It is maintained by the Workflow4Metabolomics project but open to any contributors.

Tools themselves should stick with the [IUC](https://github.com/galaxyproject/tools-iuc) (Galaxy Intergalactic Utilities Commission) [standards and best practices](https://galaxy-iuc-standards.readthedocs.io/en/latest/)

### Workflow4Metabolomics
The [Workflow4Metabolomics](http://workflow4metabolomics.org), W4M in short, is a French infrastructure offering software tool processing, analyzing and annotating metabolomics data.

## Related open source projects

### Galaxy
[Galaxy](https://galaxyproject.org/) is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses.

### Dependencies using Conda
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-xcms/README.html)

[Conda](http://conda.pydata.org/) is a package, dependency and environment management for any language—Python, R, Ruby, Lua, Scala, Java, JavaScript, C/ C++, FORTRAN, and more.


### Planemo
[Planemo](https://planemo.readthedocs.io/en/latest/readme.html) is a command-line utilities to assist in developing Galaxy

## Other informations
### Job Dynamic Destination Mapping
Some tools implement a [Job Dynamic Destination Mapping](https://docs.galaxyproject.org/en/latest/admin/jobs.html#dynamic-destination-mapping), like [xcmsSet](https://github.com/workflow4metabolomics/tools-metabolomics/tree/master/tools/xcms_xcmsset#job-dynamic-destination-mapping)


## Tool listing

### [camera_annotate](tools/camera_annotate)
The R-package CAMERA is a Collection of Algorithms for MEtabolite
pRofile Annotation. Its primary purpose is the annotation and evaluation of
LC-MS data. It includes algorithms for annotation of isotope peaks, adducts
and fragments in peak lists. Additional methods cluster mass signals that
originate from a single metabolite, based on rules for mass differences and
peak shape comparison. To use the strength of already existing programs,
CAMERA is designed to interact directly with processed peak data from the
R-package **xcms**.
**What it does?**
The CAMERA annotation procedure can be split into two parts: We want to answer the questions which peaks occur from the same molecule and secondly compute its exact mass and annotate the ion species. Therefore CAMERA annotation workflow contains following primary functions: 1. peak grouping after retention time (**groupFWHM**) 2. peak group verification with peakshape correlation (**groupCorr**) Both methods separate peaks into different groups, which we define as ”pseu- dospectra”. Those pseudospectra can consists from one up to 100 ions, de- pending on the molecules amount and ionizability. Afterwards the exposure of the ion species can be performed with: 2 1. annotation of possible isotopes (**findIsotopes**) 2. annotation of adducts and calculating hypothetical masses for the group (**findAdducts**) This workflow results in a data-frame similar to a xcms peak table, that can be easily stored in a comma separated table .csv (Excel-readable).
If you have two or more conditions, it will return a diffreport result within the annotation results.
The diffreport result shows the most significant differences between two sets of samples. Optionally
create extracted ion chromatograms for the most significant differences.


### [camera_combinexsannos](tools/camera_combinexsannos)
This function check annotations of ion species with the help of a sample from opposite ion mode.
As first step it searches for pseudospectra from the positive and the negative sample within a reten-
tion time window. For every result the m/z differences between both samples are matched against
specific rules, which are combinations from pos. and neg. ion species. As example M+H and M-H
with a m/z difference of 2.014552. If two ions matches such a difference, the ion annotations are
changed (previous annotation is wrong), confirmed or added. Returns the peaklist from one ion
mode with recalculated annotations.


### [camera_macro](tools/camera_macro)
xml macros for other camera repos

### [camera_repository_suite](tools/camera_repository_suite)
xml file describing dependencies for other camera repos

### [correlation_analysis](tools/correlation_analysis)

This tool takes as inputs either tabular table files from the metabolomic workflow (variableMetadata, dataMatrix and sampleMetadata) or a table file of your own
and can execute three different functions ("sorting", "corrdel" and "corr_matrix").

**The "sorting" function:** *used for metabolomic workflow*

 1) First of all, it sorts the data by pcgroup.
 2) It computes the mean operation of all the signal values of the metabolites by sample, and put the results in a new column "signal_moy".
 3) It finally creates a tabular output "sorted_variableMetadata.tsv".

**The "corrdel" function:** *used for metabolomic workflow*

 **For each pcgroup** of the previous sorted tabular file "sorted_table.tsv", it does the following things:
 - it computes a correlation matrix
 - it determines the metabolites which are not correlated to others from the same pcgroup based on the threshold value filled in the "Correlation threshold for pcgroup" parameter
 - the metabolites are sorted by the mean signal intensity (form the highest to the lowest), and each metabolite is tested to the previous ones in the list ; if the tested metabolite is at least correlated to one previous one, it is tagged as DEL (for "deleted", written in a column called "suppress")

It creates two additional tabular files:
 - "correlation_matrix_selected.tsv" (correlation matrix of selected metabolites only)
 - "sif_table.tsv" (for visualization in CytoScape, based on selected metabolites and "Cytoscape correlation threshold" filled value)


**The "corr_matrix" function:** *used for user table file*

 | It computes a correlation matrix named "correlation_matrix.tsv" and creates a sif file named "sif_table.tsv" (for visualization in CytoScape).


### [genform](tools/genform)
Genform generates candidate molecular formulas from high-resolution
MS data. It calculates match values (MV) that show how well candidate molecular formulas fit the MS
isotope peak distributions (MS MV) and the high-resolution MS/MS fragment peak masses (MSMS
MV). Finally it computes a combined match value from these two scores. This software can be
regarded as a further development of the ElCoCo and MolForm modules of MOLGEN-MS with a clear
specialization towards MS/MS.


### [influx_si](tools/influx_si)

Optimize free fluxes and optionaly metabolite concentrations of a given static metabolic network defined in an FTBL file to fit 13C data provided in the same FTBL file.

### [ipo](tools/ipo)
IPO.ipo4xcmsSet
A Tool for automated Optimization of XCMS Parameters

IPO.ipo4xcmsSet
A Tool for automated Optimization of XCMS Parameters


### [isoplot](tools/isoplot)
 We strongly encourage you to read the `documentation <https://isoplot.readthedocs.io/en/latest/>`_ before using Isoplot.


### [msnbase_readmsdata](tools/msnbase_readmsdata)
Reads as set of XML-based mass-spectrometry data files and
generates an MSnExp object. This function uses the functionality
provided by the ‘mzR’ package to access data and meta data in
‘mzData’, ‘mzXML’ and ‘mzML’.


### [nmr_annotation](tools/nmr_annotation)
ASICS, based on a strong statistical theory, handles automatically the metabolite identification and quantification


### [nmr_annotation2d](tools/nmr_annotation2d)
BARSA is an automatic algorithm for bi-dimensional NMR spectra annotation


### [nmr_preprocessing](tools/nmr_preprocessing)
Spectra preprocessing

These steps correspond to the following steps in the PEPS-NMR R library (https://github.com/ManonMartin/PEPSNMR):

* Group Delay suppression (First order phase correction)
* Removal of solvent residuals signal from the FID
* Apodization to increase the Signal-to-Noise ratio of the FID
* Fourier transformation
* Zero order phase correction
* Shift referencing to calibrate the spectra with internal compound referencing
* Baseline correction
* Setting of negatives values to 0

NMR Read

Nuclear Magnetic Resonance Bruker files reading (from the PEPS-NMR R package (https://github.com/ManonMartin/PEPSNMR))


### [normalization](tools/normalization)

Normalization (operation applied on each (preprocessed) individual spectrum) of preprocessed data


### [xcms_export_samplemetadata](tools/xcms_export_samplemetadata)
xcms get sampleMetadata
This tool generates a skeleton of sampleMetadata with perhaps some strange sample names which are definitely compatible with xcms and R
This sampleMetadata file have to be filled with extra information as the class, batch information and maybe conditions


### [xcms_fillpeaks](tools/xcms_fillpeaks)
xcms fillChromPeaks
**Integrate areas of missing peaks**
For each sample, identify peak groups where that sample is not
represented. For each of those peak groups, integrate the signal
in the region of that peak group and create a new peak.


### [xcms_group](tools/xcms_group)
xcms groupChromPeaks

After peak identification with xcmsSet, this tool groups the peaks which represent the same analyte across samples using overlapping m/z bins and calculation of smoothed peak distributions in chromatographic time. Allows rejection of features, which are only partially detected within the replicates of a sample class.

### [xcms_macro](tools/xcms_macro)
xml macros for other xcms repos

### [xcms_merge](tools/xcms_merge)
xcms findChromPeaks Merger
This tool allows you to run one xcms findChromPeaks process per sample in parallel and then to merge all RData images into one.
The result is then suitable for xcms groupChromPeaks.
You can provide a sampleMetadata table to attribute phenotypic values to your samples.


### [xcms_plot_chromatogram](tools/xcms_plot_chromatogram)
xcms plot chromatogram
This tool will plot Base Peak Intensity chromatogram (BPI) and Total Ion Current chromatogram (TIC) from xcms experiments.


### [xcms_refine](tools/xcms_refine)
xcms refineChromPeaks

After peak identification with xcms findChromPeaks (xcmsSet), this tool refines those peaks.
It either removes peaks that are too wide or removes peaks with too low intensity or combines peaks that are too close together.
Note well that refineChromPeaks methods will always remove feature definitions,
because a call to this method can change or remove identified chromatographic peaks, which may be part of features.
Therefore it must only be run immediately after findChromPeaks (xcmsSet).

### [xcms_repository_suite](tools/xcms_repository_suite)
xml file describing dependencies for other xcms repos

### [xcms_retcor](tools/xcms_retcor)
xcms adjustRtime
After matching peaks into groups, xcms can use those groups to identify and correct
correlated drifts in retention time from run to run. The aligned peaks can then be
used for a second pass of peak grouping which will be more accurate than the first.
The whole process can be repeated in an iterative fashion. Not all peak groups will be helpful
for identifying retention time drifts. Some groups may be missing peaks from a large
fraction of samples and thus provide an incomplete picture of the drift at that time point.
Still others may contain multiple peaks from the same sample, which is a sign of impropper grouping.


### [xcms_summary](tools/xcms_summary)
xcms process history
This tool provide a HTML summary which summarizes your analysis using the [W4M] XCMS and CAMERA tools


### [xcms_test-data](tools/xcms_test-data)
test data repo for xcms tool suit

### [xcms_xcmsset](tools/xcms_xcmsset)
xcms findChromPeaks
This tool is used for preprocessing data from multiple LC/MS files (NetCDF, mzXML and mzData formats) using the xcms_ R package.
It extracts ions from each sample independently, and using a statistical model, peaks are filtered and integrated.
A tutorial on how to perform xcms preprocessing is available as GTN_ (Galaxy Training Network).


## Historic contributors (non cited by GitHub)
- Urszula Czerwinska [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [CNRS](www.cnrs.fr)/[Sorbonne Université](http://www.sorbonne-universite.fr/) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
- Marion Landi [PFEM](https://www6.ara.inra.fr/plateforme_exploration_metabolisme) / [MetaboHUB](https://www.metabohub.fr/home.html) - [INRA](http://www.inra.fr/) - France
- Misharl Monsoor [@mmonsoor](https://github.com/mmonsoor) - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [CNRS](www.cnrs.fr)/[Sorbonne Université](http://www.sorbonne-universite.fr/) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
- Pierre Pericard [@ppericard](https://github.com/ppericard)- [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [CNRS](www.cnrs.fr)/[Sorbonne Université](http://www.sorbonne-universite.fr/) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France

=======
NMR Bucketing for Galaxy
========================

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io) [![Build Status](https://travis-ci.org/workflow4metabolomics/nmr_bucketing.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/nmr_bucketing)

Our project
-----------
The [Workflow4Metabolomics](http://workflow4metabolomics.org), W4M in short, is a French infrastructure offering software tool processing, analyzing and annotating metabolomics data. It is based on the Galaxy platform.


NMR Bucketing
-------------

Bucketing / Binning (spectra segmentation in fixed-size windows) and integration (sum of absolute intensities inside each bucket) to preprocess NMR data


Galaxy
------
Galaxy is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 

Homepage: [https://galaxyproject.org/](https://galaxyproject.org/)


Dependencies using Conda
------------------------
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io)

The main recipe: [https://github.com/bioconda/bioconda-recipes/tree/master/recipes/r-pracma](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/r-pracma)

```
#To install miniconda2
#http://conda.pydata.org/miniconda.html
#To install the needed R library using conda:
conda install r-batch r-pracma
#To set an environment:
conda create -n nmr_bucketing r-batch r-pracma`
#To activate the environment:
. activate nmr_bucketing
```

[Conda](http://conda.pydata.org/) is package manager that among many other things can be used to manage Python packages.

Travis
------
[![Build Status](https://travis-ci.org/workflow4metabolomics/nmr_bucketing.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/nmr_bucketing)

Test and Deploy with Confidence. Easily sync your GitHub projects with Travis CI and you'll be testing your code in minutes!

Historic contributors
---------------------
 - Marie Tremblay-Franco @mtremblayfr - [French Metabolomics and Fluxomics Infrastructure (MetaboHUB)](http://www.metabohub.fr/en) - [MetaToul](http://www.metatoul.fr/)
 - Marion Landi - [French Metabolomics and Fluxomics Infrastructure (MetaboHUB)](http://www.metabohub.fr/en) - [La plateforme "Exploration du Métabolisme" (PFEM, Clermont-Ferrand)](http://www6.clermont.inra.fr/plateforme_exploration_metabolisme)
 - Franck Giacomoni @fgiacomoni - [French Metabolomics and Fluxomics Infrastructure (MetaboHUB)](http://www.metabohub.fr/en) - [La plateforme "Exploration du Métabolisme" (PFEM, Clermont-Ferrand)](http://www6.clermont.inra.fr/plateforme_exploration_metabolisme)
 - Gildas Le Corguillé @lecorguille - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [UPMC](www.upmc.fr)/[CNRS](www.cnrs.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
>>>>>>> 1933f067e (add info about conda)
