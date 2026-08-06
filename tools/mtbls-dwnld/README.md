Metabolights Downloader
=======================

[![Build Status](https://travis-ci.org/workflow4metabolomics/mtbls-dwnld.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/mtbls-dwnld)

A [Metabolights](http://www.ebi.ac.uk/metabolights/) study downloader for [Galaxy](https://galaxyproject.org/), part of the [Workflow4Metabolomics](http://workflow4metabolomics.org/) project, and developed during the [PhenoMeNal](http://phenomenal-h2020.eu/home/) project.

For more information, see the galaxy tool page, help section, available inside `mtbls-dwnld_config.xml`.

## Requirements

 * Python 3.8
 * isatools 0.10.3
 * `unzip` program.
 * `wget` program.
 * `ascp` (Aspera secure copy) program, version 3.7.4 or greater, which you can download [here](http://downloads.asperasoft.com/en/downloads/62).
 
## Updates

### 4.1.4

 * Use MTBLS29 instead of MTBLS30, which has been removed from public repository, for testing.
 * Make planemo tests pass on Travis-CI.
 * Run several times ascp downloads for tests on Travis-CI, due to connection issues.

### 2.0.4

 * Test presence of downloader program.
 * Correct running of `wget`.

### 2.0.3

 * Allow to choose the downloading method (`wget` or `aspera`).
 * Ask to specify is the study to download is public or private.
 * Allow to download all assays at once and convert them to W4M format.
 * Allow to download *mzMl* and *mzData* files as collections. 

### 1.3.0

 * Allow to choose the assay to extract and convert to W4M format.
