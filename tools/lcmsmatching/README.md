LC/MS matching
==============

[![Build Status](https://travis-ci.org/workflow4metabolomics/lcmsmatching.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/lcmsmatching)

An LC/MS matching tool for [Galaxy](https://galaxyproject.org/), part of the [Workflow4Metabolomics](http://workflow4metabolomics.org/) project, and developed during the [MetaboHUB](http://www.metabohub.fr/en) project.

The two matching algorithms used in this tool have been imported from developments made at [CEA](http://www.cea.fr/english) Saclay, inside the *DSV/IBITEC-S/SPI*. They have been translated from C# to R.

For more information, see the galaxy tool page, help section, available inside `galaxy/lcmsmatching.xml`.

## lcmsmatching script

This is the script, included in this repository, that allows to run on command line an MZ matching on one of the available database types.

Please run `lcmsmatching -h` for a help page listing all options and presenting some examples.

## Dependencies

 * `R` version `3.5.1`.
 * `R` packages:
   - `getopt` >= `1.20.0`.
   - `biodb` >= `1.2.0rc2`.

## Changelog

### 4.0.2

   * Increase getopt version to 1.20.2.

### 4.0.1

   * Downgrade to Galaxy 18.05. Test in both 18.05 and 18.09.

### 4.0.0

   * Switch to biodb R library (<http://github.com/pkrog/biodb>).
   * Remove Excel and 4TabSql databases from script.
   * Remove all dynamic fields in XML (i.e.: fields computed using python scripts, like the list of chromatogaphic columns).
   * Use now a single field for in-house file databases column names, whose value is a comma separated list of key/value pairs.
   * Update Peakforest URL.

### 3.4.3

   * Returns empty match in case of NA values in mz.low and mz.high.
   * Speed up HTML output writing.

### 3.3.1

   * Correct a bug while trying to connect to Peakforest for getting the list of chromatographic columns.

### 3.3.0

   * The file database (in-house) field names are now presented in individual choice lists instead of a single text box where you had to insert a very long keys/values string.
   * The tool now tries to guess the names of the file database fields, the values of the MS mode column, and the names of the input file columns.
   * Allows to select the unit (minutes or seconds) of retention time values inside the input file, but also inside the file database (in-house).
