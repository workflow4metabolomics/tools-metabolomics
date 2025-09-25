
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.0.0] - 2025-07-10
### Added
- Specific names for the 'sampleType', 'injectionOrder', and 'batch' from sampleMetadata are now available in a dedicated parameter section.
- Addition of a sum of ions before/after plot for linear/lowess/loess methods.
- Addition of a third option in "Null values" parameter (renamed "unconsistant values") in linear/lowess/loess methods.
- linear/lowess/loess methods now handle NA in intensities and allow "blank" samples in the dataset.

### Changed
- XML optimisation using macros.
- Output name changes.
- linear/lowess/loess methods: disabling of RData output.
- linear/lowess/loess methods: split of tool-linked code and script-linked one.
- linear/lowess/loess methods: adjustments in the normalisation process to match matters linked to NA acceptance.
- linear/lowess/loess methods: better handling of special characters in IDs and column names.

## [2.2.4] - 2024-xx-xx
### Fixed
- Fixed bug for pool selection ("all_loess" methods).

## [2.2.2] - 2024-xx-xx
### Fixed
- Fixed bug for color plot ("all_loess" methods).

## [2.2.0] - 2024-xx-xx
### Added
- Specific names for the 'sampleType', 'injectionOrder', and 'batch' from sampleMetadata can be selected by the user (for compatibility with the MTBLS downloader).

## [2.1.2] - 2024-xx-xx
### Changed
- Minor modifications in config file.

## [2.1.0] - 2024-xx-xx
### Changed
- For PCA figure display only (**all_loess** options): missing values are set to the minimum value before PCA computation is performed (with svd).
- Additional running and installation tests added with planemo, conda, and travis.

### Fixed
- Variables with NA or 0 values in all reference samples are discarded before applying the **all_loess** normalization.

### Changed
- Modifications of the **all_loess_wrapper** file to handle the recent **ropls** package versions (i.e. 1.3.15 and above) which use S4 classes.