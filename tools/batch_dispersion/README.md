# Batch dispersion

Metadata
-----------

 * **@name**: Batch dispersion
 * **@galaxyID**: batch_dispersion
 * **@version**: 0.1.6+galaxy1
 * **@authors**: Original code: Brice Mulot (PFEM - UNH - INRAE) - Maintainer: Etienne Jules (PFEM - UNH - INRAE - MetaboHUB)
 * **@init date**: 2025, July
 * **@main usage**: This tool analyses dispersion of ions intensities across analytical batches of metabolomic analysis on quality controls samples.

 
Context
-----------

This tool generates convex hull plots for metabolite intensity data by injection order across different batches. This can be used to assess intensity values dispersion of ions on similar samples across batches and/or projects (QC, pool, reference materials).
 
Configuration
-----------

### Requirement:
 * R software: version = 4.5.2
 * r-ggplot2 = 3.5.2
 * r-optparse = 1.7.5
 * r-dispersionindicators = 0.1.6

Technical description
-----------

Main files:

- convex_hull_analysis.R: R function (core script)
- batch_dispersion.xml: XML wrapper (interface for Galaxy)


Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10
                     


License
-----------

 * MIT


