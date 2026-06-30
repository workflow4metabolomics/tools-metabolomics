# Table Merge

Metadata
-----------

 * **@name**: Convex dispersion
 * **@galaxyID**: convex_dispersion
 * **@version**: 0.1+galaxy1
 * **@authors**: Original code: Brice Mulot (PFEM - UNH - INRAE) - Maintainer: Etienne Jules (PFEM - UNH - INRAE - MetaboHUB)
 * **@init date**: 2025, July
 * **@main usage**: This tool displays convex hulls of metabolite intensities by injection order per batch.

 
Context
-----------

This tool generates convex hull plots for metabolite intensity data by injection order across different batches. This can be used to assess intensity values dispersion of ions on similar samples across batches and/or projects (QC, pool, reference materials).
 
Configuration
-----------

### Requirement:
 * R software: version = 4.1.2
 * r-ggplot2 = 3.3.5
 * r-optparse = 1.6.6
 * r-dplyr = 1.0.10

Technical description
-----------

Main files:

- plot_convex_hull.R: R function (core script)
- plot_convex_hull.xml: XML wrapper (interface for Galaxy)


Services provided
-----------

 * Help and support: https://community.france-bioinformatique.fr/c/workflow4metabolomics/10
                     


License
-----------

 * Cea Cnrs Inria Logiciel Libre License, version 2.1 (CECILL-2.1)


