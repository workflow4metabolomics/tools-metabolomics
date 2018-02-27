ANOVA for repeated measures statistics for Galaxy
=============================================

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io) [![Build Status](https://travis-ci.org/workflow4metabolomics/mixmodel4repeated_measures.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/mixmodel4repeated_measures)

Our project
-----------
The [Workflow4Metabolomics](http://workflow4metabolomics.org), W4M in short, is a French infrastructure offering software tool processing, analyzing and annotating metabolomics data. It is based on the Galaxy platform.


ANOVA for repeated measures
---------------------------

Analysis of variance for repeated measures using mixed model 


Galaxy
------
Galaxy is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 

Homepage: [https://galaxyproject.org/](https://galaxyproject.org/)


Dependencies using Conda
------------------------
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io)

The main recipe: [https://github.com/bioconda/bioconda-recipes/tree/master/recipes/r-lme4](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/r-lme4)

```
#To install miniconda2
#http://conda.pydata.org/miniconda.html
#To install the needed R library using conda:
conda install r-batch r-lme4 r-lmerTest r-viridisLite r-viridis r-Hmisc bioconductor-multtest
#To set an environment:
conda create -n mixmodel4repeated_measures r-batch r-lme4 r-lmerTest r-viridisLite r-viridis r-Hmisc bioconductor-multtest`
#To activate the environment:
. activate mixmodel4repeated_measures
```

[Conda](http://conda.pydata.org/) is package manager that among many other things can be used to manage Python packages.

Travis
------
[![Build Status](https://travis-ci.org/workflow4metabolomics/mixmodel4repeated_measures.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/mixmodel4repeated_measures)

Test and Deploy with Confidence. Easily sync your GitHub projects with Travis CI and you'll be testing your code in minutes!

Historic contributors
---------------------
 - Jean-François Martin - [French Metabolomics and Fluxomics Infrastructure (MetaboHUB)](http://www.metabohub.fr/en) - [MetaToul](http://www.metatoul.fr/)
 - Natacha Lenuzza - [French Metabolomics and Fluxomics Infrastructure (MetaboHUB)](http://www.metabohub.fr/en) - [CEA]
