[![R build status](https://github.com/xihaoli/MetaSTAAR/workflows/R-CMD-check/badge.svg)](https://github.com/xihaoli/MetaSTAAR/actions)
[![Build status](https://ci.appveyor.com/api/projects/status/jt95g3hy0y9rt0kg/branch/main?svg=true)](https://ci.appveyor.com/project/xihaoli/metastaar/branch/main)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# MetaSTAAR (Meta-analysis of variant-Set Test for Association using Annotation infoRmation)
This is an R package for performing MetaSTAAR procedure in whole genome sequencing studies. A lightweight package implementing MetaSTAAR with a pipeline for performing functionally-informed meta-analysis of sequencing studies is available in <a href="https://github.com/li-lab-genetics/MetaSTAARlite">**MetaSTAARlite**</a> [<a href="https://github.com/li-lab-genetics/MetaSTAARlite-tutorial">**Tutorial**</a>].
## Description
MetaSTAAR is an R package for performing Meta-analysis of variant-Set Test for Association using Annotation infoRmation (MetaSTAAR) procedure in whole genome sequencing (WGS) studies. MetaSTAAR enables functionally-informed rare variant meta-analysis of large WGS studies using an efficient, sparse matrix approach for storing summary statistic, while protecting data privacy of study participants and avoiding sharing subject-level data. MetaSTAAR accounts for relatedness and population structure of continuous and dichotomous traits, and boosts the power of rare variant meta-analysis by incorporating multiple variant functional annotations.
## Workflow Overview
![MetaSTAAR_workflow](docs/MetaSTAAR_workflow.jpg)
## Prerequisites
<a href="https://www.r-project.org">R</a> (recommended version >= 3.5.1)

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.
## Dependencies
MetaSTAAR links to R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a>, <a href="https://cran.r-project.org/web/packages/RcppArmadillo/index.html">RcppArmadillo</a> and <a href="https://github.com/xihaoli/STAAR">STAAR</a>, and also imports R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a>, <a href="https://github.com/xihaoli/STAAR">STAAR</a>, <a href="https://cran.r-project.org/web/packages/Matrix/index.html">Matrix</a>, <a href="https://cran.r-project.org/web/packages/dplyr/index.html">dplyr</a>, <a href="https://cran.r-project.org/web/packages/expm/index.html">expm</a>, <a href="https://cran.r-project.org/web/packages/MASS/index.html">MASS</a>. These dependencies should be installed before installing MetaSTAAR.
## Installation
```
library(devtools)
devtools::install_github("xihaoli/MetaSTAAR",ref="main")
```
## Usage
Please see the <a href="docs/MetaSTAAR_manual.pdf">**MetaSTAAR** user manual</a> for detailed usage of MetaSTAAR package. The scripts used to generate results in the <a href="https://doi.org/10.1038/s41588-022-01225-6">manuscript</a> are available on <a href="https://doi.org/10.5281/zenodo.6668274">_Zenodo_</a>.
## Data Availability
The whole-genome functional annotation data assembled from a variety of sources and the precomputed annotation principal components are available at the [Functional Annotation of Variant - Online Resource (FAVOR)](https://favor.genohub.org) site and [FAVOR Essential Database](https://doi.org/10.7910/DVN/1VGTJI).
## Version
The current version is 0.9.6.3 (February 5, 2024).
## Citation
If you use **MetaSTAAR** for your work, please cite:

Xihao Li, Corbin Quick, Hufeng Zhou, Sheila M. Gaynor, Yaowu Liu, Han Chen, Margaret Sunitha Selvaraj, Ryan Sun, Rounak Dey, Donna K. Arnett, Lawrence F. Bielak, Joshua C. Bis, John Blangero, Eric Boerwinkle, Donald W. Bowden, Jennifer A. Brody, Brian E. Cade, Adolfo Correa, L. Adrienne Cupples, Joanne E. Curran, Paul S. de Vries, Ravindranath Duggirala, Barry I. Freedman, Harald H. H. Göring, Xiuqing Guo, Jeffrey Haessler, Rita R. Kalyani, Charles Kooperberg, Brian G. Kral, Leslie A. Lange, Ani Manichaikul, Lisa W. Martin, Stephen T. McGarvey, Braxton D. Mitchell, May E. Montasser, Alanna C. Morrison, Take Naseri, Jeffrey R. O'Connell, Nicholette D. Palmer, Patricia A. Peyser, Bruce M. Psaty, Laura M. Raffield, Susan Redline, Alexander P. Reiner, Muagututi’a Sefuiva Reupena, Kenneth M. Rice, Stephen S. Rich, Colleen M. Sitlani, Jennifer A. Smith, Kent D. Taylor, Ramachandran S. Vasan, Cristen J. Willer, James G. Wilson, Lisa R. Yanek, Wei Zhao, NHLBI Trans-Omics for Precision Medicine (TOPMed) Consortium, TOPMed Lipids Working Group, Jerome I. Rotter, Pradeep Natarajan, Gina M. Peloso, Zilin Li, & Xihong Lin. (2023). **Powerful, scalable and resource-efficient meta-analysis of rare variant associations in large whole genome sequencing studies**. _Nature Genetics_, _55_(1), 154-164. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/36564505">36564505</a>. PMCID: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10084891/">PMC10084891</a>. DOI: <a href="https://doi.org/10.1038/s41588-022-01225-6">10.1038/s41588-022-01225-6</a>.
## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
