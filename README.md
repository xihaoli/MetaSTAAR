[![R build status](https://github.com/xihaoli/MetaSTAAR/workflows/R-CMD-check/badge.svg)](https://github.com/xihaoli/MetaSTAAR/actions)
[![Build Status](https://travis-ci.com/xihaoli/MetaSTAAR.svg?branch=main)](https://app.travis-ci.com/github/xihaoli/MetaSTAAR)
[![Build status](https://ci.appveyor.com/api/projects/status/jt95g3hy0y9rt0kg/branch/main?svg=true)](https://ci.appveyor.com/project/xihaoli/staarpipeline/branch/main)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# MetaSTAAR (Meta-analysis of variant-Set Test for Association using Annotation infoRmation)
This is an R package for performing MetaSTAAR procedure in whole-genome sequencing studies.
## Description
MetaSTAAR is an R package for performing Meta-analysis of variant-Set Test for Association using Annotation infoRmation (MetaSTAAR) procedure in whole-genome sequencing (WGS) studies.
## Workflow Overview
![MetaSTAAR_workflow](docs/MetaSTAAR_workflow.png)
## Prerequisites
<a href="https://www.r-project.org">R</a> (recommended version >= 3.5.1)

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.
## Dependencies
MetaSTAAR links to R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a>, <a href="https://cran.r-project.org/web/packages/RcppArmadillo/index.html">RcppArmadillo</a> and <a href="https://https://github.com/xihaoli/STAAR">STAAR</a>, and also imports R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a>, <a href="https://https://github.com/xihaoli/STAAR">STAAR</a>, <a href="https://cran.r-project.org/web/packages/Matrix/index.html">Matrix</a>, <a href="https://cran.r-project.org/web/packages/dplyr/index.html">dplyr</a>, <a href="https://cran.r-project.org/web/packages/expm/index.html">expm</a>, <a href="https://cran.r-project.org/web/packages/MASS/index.html">MASS</a>. These dependencies should be installed before installing MetaSTAAR.
## Installation
```
library(devtools)
devtools::install_github("xihaoli/MetaSTAAR",ref="main")
```
## Usage
Please see the <a href="docs/MetaSTAAR_manual.pdf">**MetaSTAAR** user manual</a> for detailed usage of MetaSTAAR package.
## Version
The current version is 0.9.6.1 (September 30, 2022).
## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
