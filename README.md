# MetaSTAAR (Meta-analysis of variant-Set Test for Association using Annotation infoRmation)
This is an R package for performing MetaSTAAR procedure in whole-genome sequencing studies.
## Description
MetaSTAAR is an R package for performing Meta-analysis of variant-Set Test for Association using Annotation infoRmation (MetaSTAAR) procedure in whole-genome sequencing (WGS) studies.
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
The current version is 0.9.6 (December 28, 2021).
## License
This software is licensed under GPLv3.
