% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetaSTAAR_worker_cov.R
\name{MetaSTAAR_worker_cov}
\alias{MetaSTAAR_worker_cov}
\title{Generating sparse weighted covariance file using MetaSTAARWorker}
\usage{
MetaSTAAR_worker_cov(
  genotype,
  obj_nullmodel,
  cov_maf_cutoff,
  variant_pos,
  region_midpos,
  qc_label = NULL,
  segment.size = 5e+05,
  signif.digits = 3
)
}
\arguments{
\item{genotype}{an n*p genotype matrix (dosage matrix) of the target sequence,
where n is the sample size and p is the number of genetic variants. If the input genotype matrix
is sparse (e.g. \code{dgCMatrix} format), it is assumed that it has been flipped to represent
minor allele coding.}

\item{obj_nullmodel}{an object from fitting the null model, which is the
output from either \code{\link{fit_null_glm}} function for unrelated samples or
\code{\link{fit_null_glmmkin}} function for related samples in the \code{\link{STAAR}} package.}

\item{cov_maf_cutoff}{a numeric value indicating the maximum minor allele frequency cutoff
under which the sparse weighted covariance file between variants is stored.}

\item{variant_pos}{a numeric vector of length p (listed in the same order as the columns
of \code{genotype}) indicating the position of the variants in the variant-set.}

\item{region_midpos}{a numeric value indicating the middle position of variant-set
by which the shorter edge of the rectangle is defined.}

\item{qc_label}{a vector of quality control status for each variant in \code{variant_pos}, where a pass variant
is labeled as "PASS". If \code{qc_label} is NULL, it is assumed that all variants are pass variants in the study (default = NULL).}

\item{segment.size}{a numeric value indicating the length of each segment of which
the sparse weighted covariance file is stored (default = 5e+05).}

\item{signif.digits}{an integer indicating the number of significant digits to be used
for storing the sparse weighted covariance file. If \code{signif.digits} is NULL,
it is assumed that no rounding will be performed (default = 3).}
}
\value{
\code{GTSinvG_rare}: the sparse matrix of all variants in the variant-set
whose minor allele frequency is below \code{cov_maf_cutoff} (the sparse weighted
covariance file), stored as a rectangle format.
}
\description{
The \code{MetaSTAAR_worker_cov} function takes in genotype, the object
from fitting the null model, and variant position to generate the sparse weighted
covariance file for the given variant-set as a rectangle format.
}
\references{
Li, X., et al. (2023). Powerful, scalable and resource-efficient
meta-analysis of rare variant associations in large whole genome sequencing studies.
\emph{Nature Genetics}, \emph{55}(1), 154-164.
(\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
}
