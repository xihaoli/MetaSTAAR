#' Meta-analysis of STAAR (MetaSTAAR) procedure using omnibus test
#'
#' The \code{MetaSTAAR} function takes in the object from the merged summary statistics
#' and covariance files of each individual study and functional annotation data
#' to analyze the association between a quantitative/dichotomous phenotype and
#' a variant-set by using the meta-analysis of STAAR (MetaSTAAR) procedure.
#' For each variant-set, the MetaSTAAR-O p-value is a p-value from an omnibus test
#' that aggregated SKAT-MS(1,25), SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25),
#' and ACAT-V-MS(1,1) together with p-values of each test weighted by each annotation using Cauchy method.
#' @param obj_MetaSTAAR_merge an object from merging the summary statistics
#' and covariance files from each participating study, which is the output from
#' \code{\link{MetaSTAAR_merge}}.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p),
#' where p is the number of genetic variants in the variant-set.
#' Continuous scores should be given in PHRED score scale, where the PHRED score
#' of j-th variant is defined to be -10*log10(rank(-score_j)/total) across the genome. (Binary)
#' categorical scores should be taking values 0 or 1, where 1 is functional and 0 is
#' non-functional. If not provided, MetaSTAAR will perform the
#' SKAT-MS(1,25), SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25), ACAT-V-MS(1,1)
#' and ACAT-O-MS tests (default = NULL).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @return a list with the following members:
#' @return \code{num_variant}: the number of variants with combined minor allele frequency > 0 and less than
#' \code{rare_maf_cutoff} in the given variant-set that are used for performing the
#' variant-set test using MetaSTAAR.
#' @return \code{cMAC}: the combined cumulative minor allele count of variants with
#' combined minor allele frequency > 0 and less than \code{rare_maf_cutoff} in the given variant-set.
#' @return \code{results_MetaSTAAR_O}: the MetaSTAAR-O p-value that aggregated SKAT-MS(1,25),
#' SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25), and ACAT-V-MS(1,1) together
#' with p-values of each test weighted by each annotation using Cauchy method.
#' @return \code{results_ACAT_O_MS}: the ACAT-O-MS p-value that aggregated SKAT-MS(1,25),
#' SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25), and ACAT-V-MS(1,1) using Cauchy method.
#' @return \code{results_MetaSTAAR_S_1_25}: a vector of MetaSTAAR-S(1,25) p-values,
#' including SKAT-MS(1,25) p-value weighted by MAF, the SKAT-MS(1,25)
#' p-values weighted by each annotation, and a MetaSTAAR-S(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_S_1_1}: a vector of MetaSTAAR-S(1,1) p-values,
#' including SKAT-MS(1,1) p-value weighted by MAF, the SKAT-MS(1,1)
#' p-values weighted by each annotation, and a MetaSTAAR-S(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_B_1_25}: a vector of MetaSTAAR-B(1,25) p-values,
#' including Burden-MS(1,25) p-value weighted by MAF, the Burden-MS(1,25)
#' p-values weighted by each annotation, and a MetaSTAAR-B(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_B_1_1}: a vector of MetaSTAAR-B(1,1) p-values,
#' including Burden-MS(1,1) p-value weighted by MAF, the Burden-MS(1,1)
#' p-values weighted by each annotation, and a MetaSTAAR-B(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_A_1_25}: a vector of MetaSTAAR-A(1,25) p-values,
#' including ACAT-V-MS(1,25) p-value weighted by MAF, the ACAT-V-MS(1,25)
#' p-values weighted by each annotation, and a MetaSTAAR-A(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_A_1_1}: a vector of MetaSTAAR-A(1,1) p-values,
#' including ACAT-V-MS(1,1) p-value weighted by MAF, the ACAT-V-MS(1,1)
#' p-values weighted by each annotation, and a MetaSTAAR-A(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @references Li, X., et al. (2022). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics 52}(9), 969-983.
#' (\href{https://www.nature.com/articles/s41588-020-0676-4}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics 104}(3), 410-421.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0002929719300023}{pub})
#' @export

MetaSTAAR <- function(obj_MetaSTAAR_merge,annotation_phred=NULL,rv_num_cutoff=2){

  if(length(obj_MetaSTAAR_merge$U) == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & length(obj_MetaSTAAR_merge$U) != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  if(length(obj_MetaSTAAR_merge$U) >= rv_num_cutoff){
    MAF <- obj_MetaSTAAR_merge$info$MAF

    annotation_rank <- 1 - 10^(-annotation_phred/10)

    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(dim(annotation_phred)[2] == 0){
      ## Burden, SKAT, ACAT-V
      w_B <- w_S <- as.matrix(cbind(w_1,w_2))
      w_A <- as.matrix(cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_2^2/dbeta(MAF,0.5,0.5)^2))
    }else{
      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_1 <- cbind(w_1,w_B_1)
      w_B_2 <- annotation_rank*w_2
      w_B_2 <- cbind(w_2,w_B_2)
      w_B <- cbind(w_B_1,w_B_2)
      w_B <- as.matrix(w_B)

      ## SKAT
      w_S_1 <- sqrt(annotation_rank)*w_1
      w_S_1 <- cbind(w_1,w_S_1)
      w_S_2 <- sqrt(annotation_rank)*w_2
      w_S_2 <- cbind(w_2,w_S_2)
      w_S <- cbind(w_S_1,w_S_2)
      w_S <- as.matrix(w_S)

      ## ACAT-V
      w_A_1 <- annotation_rank*w_1^2/dbeta(MAF,0.5,0.5)^2
      w_A_1 <- cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_A_1)
      w_A_2 <- annotation_rank*w_2^2/dbeta(MAF,0.5,0.5)^2
      w_A_2 <- cbind(w_2^2/dbeta(MAF,0.5,0.5)^2,w_A_2)
      w_A <- cbind(w_A_1,w_A_2)
      w_A <- as.matrix(w_A)
    }

    pvalues <- MetaSTAAR_O_SMMAT(obj_MetaSTAAR_merge$U,obj_MetaSTAAR_merge$cov,
                                 weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                 mac=obj_MetaSTAAR_merge$info$MAC)

    num_variant <- length(obj_MetaSTAAR_merge$U)
    cMAC <- sum(obj_MetaSTAAR_merge$info$MAC)
    num_annotation <- dim(annotation_phred)[2]+1
    results_MetaSTAAR_O <- CCT(pvalues)
    results_ACAT_O_MS <- CCT(pvalues[c(1,num_annotation+1,2*num_annotation+1,3*num_annotation+1,4*num_annotation+1,5*num_annotation+1)])
    pvalues_MetaSTAAR_S_1_25 <- CCT(pvalues[1:num_annotation])
    pvalues_MetaSTAAR_S_1_1 <- CCT(pvalues[(num_annotation+1):(2*num_annotation)])
    pvalues_MetaSTAAR_B_1_25 <- CCT(pvalues[(2*num_annotation+1):(3*num_annotation)])
    pvalues_MetaSTAAR_B_1_1 <- CCT(pvalues[(3*num_annotation+1):(4*num_annotation)])
    pvalues_MetaSTAAR_A_1_25 <- CCT(pvalues[(4*num_annotation+1):(5*num_annotation)])
    pvalues_MetaSTAAR_A_1_1 <- CCT(pvalues[(5*num_annotation+1):(6*num_annotation)])

    results_MetaSTAAR_S_1_25 <- c(pvalues[1:num_annotation],pvalues_MetaSTAAR_S_1_25)
    results_MetaSTAAR_S_1_25 <- data.frame(t(results_MetaSTAAR_S_1_25))

    results_MetaSTAAR_S_1_1 <- c(pvalues[(num_annotation+1):(2*num_annotation)],pvalues_MetaSTAAR_S_1_1)
    results_MetaSTAAR_S_1_1 <- data.frame(t(results_MetaSTAAR_S_1_1))

    results_MetaSTAAR_B_1_25 <- c(pvalues[(2*num_annotation+1):(3*num_annotation)],pvalues_MetaSTAAR_B_1_25)
    results_MetaSTAAR_B_1_25 <- data.frame(t(results_MetaSTAAR_B_1_25))

    results_MetaSTAAR_B_1_1 <- c(pvalues[(3*num_annotation+1):(4*num_annotation)],pvalues_MetaSTAAR_B_1_1)
    results_MetaSTAAR_B_1_1 <- data.frame(t(results_MetaSTAAR_B_1_1))

    results_MetaSTAAR_A_1_25 <- c(pvalues[(4*num_annotation+1):(5*num_annotation)],pvalues_MetaSTAAR_A_1_25)
    results_MetaSTAAR_A_1_25 <- data.frame(t(results_MetaSTAAR_A_1_25))

    results_MetaSTAAR_A_1_1 <- c(pvalues[(5*num_annotation+1):(6*num_annotation)],pvalues_MetaSTAAR_A_1_1)
    results_MetaSTAAR_A_1_1 <- data.frame(t(results_MetaSTAAR_A_1_1))

    if(dim(annotation_phred)[2] == 0){
      colnames(results_MetaSTAAR_S_1_25) <- c("SKAT-MS(1,25)","MetaSTAAR-S(1,25)")
      colnames(results_MetaSTAAR_S_1_1) <- c("SKAT-MS(1,1)","MetaSTAAR-S(1,1)")
      colnames(results_MetaSTAAR_B_1_25) <- c("Burden-MS(1,25)","MetaSTAAR-B(1,25)")
      colnames(results_MetaSTAAR_B_1_1) <- c("Burden-MS(1,1)","MetaSTAAR-B(1,1)")
      colnames(results_MetaSTAAR_A_1_25) <- c("ACAT-V-MS(1,25)","MetaSTAAR-A(1,25)")
      colnames(results_MetaSTAAR_A_1_1) <- c("ACAT-V-MS(1,1)","MetaSTAAR-A(1,1)")
    }else{
      colnames(results_MetaSTAAR_S_1_25) <- c("SKAT-MS(1,25)",
                                          paste0("SKAT-MS(1,25)-",colnames(annotation_phred)),
                                          "MetaSTAAR-S(1,25)")
      colnames(results_MetaSTAAR_S_1_1) <- c("SKAT-MS(1,1)",
                                         paste0("SKAT-MS(1,1)-",colnames(annotation_phred)),
                                         "MetaSTAAR-S(1,1)")
      colnames(results_MetaSTAAR_B_1_25) <- c("Burden-MS(1,25)",
                                          paste0("Burden-MS(1,25)-",colnames(annotation_phred)),
                                          "MetaSTAAR-B(1,25)")
      colnames(results_MetaSTAAR_B_1_1) <- c("Burden-MS(1,1)",
                                         paste0("Burden-MS(1,1)-",colnames(annotation_phred)),
                                         "MetaSTAAR-B(1,1)")
      colnames(results_MetaSTAAR_A_1_25) <- c("ACAT-V-MS(1,25)",
                                          paste0("ACAT-V-MS(1,25)-",colnames(annotation_phred)),
                                          "MetaSTAAR-A(1,25)")
      colnames(results_MetaSTAAR_A_1_1) <- c("ACAT-V-MS(1,1)",
                                         paste0("ACAT-V-MS(1,1)-",colnames(annotation_phred)),
                                         "MetaSTAAR-A(1,1)")
    }

    return(list(num_variant = num_variant,
                cMAC = cMAC,
                results_MetaSTAAR_O = results_MetaSTAAR_O,
                results_ACAT_O_MS = results_ACAT_O_MS,
                results_MetaSTAAR_S_1_25 = results_MetaSTAAR_S_1_25,
                results_MetaSTAAR_S_1_1 = results_MetaSTAAR_S_1_1,
                results_MetaSTAAR_B_1_25 = results_MetaSTAAR_B_1_25,
                results_MetaSTAAR_B_1_1 = results_MetaSTAAR_B_1_1,
                results_MetaSTAAR_A_1_25 = results_MetaSTAAR_A_1_25,
                results_MetaSTAAR_A_1_1 = results_MetaSTAAR_A_1_1))
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

