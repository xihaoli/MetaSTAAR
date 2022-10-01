#' Meta-analysis of STAAR (MetaSTAAR) procedure for conditional analysis using omnibus test
#'
#' The \code{MetaSTAAR_cond} function takes in the object from the merged conditional summary statistics
#' and covariance files of each participating study and functional annotation data
#' to analyze the conditional association between a quantitative/dichotomous phenotype and
#' a variant-set by using the meta-analysis of STAAR (MetaSTAAR) procedure,
#' adjusting for a given list of variants.
#' For each variant-set, the MetaSTAAR-O p-value is a p-value from an omnibus test that aggregated
#' conditional SKAT-MS(1,25), SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25), and ACAT-V-MS(1,1)
#' together with conditional p-values of each test weighted by each annotation using Cauchy method.
#' @param obj_MetaSTAAR_merge_cond an object from merging the conditional summary statistics
#' and covariance files from each participating study, adjusting for a given list of variants,
#' which is the output from \code{\link{MetaSTAAR_merge_cond}}.
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
#' @return \code{results_MetaSTAAR_O_cond}: the conditional MetaSTAAR-O p-value that aggregated conditional
#' SKAT-MS(1,25), SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25), and ACAT-V-MS(1,1) together
#' with conditional p-values of each test weighted by each annotation using Cauchy method.
#' @return \code{results_ACAT_O_MS_cond}: the conditional ACAT-O-MS p-value that aggregated conditional
#' SKAT-MS(1,25), SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25), and ACAT-V-MS(1,1) using Cauchy method.
#' @return \code{results_MetaSTAAR_S_1_25_cond}: a vector of conditional MetaSTAAR-S(1,25) p-values,
#' including conditional SKAT-MS(1,25) p-value weighted by MAF, the conditional SKAT-MS(1,25)
#' p-values weighted by each annotation, and a conditional MetaSTAAR-S(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_S_1_1_cond}: a vector of conditional MetaSTAAR-S(1,1) p-values,
#' including conditional SKAT-MS(1,1) p-value weighted by MAF, the conditional SKAT-MS(1,1)
#' p-values weighted by each annotation, and a conditional MetaSTAAR-S(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_B_1_25_cond}: a vector of conditional MetaSTAAR-B(1,25) p-values,
#' including conditional Burden-MS(1,25) p-value weighted by MAF, the conditional Burden-MS(1,25)
#' p-values weighted by each annotation, and a conditional MetaSTAAR-B(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_B_1_1_cond}: a vector of conditional MetaSTAAR-B(1,1) p-values,
#' including conditional Burden-MS(1,1) p-value weighted by MAF, the conditional Burden-MS(1,1)
#' p-values weighted by each annotation, and a conditional MetaSTAAR-B(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_A_1_25_cond}: a vector of conditional MetaSTAAR-A(1,25) p-values,
#' including conditional ACAT-V-MS(1,25) p-value weighted by MAF, the conditional ACAT-V-MS(1,25)
#' p-values weighted by each annotation, and a conditional MetaSTAAR-A(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_MetaSTAAR_A_1_1_cond}: a vector of conditional MetaSTAAR-A(1,1) p-values,
#' including conditional ACAT-V-MS(1,1) p-value weighted by MAF, the conditional ACAT-V-MS(1,1)
#' p-values weighted by each annotation, and a conditional MetaSTAAR-A(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics 52}(9), 969-983.
#' (\href{https://www.nature.com/articles/s41588-020-0676-4}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics 104}(3), 410-421.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0002929719300023}{pub})
#' @export

MetaSTAAR_cond <- function(obj_MetaSTAAR_merge_cond,annotation_phred=NULL,rv_num_cutoff=2){

  if(length(obj_MetaSTAAR_merge_cond$U_cond) == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & length(obj_MetaSTAAR_merge_cond$U_cond) != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  if(length(obj_MetaSTAAR_merge_cond$U_cond) >= rv_num_cutoff){
    MAF <- obj_MetaSTAAR_merge_cond$info$MAF

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

    pvalues <- MetaSTAAR_O_SMMAT(obj_MetaSTAAR_merge_cond$U_cond,obj_MetaSTAAR_merge_cond$cov_cond,
                                 weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                 mac=obj_MetaSTAAR_merge_cond$info$MAC)

    num_variant <- length(obj_MetaSTAAR_merge_cond$U_cond)
    cMAC <- sum(obj_MetaSTAAR_merge_cond$info$MAC)
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
                results_MetaSTAAR_O_cond = results_MetaSTAAR_O,
                results_ACAT_O_MS_cond = results_ACAT_O_MS,
                results_MetaSTAAR_S_1_25_cond = results_MetaSTAAR_S_1_25,
                results_MetaSTAAR_S_1_1_cond = results_MetaSTAAR_S_1_1,
                results_MetaSTAAR_B_1_25_cond = results_MetaSTAAR_B_1_25,
                results_MetaSTAAR_B_1_1_cond = results_MetaSTAAR_B_1_1,
                results_MetaSTAAR_A_1_25_cond = results_MetaSTAAR_A_1_25,
                results_MetaSTAAR_A_1_1_cond = results_MetaSTAAR_A_1_1))
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

