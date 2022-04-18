#' Effect size and standard error estimates of meta-analysis of burden test for a given variant-set
#'
#' The \code{Burden_Effect_Size_meta} function takes in the object from the merged summary statistics
#' and covariance files of each individual study and functional annotation data
#' to calculate the effect size and standard error estimates of meta-analysis of burden test for a given variant-set.
#' @param obj_MetaSTAAR_merge an object from merging the summary statistics
#' and covariance files from each individual study, which is the output from
#' \code{\link{MetaSTAAR_merge}}.
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @return a list with the following members:
#' @return \code{num_variant}: the number of variants with combined minor allele frequency > 0 and less than
#' \code{rare_maf_cutoff} in the given variant-set that are used for performing the
#' variant-set test using MetaSTAAR.
#' @return \code{cMAC}: the combined cumulative minor allele count of variants with
#' combined minor allele frequency > 0 and less than \code{rare_maf_cutoff} in the given variant-set.
#' @return \code{Burden_Score_Stat}: the score statistic of Burden-MS(1,1) for the given variant-set.
#' @return \code{Burden_SE_Score}: the standard error of \code{Burden_Score_Stat} for the given variant-set.
#' @return \code{Burden_pvalue}: the Burden-MS(1,1) p-value for the given variant-set.
#' @return \code{Burden_Est}: the effect size estimate of Burden-MS(1,1) for the given variant-set.
#' @return \code{Burden_SE_Est}: the standard error estimate of \code{Burden_Est} for the given variant-set.
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics 52}(9), 969-983.
#' (\href{https://www.nature.com/articles/s41588-020-0676-4}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics 104}(3), 410-421.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0002929719300023}{pub})
#' @export

Burden_Effect_Size_meta <- function(obj_MetaSTAAR_merge,rv_num_cutoff=2){

  if(length(obj_MetaSTAAR_merge$U) == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  if(length(obj_MetaSTAAR_merge$U) >= rv_num_cutoff){
    num_variant <- length(obj_MetaSTAAR_merge$U)
    cMAC <- sum(obj_MetaSTAAR_merge$info$MAC)

    Burden_Score_Stat <- sum(obj_MetaSTAAR_merge$U)
    Burden_Var <- sum(obj_MetaSTAAR_merge$cov)
    Burden_test_chisq <- Burden_Score_Stat^2/Burden_Var
    Burden_pvalue <- pchisq(Burden_test_chisq,df=1,lower.tail=FALSE)

    return(list(num_variant = num_variant,
                cMAC = cMAC,
                Burden_Score_Stat = Burden_Score_Stat,
                Burden_SE_Score = sqrt(Burden_Var),
                Burden_pvalue = Burden_pvalue,
                Burden_Est = Burden_Score_Stat/Burden_Var,
                Burden_SE_Est = 1/sqrt(Burden_Var)))
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

