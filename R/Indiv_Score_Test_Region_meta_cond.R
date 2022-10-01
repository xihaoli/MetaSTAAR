#' Meta-analysis of conditional score test for individual variants in a given variant-set
#'
#' The \code{Indiv_Score_Test_Region_meta_cond} function takes in the object from the
#' merged conditional summary statistics and covariance files of each participating study to analyze
#' the conditional associations between a quantitative/dichotomous phenotype and
#' all individual variants in a given variant-set by using the meta-analysis of score test,
#' adjusting for a given list of variants.
#' @param obj_MetaSTAAR_merge_cond an object from merging the conditional summary statistics
#' and covariance files from each participating study, adjusting for a given list of variants,
#' which is the output from \code{\link{MetaSTAAR_merge_cond}}.
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @return a data frame with p rows corresponding to the p genetic variants in the given variant-set
#' and three columns: \code{Score_cond} (the conditional score test statistic),
#' \code{SE_cond} (the standard error associated with the conditional score test statistic),
#' and \code{pvalue_cond} (the conditional score test p-value). If a variant in
#' the given variant-set has standard error equal to 0, the p-value will be set as 1.
#' @export

Indiv_Score_Test_Region_meta_cond <- function(obj_MetaSTAAR_merge_cond,rv_num_cutoff=2){

  if(length(obj_MetaSTAAR_merge_cond$U_cond) == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  results <- data.frame(Score_cond = rep(NA, length(obj_MetaSTAAR_merge_cond$U_cond)),
                        SE_cond = rep(NA, length(obj_MetaSTAAR_merge_cond$U_cond)),
                        pvalue_cond = rep(NA, length(obj_MetaSTAAR_merge_cond$U_cond)))

  if(length(obj_MetaSTAAR_merge_cond$U_cond) >= rv_num_cutoff){
    results[,] <- do.call(cbind,Indiv_Score_Test_meta(obj_MetaSTAAR_merge_cond$U_cond,
                                                      diag(obj_MetaSTAAR_merge_cond$cov_cond)))

    return(results)
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

