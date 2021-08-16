#' Meta-analysis of score test for individual variants in a given variant-set
#'
#' The \code{Indiv_Score_Test_Region_meta} function takes in the object from the
#' merged summary statistics and covariance files of each individual study
#' to analyze the associations between a quantitative/dichotomous phenotype and
#' all individual variants in a given variant-set by using the meta-analysis of score test.
#' @param obj_MetaSTAAR_merge an object from merging the summary statistics
#' and covariance files from each individual study, which is the output from
#' \code{\link{MetaSTAAR_merge}}.
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @return a data frame with p rows corresponding to the p genetic variants in the given variant-set
#' and three columns: \code{Score} (the score test statistic), \code{SE} (the standard error associated
#' with the score test statistic), and \code{pvalue} (the score test p-value).
#' If a variant in the given variant-set has standard error equal to 0, the p-value will be set as 1.
#' @export

Indiv_Score_Test_Region_meta <- function(obj_MetaSTAAR_merge,rv_num_cutoff=2){

  if(length(obj_MetaSTAAR_merge$U) == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  results <- data.frame(Score = rep(NA, length(obj_MetaSTAAR_merge$U)),
                        SE = rep(NA, length(obj_MetaSTAAR_merge$U)),
                        pvalue = rep(NA, length(obj_MetaSTAAR_merge$U)))

  if(length(obj_MetaSTAAR_merge$U) >= rv_num_cutoff){
    results[,] <- do.call(cbind,Indiv_Score_Test_meta(obj_MetaSTAAR_merge$U,
                                                      diag(obj_MetaSTAAR_merge$cov)))

    return(results)
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

