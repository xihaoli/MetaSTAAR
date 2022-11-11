#' Generating covariance file for conditional analysis using MetaSTAARWorker
#'
#' The \code{MetaSTAAR_worker_cov_cond} function takes in genotype, the genotype
#' of variants to be adjusted for in conditional analysis, the object
#' from fitting the null model, variant information and adjusted variant information (unique identifier)
#' to generate the conditional covariance file for the given variant-set,
#' adjusting for a given list of variants.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants. If the input genotype matrix
#' is sparse (e.g. \code{dgCMatrix} format), it is assumed that it has been flipped to represent
#' minor allele coding.
#' @param genotype_adj an n*p_adj genotype matrix (dosage matrix) of the target sequence, where n is
#' the sample size and p_adj is the number of genetic variants to be adjusted for in conditional
#' analysis (or a vector of a single variant with length n if p_adj is 1).
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples in the \code{\link{STAAR}} package.
#' @param variant_info a data frame or matrix of variant information (unique identifier)
#' with p rows (listed in the same order as the columns of \code{genotype}) and should contain
#' the following 4 columns: chromosome (chr), position (pos), reference allele (ref), and alternative allele (alt).
#' @param variant_adj_info a data frame or matrix of adjusted variant information (unique identifier)
#' with p_adj rows (listed in the same order as the rows of \code{genotype_adj}) and should contain
#' the following 4 columns: chromosome (chr), position (pos), reference allele (ref), and alternative allele (alt).
#' @return a list with the following members:
#' @return \code{GTPG_cond}: the covariance matrix between all variants in the variant-set (rows)
#' and all variants in the conditional variant-set (columns) (the covariance file
#' for conditional analysis).
#' @return \code{variant_info}: the data frame or matrix of variant information (unique identifier)
#' with p rows (listed in the same order as the rows of \code{GTPG_cond}) and 4 columns: chromosome (chr),
#' position (pos), reference allele (ref), and alternative allele (alt).
#' @return \code{variant_adj_info}: the data frame or matrix of adjusted variant information (unique identifier)
#' with p_adj rows (listed in the same order as the columns of \code{GTPG_cond}) and 4 columns: chromosome (chr),
#' position (pos), reference allele (ref), alternative allele (alt), score statistic (U), and variance (V).

MetaSTAAR_worker_cov_cond <- function(genotype,genotype_adj,obj_nullmodel,variant_info,variant_adj_info){

  if(class(genotype)[1] != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] != dim(variant_info)[1]){
    stop(paste0("Dimensions don't match for genotype and variant_info!"))
  }

  if(class(genotype_adj)[1] == "numeric" || class(genotype_adj)[1] == "integer"){
    genotype_adj <- matrix(genotype_adj, ncol=1)
  }

  if(!is.null(attr(class(genotype_adj), "package")) && attr(class(genotype_adj), "package") == "Matrix"){
    genotype_adj <- as.matrix(genotype_adj)
  }

  if(dim(genotype_adj)[2] != dim(variant_adj_info)[1]){
    stop(paste0("Dimensions don't match for genotype_adj and variant_adj_info!"))
  }

  if(dim(genotype)[1] != dim(genotype_adj)[1]){
    stop(paste0("Dimensions don't match for genotype and genotype_adj!"))
  }

  variant_info$index <- 1:dim(variant_info)[1]
  var.adj.index <- left_join(variant_adj_info,variant_info,by=c("chr"="chr",
                                                                "pos"="pos",
                                                                "ref"="ref",
                                                                "alt"="alt"))$index

  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype[,var.adj.index] <- 2 - genotype_adj
  }else{
    genotype <- matrix_flip(genotype)$Geno
    genotype[,var.adj.index] <- 2 - genotype_adj
    genotype <- as(genotype,"dgCMatrix")
  }

  genotype_adj <- 2 - genotype_adj # For faster computation

  if(obj_nullmodel$relatedness){
    if(!obj_nullmodel$sparse_kins){
      stop(paste0("Please use a sparse genetic relatedness matrix when fitting the null model!"))
    }
  }else{
    if(obj_nullmodel$family[1] == "binomial"){
      obj_nullmodel$Sigma_i <- Diagonal(x = obj_nullmodel$weights)
    }else if(obj_nullmodel$family[1] == "gaussian"){
      obj_nullmodel$Sigma_i <- Diagonal(length(obj_nullmodel$y)) / summary(obj_nullmodel)$dispersion
    }
    obj_nullmodel$Sigma_iX <- obj_nullmodel$Sigma_i %*% model.matrix(obj_nullmodel)
    obj_nullmodel$cov <- solve(t(model.matrix(obj_nullmodel)) %*% obj_nullmodel$Sigma_i %*% model.matrix(obj_nullmodel))
    obj_nullmodel$scaled.residuals <- (obj_nullmodel$y - obj_nullmodel$fitted.values) / summary(obj_nullmodel)$dispersion
  }

  GTSinvX <- as(t(genotype) %*% obj_nullmodel$Sigma_iX,"matrix")
  G_adjTSinvX <- as(t(genotype_adj) %*% obj_nullmodel$Sigma_iX,"matrix")
  GTPG_cond <- as(t(obj_nullmodel$Sigma_i %*% genotype) %*% genotype_adj,"matrix") - GTSinvX %*% obj_nullmodel$cov %*% t(G_adjTSinvX)
  V <- diag(t(genotype_adj) %*% obj_nullmodel$Sigma_i %*% genotype_adj - G_adjTSinvX %*% obj_nullmodel$cov %*% t(G_adjTSinvX))
  U <- as.vector(t(genotype_adj) %*% obj_nullmodel$scaled.residuals) # U is for alt allele count
  rm(genotype,genotype_adj)
  gc()

  return(list(GTPG_cond=GTPG_cond,
              variant_info=variant_info[,c("chr","pos","ref","alt")],
              variant_adj_info=cbind(variant_adj_info[,c("chr","pos","ref","alt")],U,V)))
}

