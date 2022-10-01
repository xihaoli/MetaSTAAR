#' Generating sparse weighted covariance file using MetaSTAAR_worker
#'
#' The \code{MetaSTAAR_worker_cov} function takes in genotype, the object
#' from fitting the null model, and variant position to generate the sparse weighted
#' covariance file for the given variant-set as a rectangle format.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants. If the input genotype matrix
#' is sparse (e.g. \code{dgCMatrix} format), it is assumed that it has been flipped to represent
#' minor allele coding.
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples in the \code{\link{STAAR}} package.
#' @param cov_maf_cutoff a numeric value indicating the maximum minor allele frequency cutoff
#' under which the sparse weighted covariance file between variants is stored.
#' @param variant_pos a numeric vector of length p (listed in the same order as the columns
#' of \code{genotype}) indicating the position of the variants in the variant-set.
#' @param region_midpos a numeric value indicating the middle position of variant-set
#' by which the shorter edge of the rectangle is defined.
#' @param qc_label a vector of quality control status for each variant in \code{variant_pos}, where a pass variant
#' is labeled as "PASS". If \code{qc_label} is NULL, it is assumed that all variants are pass variants in the study (default = NULL).
#' @param segment.size a numeric value indicating the length of each segment of which
#' the sparse weighted covariance file is stored (default = 5e+05).
#' @param signif.digits an integer indicating the number of significant digits to be used
#' for storing the sparse weighted covariance file (default = 3).
#' @return \code{GTSinvG_rare}: the sparse matrix of all variants in the variant-set
#' whose minor allele frequency is below \code{cov_maf_cutoff} (the sparse weighted
#' covariance file), stored as a rectangle format.

MetaSTAAR_worker_cov <- function(genotype,obj_nullmodel,cov_maf_cutoff,
                                 variant_pos,region_midpos,qc_label=NULL,segment.size=5e5,signif.digits=3){

  if(class(genotype)[1] != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(cov_maf_cutoff < 0 | cov_maf_cutoff > 0.5){
    stop("cov_maf_cutoff should be a number between 0 and 0.5!")
  }

  if (cov_maf_cutoff == 0.5){
    cov_maf_cutoff <- 0.5 + 1e-16
  }

  if(dim(genotype)[2] != length(variant_pos)){
    stop(paste0("Dimensions don't match for genotype and variant_pos!"))
  }

  if(!is.null(qc_label) && length(variant_pos) != length(qc_label)){
    stop(paste0("Dimensions don't match for variant_pos and qc_label!"))
  }

  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    MAF <- colMeans(genotype)/2
    if (!is.null(qc_label)){
      RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0)&(qc_label=="PASS"))
    }else{
      RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0))
    }
    Geno_rare <- genotype[,RV_label,drop=FALSE]
    Geno_rare_1 <- Geno_rare[,variant_pos[RV_label] <= region_midpos,drop=FALSE]
    Geno_rare <- as(Geno_rare,"dgCMatrix")
    Geno_rare_1 <- as(Geno_rare_1,"dgCMatrix")
  }else{
    genotype <- matrix_flip(genotype)
    MAF <- genotype$MAF
    if (!is.null(qc_label)){
      RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0)&(qc_label=="PASS"))
    }else{
      RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0))
    }
    Geno_rare <- genotype$Geno[,RV_label,drop=FALSE]
    Geno_rare_1 <- Geno_rare[,variant_pos[RV_label] <= region_midpos,drop=FALSE]
    Geno_rare <- as(Geno_rare,"dgCMatrix")
    Geno_rare_1 <- as(Geno_rare_1,"dgCMatrix")
  }

  rm(genotype,MAF)
  gc()

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

  GTSinvG_rare <- t(obj_nullmodel$Sigma_i %*% Geno_rare_1) %*% Geno_rare
  rm(Geno_rare,Geno_rare_1)
  gc()

  GTSinvG_rare <- as(GTSinvG_rare,"dgTMatrix")
  row_pos <- (variant_pos[RV_label][variant_pos[RV_label] <= region_midpos])[(GTSinvG_rare@i)+1L]
  col_pos <- (variant_pos[RV_label])[(GTSinvG_rare@j)+1L]
  remove_ind <- (col_pos - row_pos > segment.size)
  rm(row_pos,col_pos)
  gc()
  remove_ind <- (GTSinvG_rare@j < GTSinvG_rare@i) | remove_ind # Be careful about multi-allelic issue
  GTSinvG_rare@i <- GTSinvG_rare@i[!remove_ind]
  GTSinvG_rare@j <- GTSinvG_rare@j[!remove_ind]
  GTSinvG_rare@x <- GTSinvG_rare@x[!remove_ind]
  rm(remove_ind)
  gc()
  GTSinvG_rare <- as(GTSinvG_rare,"dgCMatrix")

  ### Create a version of GTSinvG_rare with rounded significant digits
  GTSinvG_rare <- signif(GTSinvG_rare, digits = signif.digits)

  ### Add row/column names to the sparse matrix
  colnames(GTSinvG_rare) <- variant_pos[RV_label]
  rownames(GTSinvG_rare) <- variant_pos[RV_label][variant_pos[RV_label] <= region_midpos]

  return(GTSinvG_rare)
}

