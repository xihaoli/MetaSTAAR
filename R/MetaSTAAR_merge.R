#' The preliminary data manipulation step for \code{MetaSTAAR}
#'
#' The \code{MetaSTAAR_merge} function takes in the summary statistics file and the sparse weighted
#' covariance file (the output from \code{\link{MetaSTAAR_worker_sumstat}} and \code{\link{MetaSTAAR_worker_cov}})
#' from each participating study and performs the preliminary data manipulation step
#' by merging them into a single unified summary statistics file and a covariance file, respectively.
#' @param chr a numeric value indicating the chromosome of the genetic region of interest.
#' @param start.loc a numeric value indicating the starting location (position) of the
#' genetic region of interest.
#' @param end.loc a numeric value indicating the ending location (position) of the
#' genetic region of interest.
#' @param study.names a character vector containing the name of each participating
#' study in the meta-analysis.
#' @param sample.sizes a numeric vector with the length of \code{study.names}
#' indicating the sample size of each study.
#' @param sumstat.dir a character vector containing the directories of the study-specific summary statistics file folders.
#' @param cov.dir a character vector containing the directories of the study-specific sparse weighted
#' covariance file folders.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param cov_maf_cutoff a numeric vector with the length of \code{study.names}
#' indicating the maximum minor allele frequency cutoffs under which the sparse weighted
#' covariance files between variants are stored.
#' @param trait a character value indicating the underlying trait of interest for
#' the meta-analysis.
#' @param segment.size a numeric value indicating the length of each segment of which
#' the summary statistics and sparse weighted covariance files are stored.
#' Note that the input value should be aligned with the input values of
#' \code{\link{MetaSTAAR_worker_cov}} (default = 5e+05).
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}
#' specified in \code{\link{MetaSTAAR_worker_sumstat}} and \code{\link{MetaSTAAR_worker_cov}}.
#' If \code{check_qc_label} is FALSE, it is assumed that no variant will be dropped (default = FALSE).
#' @return a list with the following members:
#' @return \code{info}: the merged data frame of all variants in the genetic region
#' of interest whose combined minor allele frequency is below \code{rare_maf_cutoff}, including the
#' following information (listed in the same order as \code{U} and the rows/columns of \code{cov}):
#' chromosome (chr), position (pos), reference allele (ref), alternative allele (alt),
#' combined minor allele count (MAC), and combined minor allele frequency (MAF).
#' @return \code{U} the merged score statistics vector of all variants in the genetic region
#' of interest whose combined minor allele frequency is below \code{rare_maf_cutoff}.
#' @return \code{cov} the merged covariance matrix of all variants in the genetic region
#' of interest whose combined minor allele frequency is below \code{rare_maf_cutoff}.
#' @references Li, X., et al. (2023). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}, \emph{55}(1), 154-164.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @export

MetaSTAAR_merge <- function(chr,start.loc,end.loc,study.names,sample.sizes,sumstat.dir,cov.dir,
                            rare_maf_cutoff=0.01,cov_maf_cutoff,trait,segment.size=5e5,check_qc_label=FALSE){

  segment <- floor((start.loc - 1) / segment.size) + 1
  cov_maf_cutoff[cov_maf_cutoff == 0.5] <- 0.5 + 1e-16
  if (end.loc <= segment * segment.size) {
    ### summary statistics
    sumstat.files <- paste0(sumstat.dir,"/summary.stat.",trait,".",study.names,".chr",chr,".segment",segment,".Rdata")
    sumstat.list <- lapply(sumstat.files, function(x) {
      load(file = x)
      get(ls()[ls()!= "summary.stat"])
    })
    sumstat.list <- lapply(sumstat.list, function(x) {
      if (!(is.null(x)) && !("qc_label" %in% colnames(x))){
        data.frame(x[,1:4],qc_label="PASS",x[,5:dim(x)[2]],stringsAsFactors = FALSE)
      }else{
        x
      }
    })
    position.index <- mapply(function(x,y) {
      (x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),"pos"]>=start.loc)&(x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),"pos"]<=end.loc)
    }, x = sumstat.list, y = cov_maf_cutoff, SIMPLIFY = FALSE)
    sumstat.list <- lapply(sumstat.list, function(x) {
      x[(x$pos>=start.loc)&(x$pos<=end.loc),]
    })
    sumstat.varid.list <- mapply(function(x,y) {
      x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),1:4]
    }, x = sumstat.list, y = cov_maf_cutoff, SIMPLIFY = FALSE)
    sumstat.varid.merge <- do.call("rbind",sumstat.varid.list)
    sumstat.varid.nodup <- sumstat.varid.merge[!duplicated(sumstat.varid.merge),]
    if (is.null(sumstat.varid.nodup)) {
      return(NULL)
    }else if (dim(sumstat.varid.nodup)[1] == 0) {
      return(NULL)
    }
    sumstat.merge.list <- lapply(sumstat.list, function(x) {
      if (is.null(x)) {
        cbind(sumstat.varid.nodup,qc_label=NA,alt_AC=NA,MAC=NA,MAF=NA,N=NA,U=NA,V=NA,"1"=NA)
      }else {
        left_join(sumstat.varid.nodup,x,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))
      }
    })
    sumstat.merge.list <- mapply(function(x,y) {
      x[is.na(x[,"N"]),"N"] <- y
      x[is.na(x[,"qc_label"]),"qc_label"] <- "PASS"
      x[is.na(x)] <- 0
      return(x)
    }, x = sumstat.merge.list, y = sample.sizes, SIMPLIFY = FALSE)
    sumstat.varid.nodup$index <- 1:dim(sumstat.varid.nodup)[1]
    sumstat.index.list <- lapply(sumstat.varid.list, function(x) {
      if (is.null(x)) {
        integer(0)
      }else{
        left_join(x,sumstat.varid.nodup,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))$index
      }
    })
    rm(sumstat.files,sumstat.list,sumstat.varid.list,sumstat.varid.merge)
    gc()

    ### covariance matrices
    cov.files <- paste0(cov.dir,"/GTSinvG.rare.",trait,".",study.names,".chr",chr,".segment",segment,".Rdata")
    cov.list <- lapply(cov.files, function(x) {
      load(file = x)
      get(ls()[ls()!= "cov"])
    })
    cov.list <- lapply(cov.list, function(x) {
      if (is.null(x)) {
        as(matrix(nrow=0,ncol=0),"dgCMatrix")
      }else if (dim(x)[1] == 0) {
        as(matrix(nrow=0,ncol=0),"dgCMatrix")
      }else {
        x[,1:dim(x)[1]]
      }
    })
    cov.list <- mapply(function(x,y) {
      if (length(x) == 1) {
        x <- as(matrix(x),"dgCMatrix")
      }
      x[y,y]
    }, x = cov.list, y = position.index, SIMPLIFY = FALSE)

    cov.null <- matrix(0, nrow = dim(sumstat.varid.nodup)[1], ncol = dim(sumstat.varid.nodup)[1])
    cov.merge.list <- mapply(function(x,y) {
      if (length(x) == 1) {
        # in this case y is a scalar
        cov.null[x,x] <- y
      }else if (length(x) > 1) {
        cov.null[x,x] <- as.matrix(forceSymmetric(y, uplo = "U"))
      }
      return(cov.null)
    }, x = sumstat.index.list, y = cov.list, SIMPLIFY = FALSE)

  }else if (end.loc <= (segment + 1) * segment.size) {
    ### summary statistics
    sumstat.files1 <- paste0(sumstat.dir,"/summary.stat.",trait,".",study.names,".chr",chr,".segment",segment,".Rdata")
    sumstat.list1 <- lapply(sumstat.files1, function(x) {
      load(file = x)
      get(ls()[ls()!= "summary.stat"])
    })
    sumstat.list1 <- lapply(sumstat.list1, function(x) {
      if (!(is.null(x)) && !("qc_label" %in% colnames(x))){
        data.frame(x[,1:4],qc_label="PASS",x[,5:dim(x)[2]],stringsAsFactors = FALSE)
      }else{
        x
      }
    })
    position.index1 <- mapply(function(x,y) {
      (x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),"pos"]>=start.loc)&(x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),"pos"]<=end.loc)
    }, x = sumstat.list1, y = cov_maf_cutoff, SIMPLIFY = FALSE)
    sumstat.list1 <- lapply(sumstat.list1, function(x) {
      x[(x$pos>=start.loc)&(x$pos<=end.loc),]
    })
    sumstat.varid.list1 <- mapply(function(x,y) {
      x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),1:4]
    }, x = sumstat.list1, y = cov_maf_cutoff, SIMPLIFY = FALSE)

    sumstat.files2 <- paste0(sumstat.dir,"/summary.stat.",trait,".",study.names,".chr",chr,".segment",segment+1,".Rdata")
    sumstat.list2 <- lapply(sumstat.files2, function(x) {
      load(file = x)
      get(ls()[ls()!= "summary.stat"])
    })
    sumstat.list2 <- lapply(sumstat.list2, function(x) {
      if (!(is.null(x)) && !("qc_label" %in% colnames(x))){
        data.frame(x[,1:4],qc_label="PASS",x[,5:dim(x)[2]],stringsAsFactors = FALSE)
      }else{
        x
      }
    })
    position.index2 <- mapply(function(x,y) {
      (x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),"pos"]>=start.loc)&(x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),"pos"]<=end.loc)
    }, x = sumstat.list2, y = cov_maf_cutoff, SIMPLIFY = FALSE)
    sumstat.list2 <- lapply(sumstat.list2, function(x) {
      x[(x$pos>=start.loc)&(x$pos<=end.loc),]
    })
    sumstat.varid.list2 <- mapply(function(x,y) {
      x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),1:4]
    }, x = sumstat.list2, y = cov_maf_cutoff, SIMPLIFY = FALSE)

    sumstat.list <- mapply(function(x,y) {
      rbind(x,y)
    }, x = sumstat.list1, y = sumstat.list2, SIMPLIFY = FALSE)
    sumstat.varid.list <- mapply(function(x,y) {
      rbind(x,y)
    }, x = sumstat.varid.list1, y = sumstat.varid.list2, SIMPLIFY = FALSE)
    sumstat.varid.merge <- do.call("rbind",sumstat.varid.list)
    sumstat.varid.nodup <- sumstat.varid.merge[!duplicated(sumstat.varid.merge),]
    if (is.null(sumstat.varid.nodup)) {
      return(NULL)
    }else if (dim(sumstat.varid.nodup)[1] == 0) {
      return(NULL)
    }
    sumstat.merge.list <- lapply(sumstat.list, function(x) {
      if (is.null(x)) {
        cbind(sumstat.varid.nodup,qc_label=NA,alt_AC=NA,MAC=NA,MAF=NA,N=NA,U=NA,V=NA,"1"=NA)
      }else {
        left_join(sumstat.varid.nodup,x,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))
      }
    })
    sumstat.merge.list <- mapply(function(x,y) {
      x[is.na(x[,"N"]),"N"] <- y
      x[is.na(x[,"qc_label"]),"qc_label"] <- "PASS"
      x[is.na(x)] <- 0
      return(x)
    }, x = sumstat.merge.list, y = sample.sizes, SIMPLIFY = FALSE)
    sumstat.varid.nodup$index <- 1:dim(sumstat.varid.nodup)[1]
    sumstat.index.list1 <- lapply(sumstat.varid.list1, function(x) {
      if (is.null(x)) {
        integer(0)
      }else{
        left_join(x,sumstat.varid.nodup,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))$index
      }
    })
    sumstat.index.list2 <- lapply(sumstat.varid.list2, function(x) {
      if (is.null(x)) {
        integer(0)
      }else{
        left_join(x,sumstat.varid.nodup,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))$index
      }
    })
    rm(sumstat.files1,sumstat.files2,sumstat.list,sumstat.list1,sumstat.list2,
       sumstat.varid.list,sumstat.varid.list1,sumstat.varid.list2,
       sumstat.varid.merge)
    gc()

    ### covariance matrices
    cov.files1 <- paste0(cov.dir,"/GTSinvG.rare.",trait,".",study.names,".chr",chr,".segment",segment,".Rdata")
    cov.list1 <- lapply(cov.files1, function(x) {
      load(file = x)
      get(ls()[ls()!= "cov"])
    })
    cov.list1 <- mapply(function(x,y,z) {
      if (is.null(x)) {
        x <- as(matrix(nrow=length(y),ncol=(length(y)+length(z))),"dgCMatrix")
      }
      x[y,c(y,z)]
    }, x = cov.list1, y = position.index1, z = position.index2, SIMPLIFY = FALSE)

    cov.files2 <- paste0(cov.dir,"/GTSinvG.rare.",trait,".",study.names,".chr",chr,".segment",segment+1,".Rdata")
    cov.list2 <- lapply(cov.files2, function(x) {
      load(file = x)
      get(ls()[ls()!= "cov"])
    })
    cov.list2 <- lapply(cov.list2, function(x) {
      if (is.null(x)) {
        as(matrix(nrow=0,ncol=0),"dgCMatrix")
      }else if (dim(x)[1] == 0) {
        as(matrix(nrow=0,ncol=0),"dgCMatrix")
      }else {
        x[,1:dim(x)[1]]
      }
    })
    cov.list2 <- mapply(function(x,y) {
      if (length(x) == 1) {
        x <- as(matrix(x),"dgCMatrix")
      }
      x[y,y]
    }, x = cov.list2, y = position.index2, SIMPLIFY = FALSE)

    cov.null <- matrix(0, nrow = dim(sumstat.varid.nodup)[1], ncol = dim(sumstat.varid.nodup)[1])
    cov.merge.list <- mapply(function(x,y,z,w) {
      cov.null[x,c(x,y)] <- as.matrix(z)
      cov.null[y,y] <- as.matrix(w)
      if (length(c(x,y)) > 1) {
        cov.null[c(x,y),c(x,y)] <- as.matrix(forceSymmetric(cov.null[c(x,y),c(x,y)], uplo = "U"))
      }
      return(cov.null)
    }, x = sumstat.index.list1, y = sumstat.index.list2, z = cov.list1, w = cov.list2, SIMPLIFY = FALSE)

  }

  ### select rare variant based on the input cutoff
  alt_AC.merge <- as.integer(Reduce("+",lapply(sumstat.merge.list, function(x) {x$alt_AC})))
  N.merge.nonzero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC != 0)}))
  alt_AF.merge.nonzero <- alt_AC.merge / (2 * N.merge.nonzero)
  N.merge.zero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC == 0)}))
  alt_AC.merge <- alt_AC.merge + (alt_AF.merge.nonzero > 0.5) * (2 * N.merge.zero)
  N.merge <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N}))
  MAC.merge <- pmin(alt_AC.merge, 2 * N.merge - alt_AC.merge)
  MAF.merge <- MAC.merge / (2 * N.merge)
  rv.index <- (MAF.merge<rare_maf_cutoff) & Reduce("*",mapply(function(x,y) {
    x$MAF<y
  }, x = sumstat.merge.list, y = cov_maf_cutoff, SIMPLIFY = FALSE))
  if (check_qc_label){
    rv.index <- rv.index & Reduce("*",lapply(sumstat.merge.list, function(x) {
      x$qc_label=="PASS"
    }))
  }

  info <- cbind(sumstat.varid.nodup[,c("chr","pos","ref","alt")],
                MAC=MAC.merge,MAF=MAF.merge)[rv.index,]

  U.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
    x <- x[rv.index,]
    return(x$U)
  }))

  cov.merge <- Reduce("+", mapply(function(x,y) {
    (x - as.matrix(y[,12:dim(y)[2]]) %*% t(as.matrix(y[,12:dim(y)[2]])))[rv.index, rv.index]
  }, x = cov.merge.list, y = sumstat.merge.list, SIMPLIFY = FALSE))

  rm(list=setdiff(ls(), c("info","U.merge","cov.merge")))
  gc()

  return(list(info=info,
              U=U.merge,
              cov=as.matrix(cov.merge)))
}

