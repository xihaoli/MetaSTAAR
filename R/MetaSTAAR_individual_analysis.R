#' Meta-analysis of individual variants using \code{MetaSTAAR}
#'
#' The \code{MetaSTAAR_individual_analysis} function takes in the summary statistics file
#' from each individual study (the output from \code{\link{MetaSTAAR_worker_sumstat}})
#' to analyze the associations between a quantitative/dichotomous phenotype and
#' all individual variants in the merged variant list by using the meta-analysis
#' @param chr a numeric value indicating the chromosome of the genetic region of interest.
#' @param start.loc a numeric value indicating the starting location (position) of the
#' genetic region of interest.
#' @param end.loc a numeric value indicating the ending location (position) of the
#' genetic region of interest.
#' @param study.names a character vector containing the name of each individual
#' study contributing to the meta-analysis.
#' @param sample.sizes a numeric vector with the length of \code{study.names}
#' indicating the sample size of each study.
#' @param sumstat.dir a character vector containing the directories of the study-specific summary statistics file folders.
#' @param common_mac_cutoff the cutoff of minimum combined minor allele count (inclusive) in
#' defining "common" variants.
#' @param trait a character value indicating the underlying trait of interest for
#' the meta-analysis.
#' @param segment.size a numeric value indicating the length of each segment of which
#' the summary statistics and sparse weighted covariance files are stored.
#' Note that the input value should be aligned with the inpue values of
#' \code{\link{MetaSTAAR_worker_cov}} (default = 5e+05).
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}
#' specified in \code{\link{MetaSTAAR_worker_sumstat}} and \code{\link{MetaSTAAR_worker_cov}}.
#' If \code{check_qc_label} is FALSE, it is assumed that no variant will be dropped (Default = FALSE).
#' @return a data frame with p rows corresponding to the p genetic variants in the merged variant list
#' and the following columns: chromosome (chr), position (pos), reference allele (ref), alternative allele (alt),
#' combined alternative allele count (alt_AC), combined minor allele count (MAC), combined minor allele frequency (MAF),
#' combined sample size (N), the score test p-value (p), the log score test p-value (logp), the score test statistic (Score),
#' the standard error associated with the score test statistic (Score_SE), the estimated effect size of the minor allele (Est),
#' the standard error associated with the estimated effect size (Est_se).
#' If a variant in the merged variant list has standard error equal to 0, the p-value will be set as 1.

MetaSTAAR_individual_analysis <- function(chr,start.loc,end.loc,study.names,sample.sizes,sumstat.dir,
                                          common_mac_cutoff,trait,segment.size = 5e5,check_qc_label=FALSE){

  cov_maf_cutoff <- rep(0.5 + 1e-16,length(study.names))
  segment <- floor((start.loc - 1) / segment.size) + 1
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
  }

  ### select "common" variant based on the input cutoff
  alt_AC.merge <- as.integer(Reduce("+",lapply(sumstat.merge.list, function(x) {x$alt_AC})))
  N.merge.nonzero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC != 0)}))
  alt_AF.merge.nonzero <- alt_AC.merge / (2 * N.merge.nonzero)
  N.merge.zero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC == 0)}))
  alt_AC.merge <- alt_AC.merge + (alt_AF.merge.nonzero > 0.5) * (2 * N.merge.zero)
  N.merge <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N}))
  MAC.merge <- pmin(alt_AC.merge, 2 * N.merge - alt_AC.merge)
  MAF.merge <- MAC.merge / (2 * N.merge)
  cv.index <- (MAF.merge>(common_mac_cutoff - 0.5) / (2 * N.merge)) & Reduce("*",mapply(function(x,y) {
    x$MAF<y
  }, x = sumstat.merge.list, y = cov_maf_cutoff, SIMPLIFY = FALSE))
  if (check_qc_label){
    cv.index <- cv.index & Reduce("*",lapply(sumstat.merge.list, function(x) {
      x$qc_label=="PASS"
    }))
  }

  info <- cbind(sumstat.varid.nodup[,c("chr","pos","ref","alt")],
                alt_AC=alt_AC.merge,MAC=MAC.merge,MAF=MAF.merge,N=N.merge)[cv.index,]

  U.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
    x <- x[cv.index,]
    return(x$U)
  }))

  V.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
    x <- x[cv.index,]
    return(x$V)
  }))

  p.merge <- pchisq(U.merge^2/V.merge,df=1,lower.tail=FALSE)
  logp.merge <- -pchisq(U.merge^2/V.merge,df=1,lower.tail=FALSE,log.p=TRUE)

  rm(list=setdiff(ls(), c("info","U.merge","V.merge","p.merge","logp.merge")))
  gc()

  return(data.frame(chr=info$chr,
                    pos=info$pos,
                    ref=info$ref,
                    alt=info$alt,
                    alt_AC=info$alt_AC,
                    MAC=info$MAC,
                    MAF=info$MAF,
                    N=info$N,
                    p=p.merge,
                    logp=logp.merge,
                    Score=U.merge,
                    Score_se=sqrt(V.merge),
                    Est=U.merge/V.merge,
                    Est_se=1/sqrt(V.merge)))
}
