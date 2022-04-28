
#' rankPvalue
#'
#' @param datS matrix or data.frame
#' @param columnweights NULL
#' @param na.last "keep"
#' @param ties.method "average"
#' @param calculateQvalue TRUE
#' @param pValueMethod "all"
#'
#' @return
#'
#' @examples
rankPvalue<-function(datS, columnweights = NULL,
                     na.last = "keep",
                     ties.method = "average",
                    calculateQvalue = TRUE,
                    pValueMethod = "all")
{
  no.rows = dim(datS)[[1]]
  no.cols = dim(datS)[[2]]
  if (!is.null(columnweights) & no.cols != length(columnweights))
    stop("The number of components of the vector columnweights is unequal to the number of columns of datS. Hint: consider transposing datS. ")

  if (!is.null(columnweights) ) {
    if ( min(columnweights,na.rm=TRUE)<0 )  stop("At least one component of columnweights is negative, which makes no sense. The entries should be positive numbers")
    if ( sum(is.na(columnweights))>0 )  stop("At least one component of columnweights is missing, which makes no sense. The entries should be positive numbers")
    if ( sum( columnweights)!= 1 ) {
      # warning("The entries of columnweights do not sum to 1. Therefore, they will divided by the sum. Then the resulting weights sum to 1.");
      columnweights= columnweights/sum( columnweights)
    }
  }

  if (pValueMethod != "scale") {
    percentilerank1 = function(x) {
      R1 = rank(x, ties.method = ties.method, na.last = na.last)
      (R1-.5)/max(R1, na.rm = TRUE)
    }

    datrankslow = apply(datS, 2, percentilerank1)
    if (!is.null(columnweights)) {
      datrankslow = t(t(datrankslow) * columnweights)
    }
    datSpresent = !is.na(datS) + 0
    if (!is.null(columnweights)) {
      datSpresent = t(t(datSpresent) * columnweights)
    }
    expectedsum = rowSums(datSpresent, na.rm = TRUE) *
      0.5
    varsum = rowSums(datSpresent^2, na.rm = TRUE) * 1/12
    observed.sumPercentileslow = as.numeric(rowSums(datrankslow, na.rm = TRUE))
    Zstatisticlow = (observed.sumPercentileslow - expectedsum)/sqrt(varsum)
    datrankshigh = apply(-datS, 2, percentilerank1)
    if (!is.null(columnweights)) {
      datrankshigh = t(t(datrankshigh) * columnweights)
    }
    observed.sumPercentileshigh = as.numeric(rowSums(datrankshigh, na.rm = TRUE))
    Zstatistichigh = (observed.sumPercentileshigh - expectedsum)/sqrt(varsum)
    pValueLow = pnorm((Zstatisticlow))
    pValueHigh = pnorm((Zstatistichigh))
    pValueExtreme = pmin(pValueLow, pValueHigh)
    datoutrank = data.frame(pValueExtreme, pValueLow, pValueHigh)
    if (calculateQvalue) {
      qValueLow = rep(NA, dim(datS)[[1]])
      qValueHigh = rep(NA, dim(datS)[[1]])
      qValueExtreme = rep(NA, dim(datS)[[1]])
      rest1 = !is.na(pValueLow)
      qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
      rest1 = !is.na(pValueHigh)
      qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
      rest1 = !is.na(pValueExtreme)
      qValueExtreme = pmin(qValueLow, qValueHigh)
      datq = data.frame(qValueExtreme, qValueLow, qValueHigh)
      datoutrank = data.frame(datoutrank, datq)
      names(datoutrank) = paste(names(datoutrank), "Rank",
                                sep = "")
    }
  }
  if (pValueMethod != "rank") {
    datSpresent = !is.na(datS) + 0
    scaled.datS = scale(datS)
    if (!is.null(columnweights)) {
      scaled.datS = t(t(scaled.datS) * columnweights)
      datSpresent = t(t(datSpresent) * columnweights)
    }
    expected.value = rep(0, no.rows)
    varsum = rowSums(datSpresent^2) * 1
    observed.sumScaleddatS = as.numeric(rowSums(scaled.datS, na.rm = TRUE))
    Zstatisticlow = (observed.sumScaleddatS - expected.value)/sqrt(varsum)
    scaled.minusdatS = scale(-datS)
    if (!is.null(columnweights)) {
      scaled.minusdatS = t(t(scaled.minusdatS) * columnweights)
    }
    observed.sumScaledminusdatS = as.numeric(rowSums(scaled.minusdatS, na.rm = TRUE))
    Zstatistichigh = (observed.sumScaledminusdatS - expected.value)/sqrt(varsum)
    pValueLow = pnorm((Zstatisticlow))
    pValueHigh = pnorm((Zstatistichigh))
    pValueExtreme = 2 * pnorm(-abs(Zstatisticlow))
    datoutscale = data.frame(pValueExtreme, pValueLow, pValueHigh)
    if (calculateQvalue) {
      qValueLow = rep(NA, dim(datS)[[1]])
      qValueHigh = rep(NA, dim(datS)[[1]])
      qValueExtreme = rep(NA, dim(datS)[[1]])
      rest1 = !is.na(pValueLow)
      qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
      rest1 = !is.na(pValueHigh)
      qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
      rest1 = !is.na(pValueExtreme)
      qValueExtreme[rest1] = qvalue(pValueExtreme[rest1])$qvalues
      datq = data.frame(qValueExtreme, qValueLow, qValueHigh)
      datoutscale = data.frame(datoutscale, datq)
    }
    names(datoutscale) = paste(names(datoutscale), "Scale",
                               sep = "")
  }
  if (pValueMethod == "rank") {
    datout = datoutrank
  }
  if (pValueMethod == "scale") {
    datout = datoutscale
  }
  if (pValueMethod != "rank" & pValueMethod != "scale")
    datout = data.frame(datoutrank, datoutscale)
  datout
} # End of function



#' papply
#' @description wrapper around different mclapply mechanisms
#' Abstracts out mclapply implementation, and defaults to lapply when only one core is requested (helps with debugging)
#' @param ... parameters to pass to lapply, mclapply, bplapply, etc.
#' @param n.cores default is 1.
#'
#' @return
#'
papply <- function(..., n.cores = detectCores()) {
  if (n.cores > 1) {
    parallel::mclapply(..., mc.cores = n.cores)
  } else {
    lapply(...)
  }
}



#' bh.adjust
#' @description BH P-value adjustment with a log option
#' @param x data
#' @param log default is FALSE
#'
#' @return
#'
bh.adjust <- function(x, log = FALSE) {
  nai <- which(!is.na(x))
  ox <- x
  x <- x[nai]
  id <- order(x, decreasing = FALSE)
  if (log) {
    q <- x[id] + log(length(x) / seq_along(x))
  } else {
    q <- x[id] * length(x) / seq_along(x)
  }
  a <- rev(cummin(rev(q)))[order(id)]
  ox[nai] <- a
  ox
}

#' Snps to Genes
#' @description Maps SNPs to their nearest genes within TSS windows.
#'
#' @note updated GenomicRanges to version 1.30.1 on 01/30/2018
#' Requires ignore.strand=TRUE param to properly run distanceToNearest()
#' given unknown strand assignment for SNP array-based genotyping info
#' @param snp (char) PLINK bim file (only CHR, SNP, and BP columns considered).
#' @param refGene (char) refseq table with header.
#' @param marg (integer) region upstream and downstream(default=0).
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(gtf_df)
#' snp_gene_df <- Snp2Gene(snp = as.data.frame(Pagwas$gwas_data), refGene = gtf_df, marg = 10000)
Snp2Gene <- function(snp, refGene, marg = 10000) {
  snp_GR <- GenomicRanges::GRanges(snp[, "chrom"],
                                   IRanges::IRanges(as.numeric(snp[, "pos"]),
                                                    as.numeric(snp[, "pos"])),
    name = snp[, "rsid"]
  )

  gene_GR <- GenomicRanges::GRanges(refGene[, "chrom"],
                                    IRanges::IRanges(
      refGene[, "start"] + 1 - marg,
      refGene[, "start"] + marg
    ),
    name = refGene[, "label"]
  )

  start_GR <- GenomicRanges::resize(gene_GR, fix = "start", width = 1L)
  end_GR <- GenomicRanges::resize(gene_GR, fix = "end", width = 1L)

  cat("* Computing distances between SNPs and genes\n")
  # SNPs inside the gene domain
  d0 <- GenomicRanges::distanceToNearest(snp_GR, gene_GR, ignore.strand = TRUE)
  dbit <- d0@elementMetadata$distance
  d_in <- data.frame(
    queryHits = d0@from,
    subjectHits = d0@to,
    distance = dbit
  )

  idx <- which(dbit == 0)

  out1 <- cbind(
    snp_GR$name[d_in$queryHits[idx]],
    gene_GR$name[d_in$subjectHits[idx]],
    as.numeric(snp[d_in$queryHits[idx], "pos"]),
    d_in$distance[idx]
  )

  # SNPs outside gene domain
  idx <- which(dbit > 0) # not in a gene
  snp2_GR <- snp_GR[idx]
  cat(sprintf("%i SNPs not inside gene domain\n", length(idx)))
  #  rm(snp_GR, snps, idx)

  cat("** Computing distance to domain starts\n")
  d1 <- GenomicRanges::distanceToNearest(snp2_GR, start_GR, ignore.strand = TRUE)
  dbit <- d1@elementMetadata$distance
  d_start <- cbind(d1@from, d1@to, dbit)
  colnames(d_start)[1:2] <- c("queryHits", "subjectHits")

  cat("** Computing distance to domain ends\n")
  d2 <- GenomicRanges::distanceToNearest(snp2_GR, end_GR, ignore.strand = TRUE)
  dbit <- d2@elementMetadata$distance
  d_end <- cbind(d2@from, d2@to, dbit)
  colnames(d_end)[1:2] <- c("queryHits", "subjectHits")

  if (all.equal(d_start[, 1], d_end[, 1]) != TRUE) {
    cat("d_start and d_end have different indexing. You need some other ")
    cat("way to include all SNPs\n")
  }

  d_both <- merge(d_start, d_end, by.x = "queryHits", by.y = "queryHits")

  # Start is closer than end
  idx <- which(d_both$dbit.x < d_both$dbit.y)

  out2 <- cbind(
    snp2_GR$name[d_both$queryHits[idx]],
    start_GR$name[d_both$subjectHits.x[idx]],
    as.numeric(snp[d_both$queryHits[idx], "pos"]),
    d_both$dbit.x[idx]
  )
  idx <- setdiff(seq_len(nrow(d_both)) , idx)
  out3 <- cbind(
    snp2_GR$name[d_both$queryHits[idx]],
    end_GR$name[d_both$subjectHits.y[idx]],
    as.numeric(snp[d_both$queryHits[idx], "pos"]),
    d_both$dbit.y[idx]
  )
  out <- as.data.frame(rbind(out1, out2, out3))
  colnames(out) <- c("rsid", "label", "pos", "Disstance")

  # out<-out[out$Disstance==0,]
  return(out)
}


#' Convert a large sparse matrix into a dense matrix without errors
#' @description
#' Avoid the following error
#' . Error in asMethod(object) :
#' . Cholmod error 'problem too large' at file ../Core/cholmod_dense.c,
#'
#' @param mat a big sparse matrix of a type coercible to dense Matrix::Matrix
#'
#' @return a matrix
#'
as_matrix <- function(mat){
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}


#' various association measures between sparse matrices
#'
#' @description Pearson correlation matrix between columns of X, Y
#' modification the code from:
#' http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
#' https://github.com/cysouw/qlcMatrix/blob/master/R/assoc.R
#' @note Note that results larger than 1e4 x 1e4 will become very slow, because the resulting matrix is not sparse anymore.
#'
#' @param X matrix
#' @param Y matrix or vector
#'
#' @return
#'

corSparse <- function(X, Y) {
  #X <-as_matrix(X)
  n <- nrow(X)
  muX <- colMeans(X)

  stopifnot( nrow(X) == nrow(Y) )

  muY <- colMeans(Y)
  covmat <- (as.matrix(crossprod(X,Y)) - n*tcrossprod(muX,muY) ) / (n-1)
  sdvecX <- sqrt( (colSums(X^2) - n*muX^2) / (n-1) )
  sdvecY <- sqrt( (colSums(Y^2) - n*muY^2) / (n-1) )
  cormat <- covmat/tcrossprod(sdvecX,sdvecY)
  #cormat[is.nan(cormat),1]<-0
  return(cormat)
}




