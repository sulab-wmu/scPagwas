
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
#' @export
#' @examples
rankPvalue<-function(datS, columnweights = NULL,
                     na.last = "keep",
                     ties.method = "average",
                    calculateQvalue = T,
                    pValueMethod = "rank")
{
  no.rows = dim(datS)[[1]]
  no.cols = dim(datS)[[2]]

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

    observed.sumPercentileshigh = as.numeric(rowSums(datrankshigh, na.rm = TRUE))
    Zstatistichigh = (observed.sumPercentileshigh - expectedsum)/sqrt(varsum)
    pValueLow = pnorm((Zstatisticlow))
    pValueHigh = pnorm((Zstatistichigh))
    datoutrank = data.frame(pValueLow, pValueHigh)
    if (calculateQvalue) {
      qValueLow = rep(NA, dim(datS)[[1]])
      qValueHigh = rep(NA, dim(datS)[[1]])
      #qValueExtreme = rep(NA, dim(datS)[[1]])
      rest1 = !is.na(pValueLow)
      qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
      rest1 = !is.na(pValueHigh)
      qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
      datq = data.frame(qValueLow, qValueHigh)
      datoutrank = data.frame(datoutrank, datq)
      names(datoutrank) = paste(names(datoutrank), "Rank",
                                sep = "")
    }
  }
  if (pValueMethod != "rank") {
    datSpresent = !is.na(datS) + 0
    scaled.datS = scale(datS)

    expected.value = rep(0, no.rows)
    varsum = rowSums(datSpresent^2) * 1
    observed.sumScaleddatS = as.numeric(rowSums(scaled.datS, na.rm = TRUE))
    Zstatisticlow = (observed.sumScaleddatS - expected.value)/sqrt(varsum)
    scaled.minusdatS = scale(-datS)

    observed.sumScaledminusdatS = as.numeric(rowSums(scaled.minusdatS, na.rm = TRUE))
    Zstatistichigh = (observed.sumScaledminusdatS - expected.value)/sqrt(varsum)
    pValueLow = pnorm((Zstatisticlow))
    pValueHigh = pnorm((Zstatistichigh))
    datoutscale = data.frame(pValueLow, pValueHigh)
    if (calculateQvalue) {
      qValueLow = rep(NA, dim(datS)[[1]])
      qValueHigh = rep(NA, dim(datS)[[1]])
      #qValueExtreme = rep(NA, dim(datS)[[1]])
      rest1 = !is.na(pValueLow)
      qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
      rest1 = !is.na(pValueHigh)
      qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
      datq = data.frame(qValueLow, qValueHigh)
      datoutscale = data.frame(datoutscale, datq)
    }
    names(datoutscale) = paste(names(datoutscale), "Scale",sep = "")
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


#' @title
#' Estimate the q-values for a given set of p-values
#'
#' @description
#' Estimate the q-values for a given set of p-values.  The q-value of a
#' test measures the proportion of false positives incurred (called the
#' false discovery rate) when that particular test is called significant.
#'
#' @details
#' The function \code{\link{pi0est}} is called internally and calculates the estimate of \eqn{\pi_0}{pi_0},
#' the proportion of true null hypotheses. The function \code{\link{lfdr}} is also called internally and
#' calculates the estimated local FDR values.  Arguments for these functions can be included via \code{...} and
#' will be utilized in the internal calls made in \code{\link{qvalue}}. See \url{http://genomine.org/papers/Storey_FDR_2011.pdf}
#' for a brief introduction to FDRs and q-values.
#'
#' @param p A vector of p-values (only necessary input).
#' @param fdr.level A level at which to control the FDR. Must be in (0,1]. Optional; if this is
#' selected, a vector of TRUE and FALSE is returned that specifies
#' whether each q-value is less than fdr.level or not.
#' @param pfdr An indicator of whether it is desired to make the
#' estimate more robust for small p-values and a direct finite sample estimate of pFDR -- optional.
#' @param lfdr.out If TRUE then local false discovery rates are returned. Default is TRUE.
#' @param pi0 It is recommended to not input an estimate of pi0. Experienced users can use their own methodology to estimate
#' the proportion of true nulls or set it equal to 1 for the BH procedure.
#' @param \ldots Additional arguments passed to \code{\link{pi0est}} and \code{\link{lfdr}}.
#'
#'
#' @return
#' @references
#' Storey JD. (2011) False discovery rates. In \emph{International Encyclopedia of Statistical Science}. \cr
#' \url{http://genomine.org/papers/Storey_FDR_2011.pdf} \cr
#' \url{http://www.springer.com/statistics/book/978-3-642-04897-5}
#'
#' @examples
#' # import data
#' data(hedenfalk)
#' p <- hedenfalk$p
#'
#' # get q-value object
#' qobj <- qvalue(p)
#' plot(qobj)
#' hist(qobj)
#'
#' # options available
#' qobj <- qvalue(p, lambda=0.5, pfdr=TRUE)
#' qobj <- qvalue(p, fdr.level=0.05, pi0.method="bootstrap", adj=1.2)
#'
#' @author John D. Storey
#' @keywords qvalue
#' @aliases qvalue
#' @import splines ggplot2 reshape2
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @export
qvalue <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...) {
  # Argument checks
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  }

  # Calculate pi0 estimate
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1)  {
      pi0s = list()
      pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }

  # Calculate q-value estimates
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 - (1 - p[o]) ^ m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  }
  qvals_out[rm_na] <- qvals
  # Calculate local FDR estimates
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }

  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda,
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}


#' @title Proportion of true null p-values
#' @description
#' Estimates the proportion of true null p-values, i.e., those following the Uniform(0,1) distribution.
#'
#' @param p A vector of p-values (only necessary input).
#' @param lambda The value of the tuning parameter to estimate
#' \eqn{\pi_0}{pi_0}. Must be in [0,1). Optional, see Storey (2002).
#' @param pi0.method Either "smoother" or "bootstrap"; the method for
#' automatically choosing tuning parameter in the estimation of \eqn{\pi_0}{pi_0},
#' the proportion of true null hypotheses.
#' @param smooth.df Number of degrees-of-freedom to use when estimating \eqn{\pi_0}{pi_0}
#' with a smoother. Optional.
#' @param smooth.log.pi0 If TRUE and \code{pi0.method} = "smoother", \eqn{\pi_0}{pi_0} will be
#' estimated by applying a smoother to a scatterplot of \eqn{\log(\pi_0)}{log(pi_0)} estimates
#' against the tuning parameter \eqn{\lambda}{lambda}. Optional.
#' @param \dots Arguments passed from \code{\link{qvalue}} function.
#'
#' @details
#' If no options are selected, then the method used to estimate \eqn{\pi_0}{pi_0} is
#' the smoother method described in Storey and Tibshirani (2003). The
#' bootstrap method is described in Storey, Taylor & Siegmund (2004). A closed form solution of the
#' bootstrap method is used in the package and is significantly faster.
#'
#' @return Returns a list:
#' \item{pi0}{A numeric that is the estimated proportion
#' of true null p-values.}
#' \item{pi0.lambda}{A vector of the proportion of null
#' values at the \eqn{\lambda}{lambda} values (see vignette).}
#' \item{lambda}{A vector of \eqn{\lambda}{lambda} value(s) utilized in calculating \code{pi0.lambda}.}
#' \item{pi0.smooth}{A vector of fitted values from the
#' smoother fit to the \eqn{\pi_0}{pi_0} estimates at each \code{lambda} value
#' (pi0.method="bootstrap" returns NULL).}
#'
#' @references
#' Storey JD. (2002) A direct approach to false discovery rates. Journal
#' of the Royal Statistical Society, Series B, 64: 479-498. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00346/abstract}
#'
#' Storey JD and Tibshirani R. (2003) Statistical significance for
#' genome-wide experiments. Proceedings of the National Academy of Sciences,
#' 100: 9440-9445. \cr
#" \url{http://www.pnas.org/content/100/16/9440.full}
#'
#' Storey JD. (2003) The positive false discovery rate: A Bayesian
#' interpretation and the q-value. Annals of Statistics, 31: 2013-2035. \cr
#' \url{http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdf_1&handle=euclid.aos/1074290335}
#'
#' Storey JD, Taylor JE, and Siegmund D. (2004) Strong control,
#' conservative point estimation, and simultaneous conservative
#' consistency of false discovery rates: A unified approach. Journal of
#' the Royal Statistical Society, Series B, 66: 187-205. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2004.00439.x/abstract}
#'
#' Storey JD. (2011) False discovery rates. In \emph{International Encyclopedia of Statistical Science}. \cr
#' \url{http://genomine.org/papers/Storey_FDR_2011.pdf} \cr
#' \url{http://www.springer.com/statistics/book/978-3-642-04897-5}
#'
#' @examples
#' # import data
#' data(hedenfalk)
#' p <- hedenfalk$p
#'
#' # proportion of null p-values
#' nullRatio <- pi0est(p)
#' nullRatioS <- pi0est(p, lambda=seq(0.40, 0.95, 0.05), smooth.log.pi0="TRUE")
#' nullRatioM <- pi0est(p, pi0.method="bootstrap")
#'
#' # check behavior of estimate over lambda
#' # also, pi0est arguments can be passed to qvalue
#' qobj = qvalue(p, lambda=seq(0.05, 0.95, 0.1), smooth.log.pi0="TRUE")
#' hist(qobj)
#' plot(qobj)
#'
#' @author John D. Storey
#' @seealso \code{\link{qvalue}}
#' @keywords pi0est, proportion true nulls
#' @aliases pi0est
#' @export
pi0est <- function(p, lambda = seq(0.05,0.95,0.05), pi0.method = c("smoother", "bootstrap"),
                   smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
  # Check input arguments
  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method = match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda) # guard against user input

  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  } else if (ll > 1 && ll < 4) {
    stop(sprintf(paste("ERROR:", paste("length(lambda)=", ll, ".", sep=""),
                       "If length of lambda greater than 1,",
                       "you need at least 4 values.")))
  } else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  # Determines pi0
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  } else {
    ind <- length(lambda):1
    pi0 <- cumsum(tabulate(findInterval(p, vec=lambda))[ind]) / (length(p) * (1-lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    # Smoother method approximation
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      } else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    } else if (pi0.method == "bootstrap") {
      # Bootstrap method closed form solution by David Robinson
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W / (m ^ 2 * (1 - lambda) ^ 2)) * (1 - W / m) + (pi0 - minpi0) ^ 2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    } else {
      stop('ERROR: pi0.method must be one of "smoother" or "bootstrap".')
    }
  }
  if (pi0 <= 0) {
    stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda,
              lambda = lambda, pi0.smooth = pi0Smooth))
}

#' @title Estimate local False Discovery Rate (FDR)
#'
#' @description
#' Estimate the local FDR values from p-values.
#'
#' @param p A vector of p-values (only necessary input).
#' @param pi0 Estimated proportion of true null p-values. If NULL, then \code{\link{pi0est}} is called.
#' @param trunc If TRUE, local FDR values >1 are set to 1. Default is TRUE.
#' @param monotone If TRUE, local FDR values are non-decreasing with increasing p-values. Default is TRUE; this is recommended.
#' @param transf Either a "probit" or "logit" transformation is applied to the p-values so that a local FDR estimate can be formed that
#' does not involve edge effects of the [0,1] interval in which the p-values lie.
#' @param adj Numeric value that is applied as a multiple of the smoothing bandwidth used in the density estimation. Default is \code{adj=1.0}.
#' @param eps Numeric value that is threshold for the tails of the empirical p-value distribution. Default is 10^-8.
#' @param \ldots Additional arguments, passed to \code{\link{pi0est}}.
#'
#' @details It is assumed that null p-values follow a Uniform(0,1) distribution.
#'   The estimated proportion of true null hypotheses \eqn{\hat{\pi}_0}{pi_0} is either
#'   a user-provided value or the value calculated via \code{\link{pi0est}}.
#'   This function works by forming an estimate of the marginal density of the
#'   observed p-values, say \eqn{\hat{f}(p)}{f(p)}. Then the local FDR is estimated as
#'   \eqn{{\rm lFDR}(p) = \hat{\pi}_0/\hat{f}(p)}{lFDR(p) = pi_0/f(p)}, with
#'   adjustments for monotonicity and to guarantee that \eqn{{\rm lFDR}(p) \leq
#'   1}{lFDR(p) <= 1}.  See the Storey (2011) reference below for a concise
#'   mathematical definition of local FDR.
#'
#' @return
#' A vector of estimated local FDR values, with each entry corresponding to the entries of the input p-value vector \code{p}.
#'
#' @references
#' Efron B, Tibshirani R, Storey JD, and Tisher V. (2001) Empirical Bayes analysis
#' of a microarray experiment. Journal of the American Statistical Association, 96: 1151-1160. \cr
#' \url{http://www.tandfonline.com/doi/abs/10.1198/016214501753382129}
#'
#' Storey JD. (2003) The positive false discovery rate: A Bayesian
#' interpretation and the q-value. Annals of Statistics, 31: 2013-2035. \cr
#' \url{http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdf_1&handle=euclid.aos/1074290335}
#'
#' Storey JD. (2011) False discovery rates. In \emph{International Encyclopedia of Statistical Science}. \cr
#' \url{http://genomine.org/papers/Storey_FDR_2011.pdf} \cr
#' \url{http://www.springer.com/statistics/book/978-3-642-04897-5}
#'
#' @examples
#' # import data
#' data(hedenfalk)
#' p <- hedenfalk$p
#' lfdrVals <- lfdr(p)
#'
#' # plot local FDR values
#' qobj = qvalue(p)
#' hist(qobj)
#'
#' @author John D. Storey
#' @seealso \code{\link{qvalue}}, \code{\link{pi0est}}, \code{\link{hist.qvalue}}
#' @aliases lfdr
#' @keywords local False Discovery Rate, lfdr
#' @export
lfdr <- function(p, pi0 = NULL, trunc = TRUE, monotone = TRUE,
                 transf = c("probit", "logit"), adj = 1.5, eps = 10 ^ -8, ...) {
  # Check inputs
  lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range [0,1].")
  } else if (is.null(pi0)) {
    pi0 <- pi0est(p, ...)$pi0
  }
  n <- length(p)
  transf <- match.arg(transf)
  # Local FDR method for both probit and logit transformations
  if (transf == "probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)
    x <- qnorm(p)
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0 * dnorm(x) / y
  } else {
    x <- log((p + eps) / (1 - p + eps))
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1 + exp(x)) ^ 2
    lfdr <- (pi0 * dx) / y
  }
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {
    o <- order(p, decreasing = FALSE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  return(lfdr_out)
}

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
    IRanges::IRanges(
      as.numeric(snp[, "pos"]),
      as.numeric(snp[, "pos"])
    ),
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
  idx <- setdiff(seq_len(nrow(d_both)), idx)
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
as_matrix <- function(mat) {
  tmp <- matrix(data = 0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i + 1
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1]) + 1
  val <- mat@x
  for (i in seq_along(val)) {
    tmp[row_pos[i], col_pos[i]] <- val[i]
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
  # X <-as_matrix(X)
  n <- nrow(X)
  muX <- colMeans(X)

  stopifnot(nrow(X) == nrow(Y))

  muY <- colMeans(Y)
  covmat <- (as.matrix(crossprod(X, Y)) - n * tcrossprod(muX, muY))
  sdvecX <- sqrt((colSums(X^2) - n * muX^2))
  sdvecY <- sqrt((colSums(Y^2) - n * muY^2))
  cormat <- covmat / tcrossprod(sdvecX, sdvecY)
  # cormat[is.nan(cormat),1]<-0
  return(cormat)
}
