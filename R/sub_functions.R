
#' bh.adjust
#' @description BH P-value adjustment with a log option
#' @param x data
#' @param log default is FALSE
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



#' Convert a large sparse matrix into a dense matrix without errors
#' @description
#' Avoid the following error
#' . Error in asMethod(object) :
#' . Cholmod error 'problem too large' at file ../Core/cholmod_dense.c,
#'
#' @param mat a big sparse matrix of a type coercible to dense Matrix::Matrix
#'
#' @return a matrix
#' @export


as_matrix <- function(mat){
  Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){

  int k = z.size() ;

  IntegerMatrix  mat(nrows, ncols);

  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }

  return mat;
}
'
  )
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])

  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}


#' various association measures between sparse matrices
#'
#' @description Pearson correlation matrix between columns of X, Y
#' modification the code from:
#' http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-
#' over-a-sparse-matrix-in-r
#' https://github.com/cysouw/qlcMatrix/blob/master/R/assoc.R
#' @note Note that results larger than 1e4 x 1e4 will become very slow,
#'  because the resulting matrix is not sparse anymore.
#'
#' @param X matrix
#' @param Y matrix or vector
#'
#' @return correlation cofficient
#' @export

corSparse <- function(X, Y) {
  n <- nrow(X)
  muX <- colMeans(X)

  stopifnot(nrow(X) == nrow(Y))

  muY <- colMeans(Y)
  covmat <- (as.matrix(crossprod(X, Y)) - n * tcrossprod(muX, muY))
  sdvecX <- sqrt((colSums(X^2) - n * muX^2))
  sdvecY <- sqrt((colSums(Y^2) - n * muY^2))
  cormat <- covmat / tcrossprod(sdvecX, sdvecY)
  return(cormat)
}
