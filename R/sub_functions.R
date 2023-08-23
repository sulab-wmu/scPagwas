
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
#' the code comes from https://xuzhougeng.top/archives/R-Sparse-Matrix-Note
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

##  random -- An R interface to the random.org service
##
##  Copyright (C) 2006 - 2015  Dirk Eddelbuettel <edd@debian.org>
##
##  This file is part of random
##
##  random is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 2 of the License, or
##  (at your option) any later version.
##
##  random is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.

getConnection <- function(urltxt, ...) {
  if (getRversion() >= "3.2.0" && capabilities()["libcurl"]) {
    url(urltxt, ..., method="libcurl")
  } else {
    curl(urltxt, ...)
  }
}

closeConnection <- function(con) {
  if (getRversion() >= "3.2.0" && capabilities()["libcurl"]) {
    close(con)
  } else {
    ## nothing to be done for curl()
  }
}

randomStrings <- function(n=10, len=5, digits=TRUE,
                          upperalpha=TRUE, loweralpha=TRUE,
                          unique=TRUE, check=TRUE) {
  if (n < 1 || n > 1e4)
    stop("Random string requests must be between 1 and 10,000 numbers")
  if (len < 1 || len > 20)
    stop("Random string length must be between 1 and 20")
  if (class(digits)!="logical" || class(upperalpha)!="logical" ||
      class(loweralpha)!="logical" || class(unique)!="logical")
    stop("The 'digits', '(lower|upper)alpha' and 'unique' arguments has to be logical")
  if ( !digits && !upperalpha && !loweralpha)
    stop("The 'digits', 'loweralpha' and 'loweralpha' cannot all be false at the same time")
  if (check && !quotaCheck())
    stop("random.org suggests to wait until tomorrow")
  urlbase <- "https://www.random.org/strings/"
  urltxt <- paste(urlbase,
                  "?num=", n,
                  "&len=", len,
                  "&digits=", ifelse(digits, "on", "off"),
                  "&upperalpha=", ifelse(upperalpha, "on", "off"),
                  "&loweralpha=", ifelse(loweralpha, "on", "off"),
                  "&unique=", ifelse(unique, "on", "off"),
                  "&format=plain",
                  "&rnd=new",
                  sep="")
  con <- getConnection(urltxt, open="r")
  randStrings <- as.matrix(read.table(con))
  closeConnection(con)
  return(randStrings)
}

