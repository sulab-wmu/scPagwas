% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sub_functions.R
\name{as_matrix}
\alias{as_matrix}
\title{Convert a large sparse matrix into a dense matrix without errors}
\usage{
as_matrix(mat)
}
\arguments{
\item{mat}{a big sparse matrix of a type coercible to dense Matrix::Matrix}
}
\value{
a matrix
}
\description{
Avoid the following error
. Error in asMethod(object) :
. Cholmod error 'problem too large' at file ../Core/cholmod_dense.c,
the code comes from https://xuzhougeng.top/archives/R-Sparse-Matrix-Note
}
