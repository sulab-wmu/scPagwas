
#' link_scCell_pwpca_block
#' @description link the single cell pca score and expression for each
#' pathway genes for each block(single cell)
#' construct the final matrix for regression in blocks
#' @param Pagwas Pagwas format, deault is NULL.
#' @param n.cores Parallel cores,default is 1. use detectCores() to
#' check the cores in computer.
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' # Pagwas should have inhibit data
#' Pagwas <- link_pwpca_block(Pagwas)
link_scCell_pwpca_block <- function(Pagwas, n.cores = 1) {
  scexpr <- Seurat::GetAssayData(object = Single_data, slot = "data")
  SOAR::Store(Single_data)

  if (is.null(Pathway_ld_gwas_data)) {
    message("* no loaded Pathway_ld_gwas data")
    return(Pagwas)
  }

  if (dim(scexpr)[2] == dim(pca_scCell_mat)[2]) {
    colnames(scexpr) <- colnames(pca_scCell_mat)
  } else {
    a <- intersect(colnames(scexpr), colnames(pca_scCell_mat))
    if (length(a) == 0) {
      stop("* please check the colnames of Singlecell data.There should have no specific symbol")
    }
    scexpr <- scexpr[, a]
    pca_scCell_mat <- pca_scCell_mat[, a]
  }
  #rm(a)


  message("* Merging pathway score and expression information about blocks in ", length(Pathway_ld_gwas_data), " pathways")
  pb <- txtProgressBar(style = 3)

  Pathway_ld_gwas_data <- papply(Pathway_ld_gwas_data, function(pa_block) {
    pathway <- unique(pa_block$block_info$pathway)

    x <- pca_scCell_mat[pathway, ]

    if(length(pathway)==1){
      x <- matrix(x, nrow = 1)
      rownames(x)<-pathway
    }

    if (nrow(pa_block$snps) == 0) {
      pa_block$include_in_inference <- F
      pa_block$x <- NULL # to make sure we totally replace previous stuffs
      return(pa_block)

    }

    mg <- intersect(Pagwas$rawPathway_list[[pathway]], rownames(scexpr))
    if (length(mg) < 2) {
      return(NULL)
    }

    tryCatch(
      {
        x_FBM <- bigstatsr::as_FBM(as.matrix(scexpr[mg, ]))
      }, error = function(e) {
        x_FBM <- bigstatsr::as_FBM(as_matrix(scexpr[mg, ]))
      })

    if (length(mg) > 1) {

      x2 <- bigstatsr::big_apply(x_FBM, a.FUN = colnorm_sub, a.combine = "cbind")
      rownames(x2) <- mg
      colnames(x2) <- colnames(scexpr)
    }

    if (pa_block$n_snps > 1) {
      x2 <-x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- x[rep(1, pa_block$n_snps), ]
      rownames(x) <- pa_block$snps$rsid
      #snp_gene_df <- snp_gene_df
      rownames(snp_gene_df) <- snp_gene_df$rsid


      x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
      x3 <- x2 * x
      rm(x)
      rm(x2)
    } else {
      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label

      pa_block$n_snps <- nrow(pa_block$snps)

      x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
      rownames(x) <- pa_block$snps$rsid

      #snp_gene_df <- snp_gene_df
      rownames(snp_gene_df) <- snp_gene_df$rsid

      x <- matrix(as.numeric(x) * as.numeric(snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
      x3 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
      rm(x)
      rm(x2)
    }

    ## add the slope of eqtls for x
    rownames(x3) <- pa_block$snps$rsid

    pa_block$x <- Matrix::crossprod(t(pa_block$ld_matrix_squared), x3)
    rm(x3)
    pa_block$include_in_inference <- T

    setTxtProgressBar(pb, which(names(Pathway_ld_gwas_data) == pathway) / length(Pathway_ld_gwas_data))

    return(pa_block)
  }, n.cores = n.cores)
  close(pb)
  rm(scexpr)

  Pathway_ld_gwas_data <- Pathway_ld_gwas_data[!sapply(Pathway_ld_gwas_data, is.null)]
  SOAR::Store(Pathway_ld_gwas_data)
  SOAR::Store(pca_scCell_mat)
  SOAR::Store(snp_gene_df)

  return(Pagwas)
}



#' colnorm_sub
#'
#' @param X  data
#' @param ind index
#'
#' @return
#'
colnorm_sub <- function(X, ind) {
  if (sum(X[, ind]) == 0) {
    return(y)
  } else {
    return(X[, ind] / sum(X[, ind]))
  }
}

#' link_scCell_pwpca_block
#' @description link the single cell pca score and expression for each
#' pathway genes for each block(single cell)
#' construct the final matrix for regression in blocks
#' @param Pagwas Pagwas format, deault is NULL.
#' @param n.cores Parallel cores,default is 1. use detectCores() to
#' check the cores in computer.
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' # Pagwas should have inhibit data
#' Pagwas <- link_pwpca_block(Pagwas)
link_scCell_pwpca_block <- function(Pagwas, n.cores = 1) {
  scexpr <- Seurat::GetAssayData(object = Single_data, slot = "data")
  SOAR::Store(Single_data)


  if (is.null(Pathway_ld_gwas_data)) {
    message("* no loaded Pathway_ld_gwas data")
    return(Pagwas)
  }

  if (dim(scexpr)[2] == dim(pca_scCell_mat)[2]) {
    colnames(scexpr) <- colnames(pca_scCell_mat)
  } else {
    a <- intersect(colnames(scexpr), colnames(pca_scCell_mat))
    if (length(a) == 0) {
      stop("* please check the colnames of Singlecell data.There should have no specific symbol")
    }
    scexpr <- scexpr[, a]
    pca_scCell_mat <- pca_scCell_mat[, a]
  }

  message("* Merging pathway score and expression information about blocks in ", length(Pathway_ld_gwas_data), " pathways")
  pb <- txtProgressBar(style = 3)

  Pathway_ld_gwas_data <- papply(Pathway_ld_gwas_data, function(pa_block) {
    pathway <- unique(pa_block$block_info$pathway)

    x <- pca_scCell_mat[pathway, ]

    if(length(pathway)==1){
      x <- matrix(x, nrow = 1)
      rownames(x)<-pathway
    }

    if (nrow(pa_block$snps) == 0) {
      pa_block$include_in_inference <- F
      pa_block$x <- NULL # to make sure we totally replace previous stuffs
      return(pa_block)

    }

    mg <- intersect(Pagwas$rawPathway_list[[pathway]], rownames(scexpr))
    if (length(mg) < 2) {
      return(NULL)
    }

    tryCatch(
      {
        x_FBM <- bigstatsr::as_FBM(as.matrix(scexpr[mg, ]))
      }, error = function(e) {
        x_FBM <- bigstatsr::as_FBM(as_matrix(scexpr[mg, ]))
      })

    if (length(mg) > 1) {

      x2 <- bigstatsr::big_apply(x_FBM, a.FUN = colnorm_sub, a.combine = "cbind")
      rownames(x2) <- mg
      colnames(x2) <- colnames(scexpr)
    }

    if (pa_block$n_snps > 1) {
      x2 <-x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- x[rep(1, pa_block$n_snps), ]

      rownames(x) <- pa_block$snps$rsid
      #snp_gene_df <- snp_gene_df
      rownames(snp_gene_df) <- snp_gene_df$rsid

      x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
      x3 <- x2 * x
      rm(x)
      rm(x2)
    } else {
      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label

      pa_block$n_snps <- nrow(pa_block$snps)

      x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
      rownames(x) <- pa_block$snps$rsid

      #snp_gene_df <- snp_gene_df
      rownames(snp_gene_df) <- snp_gene_df$rsid

      x <- matrix(as.numeric(x) * as.numeric(snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
      x3 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
      rm(x)
      rm(x2)
    }

    ## add the slope of eqtls for x
    pa_block$x <- Matrix::crossprod(t(pa_block$ld_matrix_squared), x3)
    colnames(pa_block$x)<-colnames(scexpr)
    rownames(pa_block$x)<-pa_block$snps$rsid
    rm(x3)
    pa_block$include_in_inference <- T

    setTxtProgressBar(pb, which(names(Pathway_ld_gwas_data) == pathway) / length(Pathway_ld_gwas_data))

    return(pa_block)
  }, n.cores = n.cores)
  close(pb)
  rm(scexpr)

  Pathway_ld_gwas_data <- Pathway_ld_gwas_data[!sapply(Pathway_ld_gwas_data, is.null)]
  SOAR::Store(Pathway_ld_gwas_data)
  SOAR::Store(pca_scCell_mat)
  SOAR::Store(snp_gene_df)

  return(Pagwas)
}



#' colnorm_sub
#'
#' @param X  data
#' @param ind index
#'
#' @return
#'
colnorm_sub <- function(X, ind) {
  if (sum(X[, ind]) == 0) {
    return(y)
  } else {
    return(X[, ind] / sum(X[, ind]))
  }
}
