
#' link_scCell_pwpca_block
#' @description link the single cell pca score and expression for each
#' pathway genes for each block(single cell)
#' construct the final matrix for regression in blocks
#' @param Pagwas Pagwas format, deault is NULL.
#' @param n.cores Parallel cores,default is 1. use detectCores() to
#' check the cores in computer.
#' @param scPagwasSession "scPagwasSession"
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)


link_scCell_pwpca_block <- function(Pagwas, n.cores = 1,scPagwasSession="scPagwasSession") {

  Sys.setenv(R_LOCAL_CACHE=scPagwasSession)

  if (is.null(Pathway_ld_gwas_data)) {
    message("* no loaded Pathway_ld_gwas data")
    return(Pagwas)
  }

  if (dim(data_mat)[2] == dim(pca_scCell_mat)[2]) {
    colnames(data_mat) <- colnames(pca_scCell_mat)
    data_mat <- data_mat
    pca_scCell_mat <- pca_scCell_mat
  } else {
    a <- intersect(colnames(data_mat), colnames(pca_scCell_mat))
    if (length(a) == 0) {
      stop("* please check the colnames of Singlecell data.There should have no specific symbol")
    }

    data_mat <- data_mat[, a]
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

    mg <- intersect(Pagwas$rawPathway_list[[pathway]], rownames(data_mat))
    if (length(mg) < 2) {
      return(NULL)
    }

    if (length(mg) > 1) {

      x2 <- bigstatsr::big_apply(bigstatsr::as_FBM(data_mat[mg, ]), a.FUN = colnorm_sub, a.combine = "cbind")
      rownames(x2) <- mg
      #colnames(x2) <- colnames(x)
    }

    if (pa_block$n_snps > 1) {
      x2 <-x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- x[rep(1, pa_block$n_snps), ]
      rownames(x) <- pa_block$snps$rsid
      rownames(snp_gene_df) <- snp_gene_df$rsid

      x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
      x <- as(x, "dgCMatrix")
      x2 <- as(x2, "dgCMatrix")
      x3 <- x2 * x
    } else {
      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label
      x2 <- as(x2, "dgCMatrix")
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
      rownames(x) <- pa_block$snps$rsid
      #snp_gene_df <- snp_gene_df
      rownames(snp_gene_df) <- snp_gene_df$rsid
      x <- matrix(as.numeric(x) * as.numeric(snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
      x <- as(x, "dgCMatrix")
      x3 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)

    }
    rm(x)
    rm(x2)
    ## add the slope of eqtls for x
    rownames(x3) <- pa_block$snps$rsid

    pa_block$x<-Matrix::crossprod( pa_block$ld_matrix_squared,x3)
    pa_block$include_in_inference <- T

    setTxtProgressBar(pb, which(names(Pathway_ld_gwas_data) == pathway) / length(Pathway_ld_gwas_data))

    return(pa_block)
  }, n.cores = n.cores)
  close(pb)

  SC_Pathway_ld_gwas_data <- Pathway_ld_gwas_data[!sapply(Pathway_ld_gwas_data, is.null)]
  message("*** Start to store the SC_CPathway_ld_gwas_data")
  SOAR::Store(SC_Pathway_ld_gwas_data)

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
