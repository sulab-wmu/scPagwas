
#' get_Pathway_sclm
#'
#' @param pa_block single pathway blocks.
#' @param pca_scCell_mat Element in Pagwas.
#' @param data_mat Element in Pagwas.
#' @param rawPathway_list Element in Pagwas..
#' @param n.cores cores for regression
#' @param backingpath file address for bk files, no "/" in the end.
#' @param snp_gene_df Element in Pagwas.
#' @param Rns Random bk names
#' @return lm result for pathway in single cell.
get_Pathway_sclm <- function(pa_block,
                             pca_scCell_mat,
                             data_mat,
                             rawPathway_list,
                             snp_gene_df,
                             n.cores=1,
                             backingpath,
                             Rns) {

  pathway <- unique(pa_block$block_info$pathway)
  x <- matrix(pca_scCell_mat[pathway, ], nrow = 1)
  rownames(x) <- pathway

  if (pa_block$n_snps == 0) {
    pa_block$include_in_inference <- F
    pa_block$x <- NULL # to make sure we totally replace previous stuffs
    return(pa_block)
  }

  mg <- intersect(rawPathway_list[[pathway]], rownames(data_mat))
  if (length(mg) == 1) {
    x2 <- matrix(data_mat[mg, ], nrow = 1)
    x2 <- x2 / (x2 + 0.0001)
    rownames(x2) <- mg
  } else {
    x2 <- biganalytics::apply(data_mat[mg, ], 2, function(ge) {
      if (sum(ge) == 0) {
        return(rep(0, length(ge)))
      } else {
        return(ge / sum(ge))
      }
    })
    rownames(x2) <- mg
  }

  x2 <- as(x2, "CsparseMatrix")

  if (pa_block$n_snps > 1) {
    x2 <- x2[pa_block$snps$label, ]
    pa_block$n_snps <- nrow(pa_block$snps)

    x <- x[rep(1, pa_block$n_snps), ]
    rownames(x) <- pa_block$snps$rsid

    #x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
    x2 <- x2 * x
  } else {
    x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
    rownames(x2) <- pa_block$snps$label
    pa_block$n_snps <- nrow(pa_block$snps)

    x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
    rownames(x) <- pa_block$snps$rsid


    #x <- matrix(as.numeric(x) * as.numeric(
    #  snp_gene_df[pa_block$snps$rsid, "slope"]
    #), nrow = 1)
    x2 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
    x2 <- as(x2, "CsparseMatrix")
  }
  pa_block$x <- as(pa_block$ld_matrix_squared %*% x2 , "CsparseMatrix")

  pa_block$include_in_inference <- T

  noise_per_snp <- pa_block$snps$se**2

  if (!is.null(pa_block$x)) {
    if (pa_block$n_snps > 2) {
      na_elements <- is.na(pa_block$y) | apply(pa_block$x, 1, function(x) {
        any(is.na(x))
      }) | is.na(noise_per_snp)

      results <- scParameter_regression(
        Pagwas_x = pa_block$x[!na_elements, ],
        Pagwas_y = pa_block$y[!na_elements],
        noise_per_snp = noise_per_snp[!na_elements],
        Rns=Rns,
        backingpath=backingpath,
        n.cores=n.cores
      )
      results[is.na(results)] <- 0
      names(results) <- colnames(data_mat)
    } else {
      results <- NULL
    }
  } else {
    results <- NULL
  }

  return(results)
}

#' scPagwas_perform_score
#' @description Get the scPagwas score for each cells
#' @param Pagwas a list data including all Intermediate result
#' @param remove_outlier Whether to remove the outlier for scPagwas score.
#' storage, you can set a large number to save time.
#' @return Pagwas result list including scPagwas.gPAS.score
#' @export

scPagwas_perform_score <- function(Pagwas,
                                   remove_outlier = TRUE) {
  options(bigmemory.allow.dimnames = TRUE)

  Pathway_sclm_results <- Pagwas$Pathway_sclm_results

  Pathway_names <- colnames(Pathway_sclm_results)

  pathway_expr <- lapply(Pathway_names, function(pa) {
    a <- intersect(Pagwas$Pathway_list[[pa]], rownames(Pagwas$data_mat))
    if (length(a) == 0) {
      return(NULL)
    } else if (length(a) == 1) {
      return(Pagwas$data_mat[a, ])
    } else {
      b <- biganalytics::apply(Pagwas$data_mat[a, ], 2, mean)
      return(b)
    }
  })
  a<-!sapply(pathway_expr, is.null)
  pathway_expr <- data.matrix(as.data.frame(pathway_expr[a]))
  Pathway_names <- Pathway_names[a]

  colnames(pathway_expr) <- Pathway_names
  pa_exp_mat <- t(Pagwas$pca_scCell_mat[Pathway_names, ]) * pathway_expr
  rm(pathway_expr)
  Pagwas$Pathway_single_results <- Pathway_sclm_results[, Pathway_names] * pa_exp_mat
  rownames(Pagwas$Pathway_single_results) <- colnames(Pagwas$pca_scCell_mat)
  message("* Get Pathways'rankPvalue for each celltypes!")
  cl <- unique((Pagwas$Celltype_anno$annotation))

  Pagwas$Pathway_single_results <- t(data.matrix(
    Pagwas$Pathway_single_results
  ))

  Pathways_rankPvalue <- lapply(cl, function(ss) {
    print(ss)
    tt <- Pagwas$Celltype_anno$annotation == ss
    PathwayrankPvalue <- scGene_rankP(Pagwas$Pathway_single_results[, tt])
    return(PathwayrankPvalue$pValueHigh)
  })

  Pagwas$scPathways_rankPvalue <- Reduce(
    function(dtf1, dtf2) cbind(dtf1, dtf2),
    Pathways_rankPvalue
  )
  Pagwas$scPathways_rankPvalue <- as.data.frame(Pagwas$scPathways_rankPvalue)
  colnames(Pagwas$scPathways_rankPvalue) <- cl
  rownames(Pagwas$scPathways_rankPvalue) <- rownames(Pagwas$Pathway_single_results)
  rm(Pathways_rankPvalue)
  message("* Get scPgwas score for each single cell")
  scPagwas.gPAS.score <- colSums(Pagwas$Pathway_single_results)
  gc()
  if (remove_outlier) {
    Pagwas$scPagwas.gPAS.score <- scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
  }
  names(Pagwas$scPagwas.gPAS.score) <- colnames(Pagwas$pca_scCell_mat)
  return(Pagwas)
}



#' scParameter_regression
#' @description Find parameter estimates for the data.
#' In terms of settings, it can handle data from at least 10,000 cells in your computer.
#' @param Pagwas_x x parameter for lm
#' @param Pagwas_y y parameter for lm
#' @param noise_per_snp noise
#' @param n.cores cores for regression
#' @param Rns the names for bk temfiles
#' @param backingpath file address for bk files, no "/" in the end.
#'
#' @export


scParameter_regression <- function(Pagwas_x,
                                   Pagwas_y,
                                   noise_per_snp,
                                   Rns,
                                   n.cores=1,
                                   backingpath) {

    backingpath<- paste0(backingpath,"/",Rns)
    if(dim(Pagwas_x)[2] <= 20000){
      Pagwas_x<-as.matrix(Pagwas_x)

    }else if(dim(Pagwas_x)[2] > 20000){
      # 将矩阵划分为n个块（按列划分）
      n<- floor(ncol(Pagwas_x)/10000)
      split_cols <- split(1:ncol(Pagwas_x), cut(1:ncol(Pagwas_x), n, labels = FALSE))
      # 逐个块处理，并将它们合并
      Pagwas_x <- do.call("cbind", lapply(split_cols, function(cols){
        as.matrix(Pagwas_x[,cols])
      }))
    }
    liear_m <- bigstatsr::big_univLinReg(
      X=bigstatsr::as_FBM(Pagwas_x,backingfile = backingpath),
      y.train=Pagwas_y,
      covar.train = bigstatsr::covar_from_df(
        data.frame(offset(noise_per_snp))
      ),
      ncores = n.cores
    )
    parameters <- liear_m$estim
      unlink(paste0(backingpath,".bk"),recursive = TRUE)
   # rm(Pagwas_x)
  return(parameters)
}

#' scPagwas_score_filter
#' @description filter the scPagwas_score for outliers.
#' @param scPagwas_score (data.frame)
#' @export

scPagwas_score_filter <- function(scPagwas_score) {
  # remove the NAN!
  if (sum(is.nan(scPagwas_score)) > 0) {
    scPagwas_score[is.nan(scPagwas_score)] <- 0
  }
  # remove the inf values!
  if (Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == Inf)] <-
      max(scPagwas_score[-which(scPagwas_score == Inf)], na.rm = TRUE)
  }
  if (-Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == -Inf)] <-
      min(scPagwas_score[-which(scPagwas_score == -Inf)], na.rm = TRUE)
  }
  lower_bound <- stats::quantile(scPagwas_score, 0.01, na.rm = TRUE)
  upper_bound <- stats::quantile(scPagwas_score, 0.99, na.rm = TRUE)

  lower_ind <- which(scPagwas_score < lower_bound)
  upper_ind <- which(scPagwas_score > upper_bound)
  scPagwas_score[lower_ind] <- lower_bound
  scPagwas_score[upper_ind] <- upper_bound

  return(scPagwas_score)
}



#' scGet_PCC
#'
#' @param scPagwas.gPAS.score score of scPagwas pathway lm
#' @param data_mat matrix for single cell data
#' @return result list including pearson correlation coefficients(PCC)
#' @export

scGet_PCC <- function(scPagwas.gPAS.score,data_mat) {
  if (!inherits(data_mat, "matrix")) {
  data_mat <- as_matrix(data_mat)
  }
  if (all(names(scPagwas.gPAS.score) == colnames(data_mat))) {
    scPagwas.gPAS.score <- data.matrix(scPagwas.gPAS.score)
    sparse_cor <- corSparse(
      X = t(data_mat),
      Y = scPagwas.gPAS.score
    )
  } else {
    scPagwas.gPAS.score <- data.matrix(scPagwas.gPAS.score)
    sparse_cor <- corSparse(
      X = t(data_mat),
      Y = scPagwas.gPAS.score
    )
  }
  rownames(sparse_cor) <- rownames(data_mat)
  colnames(sparse_cor) <- "PCC"
  sparse_cor[is.nan(sparse_cor)] <- 0
  #
  return(sparse_cor)
}


#' pa_meanexpr
#'
#' @param Pagwas
#'
#' @return Pagwas
#' @export
#'
#' @examples
pa_meanexpr <- function(Pagwas){
  Pathway_names <- names(Pagwas$Pathway_list)
  pathway_expr <- lapply(Pathway_names, function(pa) {
    a <- intersect(Pagwas$Pathway_list[[pa]], rownames(Pagwas$data_mat))
    if (length(a) == 0) {
      return(NULL)
    } else if (length(a) == 1) {
      return(Pagwas$data_mat[a, ])
    } else {
      b <- biganalytics::apply(Pagwas$data_mat[a, ], 2, mean)
      return(b)
    }
  })

  a<-!sapply(pathway_expr, is.null)
  pathway_expr <- data.matrix(as.data.frame(pathway_expr[a]))
  Pathway_names <- Pathway_names[a]

  colnames(pathway_expr) <- Pathway_names
  pa_exp_mat <- t(Pagwas$pca_scCell_mat[Pathway_names, ]) * pathway_expr
  Pagwas$pa_exp_mat<-pa_exp_mat
  return(Pagwas)
}


#' Corr_Random
#' Split the large-scale single-cell data and perform random computations of correlation. Multiple iterations of random computations ensure the accuracy of the results.
#' @param data_mat the expression matrix for single cell.
#' @param scPagwas.gPAS.score the gPas score for scPagwas
#' @param seed random seed.
#' @param random whether to random, default is TRUE.
#' @param Nrandom the number of random.
#' @param Nselect the number of single cells for each random.
#'
#' @return
#' @export
#'
#' @examples
Corr_Random<-function(data_mat,scPagwas.gPAS.score,seed=1234,random=T,Nrandom=5,Nselect=10000){
  set.seed(seed)
  merg_m<-apply(matrix(seq(Nrandom)),1,function(x){
    print(x)
    ss=sample(seq(length(scPagwas.gPAS.score)),Nselect)
    PCC <- scGet_PCC(scPagwas.gPAS.score=scPagwas.gPAS.score[ss],data_mat=data_mat[,ss])
    return(PCC)
  })
  meancor<-apply(merg_m,1,mean)
  names(meancor)<-rownames(data_mat)
  return(meancor)
}


#' Merge_gPas
#'
#' @param Pagwas
#' @param pmat_merge the merged sclm matrix
#'
#' @return
#' @export
#'
#' @examples
Merge_gPas <-function(Pagwas,pmat_merge){
  #判断Pagwas中是否包含pa_exp_mat
  if(!exists("pa_exp_mat",where = Pagwas)){
    Pagwas<-pa_meanexpr(Pagwas)
  }
  pa_exp_mat <- t(Pagwas$pca_scCell_mat) * Pagwas$pa_exp_mat
  Pathway_names<-intersect(colnames(pmat_merge),colnames(Pagwas$pa_exp_mat))
  Pathway_single_results <- pmat_merge[, Pathway_names] * pa_exp_mat[, Pathway_names]
  message("* Get scPgwas score for each single cell")
  scPagwas.gPAS.score <- rowSums(Pathway_single_results)
  scPagwas.gPAS.score <- scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
  return(scPagwas.gPAS.score)
}


#' Random_PCC
#'
#' @param gPas gPas score
#' @param datamat data matrix for single cell
#' @param seed random seed
#' @param random_times random times
#' @param select_num the cells for each random
#'
#' @return
#' @export
#'
#' @examples
Random_PCC<- function(gPas,datamat,seed=1234,random_times=100,select_num=10000){
  message("Start to run the corSparse function in random!")
  set.seed(seed)
  sparse_cor_list<-list()
  for (j in 1:random_times) {
    print(paste0("Randome Times: ",j))
    index<-sample(1:ncol(datamat),select_num)
    gPas_select <- data.matrix(gPas[index])
    sparse_cor <- scPagwas::corSparse(
      X = t(datamat[,index]),
      Y = gPas_select
    )
    colnames(sparse_cor) <- "PCC"
    sparse_cor[is.nan(sparse_cor)] <- 0
    sparse_cor[is.na(sparse_cor)] <- 0
    sparse_cor_list[[j]]<-unlist(sparse_cor[,1])
  }
  sparse_cor_list<-as.data.frame(sparse_cor_list)
  sparse_cor<- apply(sparse_cor_list,1,function(x) mean(x, na.rm = TRUE))
  return(data.matrix(sparse_cor))
}

#' scGet_PCC2
#'
#' @param Pagwas a result list including scPagwas.gPAS.score and data_mat
#' @param random_times the random times
#'
#' @return
#' @export
#'
#' @examples
#'
scGet_PCC2<-function(Pagwas,random_times=100){
  scPagwas.gPAS.score <- Pagwas$scPagwas.gPAS.score
  pcc <- Random_PCC(gPas=Pagwas$scPagwas.gPAS.score,datamat=Pagwas$data_mat,
                    seed=1234,random_times=random_times,
                    select_num=floor(ncol(Pagwas$Pathway_single_results)/5))
  colnames(pcc)<-'PCC'

  # 提取前5%的scPagwas.gPAS.score
  top5_idx <- order(scPagwas.gPAS.score, decreasing = TRUE)[1:(0.05 * length(scPagwas.gPAS.score))]

  # 根据索引获取名称
  top5_cellnames <- names(scPagwas.gPAS.score)[top5_idx]

  # 在scPagwas.gPAS.score的名称中判断是否在top5_cellnames中
  top5_idx2 <- names(scPagwas.gPAS.score) %in% top5_cellnames

  # 对Pagwas$data_mat的每一行应用wilcox.test并进行bonferroni校正
  p_list <- apply(Pagwas$data_mat, 1, function(x) {
    p.adjust(wilcox.test(x[top5_idx2], x[!top5_idx2], alternative = "greater")$p.value, method = "bonferroni")
  })

  # 创建pcc数据框并添加pvalue
  pcc <- data.frame(pcc, pvalue = p_list)

  # 对pvalue进行bonferroni校正
  pcc <- within(pcc, {
    adj_pvalue <- p.adjust(pvalue, method = "bonferroni")
    adj_logp <- -log10(adj_pvalue)
    adj_logp <- ifelse(is.finite(-log10(adj_pvalue)), -log10(adj_pvalue), max(adj_logp, na.rm = TRUE) + 1)
  })

  pcc$weight_pcc <- pcc$adj_logp * pcc$PCC
  Pagwas$PCC <- pcc
  return(Pagwas)
}
