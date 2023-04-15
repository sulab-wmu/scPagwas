#' @title Get_CorrectBg_p
#' @description Get_CorrectBg_p
#' @param Single_data Seurat object
#' @param scPagwas.TRS.Score raw score
#' @param iters_singlecell number of iterations
#' @param n_topgenes number of genes
#' @param scPagwas_topgenes gene list
#' @return correct_pdf
#' @export
Get_CorrectBg_p<-function(Single_data,scPagwas.TRS.Score, iters_singlecell,n_topgenes,scPagwas_topgenes){
    gene_matrix <- GetAssayData(Single_data, slot = "data", assay = 'RNA')
    mat_ctrl_raw_score <- matrix(0, nrow =ncol(gene_matrix), ncol = iters_singlecell)
    #随机100次抽取基因集，计算每个样本的control分数
    dic_ctrl_list <- list()
    pb <- txtProgressBar(style = 3)
    for (i in 1:iters_singlecell) {
        #随机抽取基因集
        #print(i)
        set.seed(i)
        dic_ctrl_list[[i]] <- sample(rownames(Single_data), n_topgenes)
        #计算每个样本的control分数
        Single_data<-Seurat::AddModuleScore(Single_data,
            assay = 'RNA',
            list(dic_ctrl_list[[i]]),
            name = c("contr_genes")
            )
        #将每个样本的control分数放入矩阵中
        mat_ctrl_raw_score[, i] <-Single_data$contr_genes1
        Single_data$contr_genes1<-NULL
        setTxtProgressBar(pb, i/iters_singlecell)
    }
    close(pb)
    #获得基因表达的方差
    genes<-intersect(rownames(Single_data),rownames(gene_matrix))
    scPagwas_topgenes<-intersect(scPagwas_topgenes,genes)
    gene_matrix<-gene_matrix[genes,]
    df_gene <- data.frame(
    gene = genes,
    var = apply(gene_matrix,1,var)
    )
    rownames(df_gene) <- df_gene$gene
    v_var_ratio_c2t <- rep(1, iters_singlecell)
    # For mean_var matched control genes and raw scores computed as weighted average,
    # estimate variance ratio assuming independence.
    for (i_ctrl in 1:iters_singlecell) {
    v_var_ratio_c2t[i_ctrl] <- sum(df_gene[dic_ctrl_list[[i_ctrl]], "var"])
    }
    v_var_ratio_c2t <- v_var_ratio_c2t / sum(df_gene[scPagwas_topgenes, "var"])

    correct_pdf<-correct_background(
            scPagwas.TRS.Score,
            mat_ctrl_raw_score,
            v_var_ratio_c2t
        )
    rownames(correct_pdf)<-colnames(Single_data)
    return(correct_pdf)
}

#' @title correct_background
#' @description correct_background
#' @param scPagwas.TRS.Score raw score
#' @param mat_ctrl_raw_score control raw score
#' @param v_var_ratio_c2t variance ratio
#' @return correct_pdf
#' @export
correct_background <- function(scPagwas.TRS.Score, mat_ctrl_raw_score, v_var_ratio_c2t) {
    ind_zero_score <- scPagwas.TRS.Score == 0
    ind_zero_ctrl_score <- mat_ctrl_raw_score == 0
    #
    scPagwas.TRS.Score <- scPagwas.TRS.Score - mean(scPagwas.TRS.Score)
    mat_ctrl_raw_score <- mat_ctrl_raw_score - colMeans(mat_ctrl_raw_score)
    mat_ctrl_raw_score <- mat_ctrl_raw_score / sqrt(v_var_ratio_c2t)
    v_mean <- rowMeans(mat_ctrl_raw_score)
    v_std <- apply(mat_ctrl_raw_score,1,sd)
    v_norm_score <- scPagwas.TRS.Score
    v_norm_score <- (v_norm_score - v_mean) / v_std
    mat_ctrl_norm_score <- t(t(mat_ctrl_raw_score) - v_mean) / v_std
    v_norm_score <- v_norm_score - mean(v_norm_score)
    mat_ctrl_norm_score <- mat_ctrl_norm_score - colMeans(mat_ctrl_norm_score)
    norm_score_min <- min(min(v_norm_score), min(mat_ctrl_norm_score))
    v_norm_score[ind_zero_score] <- norm_score_min - 1e-3
    mat_ctrl_norm_score[ind_zero_ctrl_score] <- norm_score_min
    flatten_mat_ctrl_norm_score <- as.vector(t(mat_ctrl_norm_score))
        pooled_p = get_p_from_empi_null(v_norm_score, flatten_mat_ctrl_norm_score)
    #校正pvalue
    adj_p <- p.adjust(pooled_p, method = "BH")

    nlog10_pooled_p <- -log10(pooled_p)
    pooled_z <- -qnorm(pooled_p)
    pooled_z <- pmin(pmax(pooled_z, -10), 10)
    correct_pdf<-data.frame(
        raw_score = scPagwas.TRS.Score,
        pooled_p = pooled_p,
        adj_p = adj_p,
        nlog10_pooled_p = nlog10_pooled_p,
        pooled_z = pooled_z
    )
  return(correct_pdf)
}

#' @title get_p_from_empi_null
#' @description get_p_from_empi_null
#' @param v_t raw score
#' @param v_t_null control raw score
#' @return p
#' @export
get_p_from_empi_null <- function(v_t, v_t_null) {
  v_pos <- sapply(v_t, function(x) sum(v_t_null <= x))
  v_p <- (length(v_t_null) - v_pos + 1) / (length(v_t_null) + 1)
  return(v_p)
}

#' @title Merge_celltype_p
#' @description Merge_celltype_p
#' @param single_p single p
#' @param celltype celltype
#' @return celltype_p
#' @export
Merge_celltype_p<-function(single_p,celltype){
    celltype_p<-data.frame(celltype=celltype,pvalue=single_p)
    celltype_p<-aggregate(pvalue~celltype,celltype_p,p_merge)
    return(celltype_p)
}
#' @title p_merge
#' @description p_merge
#' @param pvalues pvalues
#' @return p_total
#' @export
p_merge<-function(pvalues){
zvalues <- -sqrt(2) * qnorm(pvalues/2)
ztotal <- mean(zvalues)
p_total <- pnorm(-abs(ztotal))
return(p_total)
}
