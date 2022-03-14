
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

#' mean_expr
#' @description To get the mean expression for single cell data.
#' @param Pagwas Pagwas data list
#'
#' @return
#'
#' @examples
#' library(scPagwas)
#' Pagwas$merge_scexpr <- mean_expr(Pagwas)
#'
mean_expr <- function(Pagwas){
  datamat<-Seurat::GetAssayData(object =Pagwas$Single_data, slot = "data")
  merge_scexpr<-lapply(unique(Pagwas$Celltype_anno$annotation),function(x){
    an<-Pagwas$Celltype_anno$cellnames[Pagwas$Celltype_anno$annotation %in% x]
    tryCatch(
      {
        datamat<-as.matrix(datamat)
      }, error = function(e) {
        datamat<-as_matrix(datamat)

      })

    return(apply(datamat[,an],1,mean))
  })
  rm(datamat)
  merge_scexpr<-as.data.frame(merge_scexpr)
  colnames(merge_scexpr)<-unique(Pagwas$Celltype_anno$annotation)
  return(merge_scexpr)
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
Snp2Gene <- function(snp, refGene, marg = 5000) {
  snp_GR <- GenomicRanges::GRanges(snp[, "chrom"],
                                   IRanges::IRanges(as.numeric(snp[, "pos"]), as.numeric(snp[, "pos"])),
    name = snp[, "rsid"]
  )

  gene_GR <- GenomicRanges::GRanges(refGene[, "chrom"],
    IRanges(
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


#' as_matrix
#'
#' @param mat matrix
#'
#' @return
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
