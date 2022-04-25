
#library(gtools)
#library(verification)
#library(doParallel)
#library(foreach)
#library(magrittr)

scPagwas_celltype_enrichment<- function(Pagwas){

celltypelist<-tapply(Pagwas$Celltype_anno$cellnames, Pagwas$Celltype_anno$annotation, function(x){
  return(x)
})

m<-data.matrix(Pagwas$scPagwas_score)

## number of sample columns
Ns <- 1
## number of genes/ptm sites
Ng <- nrow(m)
## Extract input file name
#input.file.prefix <-  sub('.*/(.*)\\.gct$', '\\1', input.ds)

## ####################################################
##            Sample normalization
## -take care of missing values already here
## ####################################################

if (sample.norm.type == "rank") {
  m <- apply(m, 2, function(x){
    x.rank=rep(NA, length(x))
    keep.idx=which(!is.na(x))
    x.rank[keep.idx ]=rank(x[keep.idx], ties.method="average")
    return(x.rank)
  })
  m <- 10000*m/Ng

} else if (sample.norm.type == "log.rank") {
  m <- apply(m, 2, function(x){
    x.rank=rep(NA, length(x))
    keep.idx=which(!is.na(x))
    x.rank[keep.idx ]=rank(x[keep.idx], ties.method="average")
    return(x.rank)
  })
  m <- log(10000*m/Ng + exp(1))

} else if (sample.norm.type == "log") {
  m[m < 1] <- 1
  m <- log(m + exp(1))

} ##else if (sample.norm.type == "none") {

be<-Pagwas$bootstrap_results$bias_corrected_estimate[-1]
correl.be <- (be - mean(be))/sd(be)
correl.be <-(correl.be-min(correl.be))/(max(correl.be)-min(correl.be))

tmp<-lapply(names(celltypelist), function(i){
  weight<-correl.be[i]
  celltypel<-celltypelist[[i]]
  print(i)
  OPAM <- project.geneset (data.array = m,
                           gene.names = names(Pagwas$scPagwas_score),#矩阵的基因名字
                           n.cols = Ns,
                           n.rows= Ng,
                           weight =weight ,
                           statistic = "area.under.RES",
                           gene.set = celltypel,#基因集合的基因名
                           nperm = 1000,
                           correl.type = "rank",
                           gene.set.direction =  NULL,
                           min.overlap = 10,
                           size.G.current =length(celltypel) )
  return(OPAM)
})

## #######################################
## extract scores and pvalues
##  and generate matrices
tmp.pval <- lapply(tmp, function(x)x$p.val.vector)
pval.matrix <- matrix(unlist(tmp.pval), byrow=T, nrow=Ng)

if (output.score.type == "ES"){
  tmp.es <- lapply(tmp, function(x)x$ES.vector)
  score.matrix <- matrix(unlist(tmp.es), byrow=T, nrow=Ng)
}
if (output.score.type == "NES"){
  tmp.nes <- lapply(tmp, function(x)x$NES.vector)
  score.matrix <- matrix(unlist(tmp.nes), byrow=T, nrow=N.gs)
}

  fdr.matrix <- matrix ( p.adjust(unlist (pval.matrix.2), method='fdr'),
                           ncol=ncol(pval.matrix.2))
  locs<-1:length(tmp)

  sample.descs.tmp <- data.frame(signature.score=as.numeric(score.matrix[locs, ]),
                                 signature.pvalue=as.numeric(pval.matrix[locs, ]),
                                 signature.fdr.pvalue=as.numeric( fdr.matrix[locs, ]),
                                 stringsAsFactors=F)

 # fn <- paste(getwd(),'/signature_gct/', make.names(sig.name), sep='')
  #if((nchar(fn)+15) > max.nchar.file.path){
  #  fn <- paste( unlist(strsplit(fn,''))[1:(max.nchar.file.path-15)], collapse='')
  #  cat(nchar(fn), ' ', fn ,'\n')
 # }
 # write.gct(gct.tmp, ofile=sub('.*/', 'signature_gct/', fn), appenddim = T)

  Pagwas$celltypes_enrichment_result<-sample.descs.tmp
  return(Pagwas)

 }
## ####################################################
##
##   function executed per gene set
##
## ####################################################
project.geneset <- function (data.array,
                             gene.names,
                             n.cols,
                             n.rows,
                             weight = 0,
                             statistic = "Kolmogorov-Smirnov",   ## alternatives: "Kolmogorov-Smirnov", "area.under.RES"
                             gene.set,
                             nperm = 200,
                             correl.type  = "rank",              ## "rank", "z.score", "symm.rank"
                             gene.set.direction=NULL,             ## direction of regulation; has to be in same order than 'gene.set'
                             min.overlap,
                             size.G.current
) {

  ## end function 'gsea.score'

  ## fast implementation of GSEA)
  ES.vector <- NES.vector <- p.val.vector <- rep(NA, n.cols)
  correl.vector <- rep(NA, n.rows)

  ## vectors to return number of overlapping genes/sites
  OL.vector <- OL.numb.vector <- OL.perc.vector <- rep(NA, n.cols)

  ## Compute ES score for signatures in each sample
  phi <- array(NA, c(n.cols, nperm))

  ## list to store random walk accross samples
  random.walk <- vector('list', n.cols)

  ## locations of gene set in input data (before ranking/ordering)
  ## 'gene.names' is in the same order as the input data
  gene.names.all <- gene.names
  gene.set.all <- gene.set
  gene.set.size <- length(gene.set)
  gene.set.direction.all <- gene.set.direction
  n.rows.all <- n.rows

  ## loop over columns/samples
 # for (sample.index in 1:n.cols) {
  sample.index<-1
    ## ######################################################
    ##         handle missing values appropriatly
    ## ######################################################

    ## reset values
    gene.names <- gene.names.all
    gene.set <- gene.set.all
    gene.set.direction <- gene.set.direction.all
    n.rows <- n.rows.all

    ## extract expression values of current sample
    data.expr <- data.array[ , sample.index]

    ## index of missing values
    na.idx <- which( is.na(data.expr) | is.infinite(data.expr) )

    ## if there are missing values ...
    if( length(na.idx) > 0){

      ## update input data + gene names
      ## those are in the same order
      data.expr <- data.expr[-na.idx]
      gene.names <- gene.names[-na.idx]

      ## update numbers
      n.rows=length(data.expr)

      ## extract gene set members present in data
      gene.set.no.na <- which(gene.set %in% gene.names)

      ## update gene sets + gene set directions
      ## those are in the same order
      gene.set <- gene.set[ gene.set.no.na ]
      gene.set.direction <- gene.set.direction[ gene.set.no.na ]

      ##cat('gene set overlap:', length(gene.set), '\n', file=log.file, append=T)
      #cat('gene set overlap:', sum(gene.names %in% gene.set), '\n', file=log.file, append=T)
    } ## end if missing values are present

    ## ###############################################
    ## if there is NOT sufficient overlap...

      ## locations of gene set in input data (before ranking/ordering)
      ## 'gene.names' is in the same order as the input data
      ## gene.set2 <- match(gene.set, gene.names)
      ## gene.set2 <- which( gene.names %in% gene.set ) ## to work with redundant gene names, e.g. phospho-data
      ## this takes care about redundant gene lists, the approach above does not return the locations of 'gene.set' in 'gene.names' but the indices of 'gene.names' present in 'gene.set'
      gene.set2 <- as.numeric( unlist(sapply( gene.set, function(x) which(gene.names == x) )))

      # save(gene.set2, gene.set, gene.names, file='tmp.RData')

      ## order of ranks, list is now ordered, elements are locations of the ranks in
      ## original data,
      gene.list <- order(data.expr, decreasing=T)

      ## ##############################################
      ## calculate enrichment score
      GSEA.results <- gsea.score (ordered.gene.list=gene.list, gene.set2=gene.set2,
                                  statistic =statistic,
                                  weight=weight,
                                  n.rows=n.rows,
                                  correl.type=correl.type,
                                  gene.set.direction,
                                  data.expr=data.expr)
      ES.vector[sample.index] <- GSEA.results$ES

      ## overlap between gene set and data

        OL.vector[ sample.index ] <- paste( unique( gene.names[gene.list][GSEA.results$indicator ]) , collapse='|')
        OL.numb.vector[sample.index ] <- length( unique( gene.names[gene.list][GSEA.results$indicator ]) )
        OL.perc.vector[sample.index ] <- round(100*( OL.numb.vector[ sample.index ]/size.G.current  ),1)


      ## store the gsea results to
      random.walk[[sample.index]] <- GSEA.results

      ## ##############################################
      ## no permutations: - ES and NES are the same
      ##                  - all p-values are 1

        ES.tmp = sapply(1:nperm,  function(x) gsea.score(sample(1:n.rows),
                                                         statistic =statistic,
                                                         gene.set2, weight,
                                                         n.rows, correl.type=correl.type,
                                                         gene.set.direction,
                                                         data.expr)$ES)
        phi[sample.index, ] <- unlist(ES.tmp)

        ## #######################################################
        ## calculate NES and p-values
        if (ES.vector[sample.index] >= 0) {
          pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
          if (length(pos.phi) == 0) pos.phi <- 0.5
          pos.m <- mean(pos.phi)
          NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
          s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
          p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)

        } else {
          neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
          if (length(neg.phi) == 0) neg.phi <- 0.5
          neg.m <- mean(neg.phi)
          NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
          s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
          p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
        }
      ## end do permutations


  return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))

} ## end function 'project.geneset'


## #############################################################################
##
## function to calculate GSEA enrichment score
## - apply correlation scheme and weighting
## - calculate ES
## ############################################################################
gsea.score <- function (ordered.gene.list,
                        gene.set2,
                        statistic,
                        weight,
                        n.rows,
                        correl.type,
                        gene.set.direction,
                        data.expr) {

  ##################################################################
  ## function to calculate ES score
  score <- function(max.ES, min.ES, RES, gaps, valleys, statistic){
    ## KM
    if( statistic == "Kolmogorov-Smirnov" ){
      if( max.ES > -min.ES ){
        ES <- signif(max.ES, digits=5)
        arg.ES <- which.max(RES)
      } else{
        ES <- signif(min.ES, digits=5)
        arg.ES <- which.min(RES)
      }
    }
    ## AUC
    if( statistic == "area.under.RES"){
      if( max.ES > -min.ES ){
        arg.ES <- which.max(RES)
      } else{
        arg.ES <- which.min(RES)
      }
      gaps = gaps+1
      RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES) - c(valleys,0) ) * (gaps)
      ES = sum(RES)
    }
    return(list(RES=RES, ES=ES, arg.ES=arg.ES))
  } ## end function score


  ## #######################################
  ## weighting
  ## #######################################
  if (weight == 0) {

    correl.vector <- rep(1, n.rows)

  } else if (weight > 0) {
    ## if weighting is used (weight > 0), bring
    ## 'correl.vector' into the same order
    ## as the ordered gene list
    if (correl.type == "rank") {
      ##correl.vector <- data.array[ordered.gene.list, sample.index]
      correl.vector <- data.expr[ordered.gene.list]

    } else if (correl.type == "symm.rank") {
      ##correl.vector <- data.array[ordered.gene.list, sample.index]
      correl.vector <- data.expr[ordered.gene.list]

      correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)],
                              correl.vector,
                              correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)])
    } else if (correl.type == "z.score") {
      ##x <- data.array[ordered.gene.list, sample.index]
      x <- data.expr[ordered.gene.list]
      correl.vector <- (x - mean(x))/sd(x)
    }
  }

  ## length of gene list - equals number of rows in input matrix
  N = length(ordered.gene.list)

    Nh <- length(gene.set2)
    Nm <-  N - Nh

    ## #####################################
    ## match gene set to data
    tag.indicator <- sign(match(ordered.gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
    ## positions of gene set in ordered gene list
    ind = which(tag.indicator==1)
    ## 'correl.vector' is now the size of 'gene.set2'
    correl.vector <- abs(correl.vector[ind])^weight
    ## sum of weights
    sum.correl = sum(correl.vector)

    #########################################
    ## determine peaks and valleys
    ## divide correl vector by sum of weights
    up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
    gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
    down = gaps/Nm

    RES = cumsum(c(up,up[Nh])-down)
    valleys = RES[1:Nh]-up

    max.ES = max(RES)
    min.ES = min(valleys)

    ## calculate final score
    score.res <- score(max.ES, min.ES, RES[1:Nh], gaps, valleys,
                       statistic=statistic)

    ES <- score.res$ES
    arg.ES <- score.res$arg.ES
    RES <- score.res$RES

    gsea.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=up, step.down=1/Nm)


  return (gsea.results)
}


