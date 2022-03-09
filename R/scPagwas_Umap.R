
#' scPagwas_Visualization
#' @description Visualize the scPagwas score in Umap and Tsne.
#'
#' @param scPagwas
#' @param Single_data Single_data in seruat format ,the save with scPagwas_main(), you'd better to run reduction of UMAP AND TSNE
#' @param Reduction (logical) default is FALSE. Whether to run the Reduction for Single_data.If you are do it before,ignore it.
#' @param assay (character)"RNA" or "SCT", It depens on the Single_data.
#' @param files (character)the file for save the figures.
#' @param FigureType (character)"tsne" or "umap
#' @param cellpercent (numeric),default is 0.1, Threshold for pecent(<1) of Positive cells for level of scPagwas_score.
#' @param width (numeric)Figure width
#' @param height (numeric)Figure height
#' @param size (numeric)size for scatters
#' @param title (character)Figure title
#' @param lowColor (character)Color for low scPagwas score
#' @param highColor (character)Color for high scPagwas score
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(scRNAexample)
#' scPagwas_Visualization(scPagwas=Pagwas,
#'                        Single_data=scRNAexample,
#'                        Reduction=FALSE,
#'                        assay="SCT",
#'                        cellpercent= 0.1,
#'                        files="scpagwas_pbc_False",
#'                        FigureType="tsne",
#'                        width=7,
#'                        height=7,
#'                        lowColor="#000957",highColor="#EBE645",
#'                        size=0.5,
#'                        title="scPagwas_score"
#'                                  )
#'
scPagwas_Visualization<-function(scPagwas=Pagwas,
                                 Single_data,
                                 Reduction=FALSE,
                                 assay="SCT",
                                 cellpercent= 0.1,
                                 files="scpagwas_pbc_False",
                                 FigureType="tsne",
                                 width=7,
                                 height=7,
                                 lowColor="#000957",highColor="#EBE645",
                                 size=0.5,
                                 title="scPagwas_score"
){


   suppressMessages(require(ggnewscale))
   suppressMessages(require(ggtext))
   suppressMessages(require(ggrepel))

  if(!is.null(scPagwas$scPagwas_score)){
  scPagwas_score<-scPagwas$scPagwas_score
  }else{
    stop("ERROR: scPagwas_score is NULL. scPagwas_score can be calsulated by scPagwas_perform_score function!")
  }

  if(sum(colnames(Single_data) %in% scPagwas$sclm_results$cellnames) != length(scPagwas$sclm_results$cellnames)){
    message("There is blank , '-' or '+' within cell names!")
    colnames(Single_data) <- stringr::str_replace_all(colnames(Single_data)," ",".")
    colnames(Single_data) <- stringr::str_replace_all(colnames(Single_data),"\\+",".")
    colnames(Single_data) <- stringr::str_replace_all(colnames(Single_data),"-",".")
  }

  cell_id<-intersect(colnames(Single_data) %in% scPagwas_score$cellid)
  Single_data <- Single_data[,cell_id]
  scPagwas_score<-scPagwas_score[cell_id,]
  Single_data$scPagwas_score <- scPagwas_score$scPagwas_score

 if(Reduction){

    Single_data <- RunPCA(object = Single_data, assay = assay, npcs = 50)
    if(FigureType=="tsne"){
      Single_data <- RunTSNE(object = Single_data,assay = assay, reduction = "pca", dims = 1:50)
    }else if(FigureType=="umap"){
      Single_data <- RunUMAP(object = Single_data, assay =assay, reduction = "pca", dims = 1:50)
    }else{
      stop("ERROR:: The Reduction is TRUE, but no correct FigureType is choose, either tsne or umap")
    }

  }

  if(!file.exists(paste0("./",files))){
    dir.create(files)
  }

    #num<-length(unique(as.vector(Idents(Single_data))))
    Single_data<-Single_data[,!is.na(Single_data$scPagwas_score)]

    thre<- sort(Single_data$scPagwas_score,decreasing=T)[nrow(Single_data)*cellpercent]

  if(FigureType=="umap"){

     all_fortify_can<-fortify.Seurat.umap(Single_data)
     plot_scPagwas_score <- ggscatter(all_fortify_can, x = "UMAP_1", y = "UMAP_2",
                            color = "scPagwas_score",fill= "scPagwas_score",size=size,title=title,
                            repel = TRUE)+ umap_theme()+
                            scale_fill_gradient(low=lowColor,high=highColor)+
                            scale_color_gradient(low=lowColor,high=highColor)+
                            theme(aspect.ratio=1)

   pdf(file=paste0("./",files,"/scPagwas_score_umap.pdf"),height=height,width=width)
   print(plot_scPagwas_score)
   dev.off()

   plots_sigp1 <- ggplot() +
                  geom_point(data = all_fortify_can[all_fortify_can$scPagwas_score <= thre,],
                             aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8, color = "gray") +
                  umap_theme() +
                  new_scale_color() +
                  geom_point(data = all_fortify_can[all_fortify_can$scPagwas_score > thre,],
                             aes(x = UMAP_1, y = UMAP_2),color= "#F90716", size = .2) +
                  umap_theme()+
                  new_scale_color()+
                  ggtitle(paste0("Top ",cellpercent*100 ,"% cells"))

  pdf(file=paste0("./",files,"/scPagwas_TOP",cellpercent,"_umap.pdf"),height=height,width=width)
  print(plots_sigp1)
  dev.off()


  }


if(FigureType=="tsne"){

   all_fortify_can<-fortify.Seurat.tsne(Single_data)

   plot_scPagwas_score <- ggscatter(all_fortify_can, x = "TSNE_1", y = "TSNE_2",
                                    color = "scPagwas_score",fill= "scPagwas_score",size=0.5,
                                    repel = TRUE)+ umap_theme()+
                          scale_fill_gradient(low=lowColor,high=highColor)+
                          scale_color_gradient(low=lowColor,high=highColor)+
                          ggtitle("scPagwas_score") +theme(aspect.ratio=1)

   pdf(file=paste0("./",files,"/scPagwas_score_tsne.pdf"),height=height,width=width)
   print(plot_scPagwas_score)
   dev.off()

   plots_sigp1 <- ggplot() +
                  geom_point(data = all_fortify_can[all_fortify_can$scPagwas_score <= thre,],
                             aes(x = TSNE_1, y = TSNE_2), size = 0.2, alpha = 0.8, color = "gray") +
                  umap_theme() +
                  new_scale_color() +
                  geom_point(data = all_fortify_can[all_fortify_can$scPagwas_score > thre,],
                             aes(x = TSNE_1, y = TSNE_2),color= "#F90716", size = .2) +
                  umap_theme()+
                  new_scale_color()+
                  ggtitle(paste0("Top ",cellpercent*100, "% cells"))

   pdf(file=paste0("./",files,"/scPagwas_TOP",cellpercent,"_tsne.pdf"),height=height,width=width)
   print(plots_sigp1)
   dev.off()

  }
 # return(all_fortify_can)
}


#' umap_theme
#' @description ggplot2 theme for umap plot
#' @return

umap_theme <- function(){
  theme_grey() %+replace%
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank())
}




#' fortify.Seurat.umap
#' @description set data frame to ggplot
#' @param x Seurat format
#'
#' @return

fortify.Seurat.umap <- function(x){
  xy1 <- as.data.frame(Embeddings(x, reduction = "umap"))
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)

  return(cbind(xy1, as.data.frame(x@meta.data)))
}


#' fortify.Seurat.tsne
#' @description set data frame to ggplot
#' @param x seruat
#'
#' @return

fortify.Seurat.tsne <- function(x){
 xy2 <- as.data.frame(Embeddings(x, reduction = "tsne"))
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)

  return(cbind(xy2, as.data.frame(x@meta.data)))
}
