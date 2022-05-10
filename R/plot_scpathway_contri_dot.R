#' plot_scpathway_dot
#'
#' @param Pagwas
#' @param topn_path_celltype number specific pahtways to select for each celltypes.
#' @param filter_p 0.05 threshold for pvalues
#' @param display_max_sizes Boolean : Display max shape size behind each shape ? (Default=TRUE)
#' @param size_var
#' #'   If numeric : Column/List index which control shape sizes. This column/element has to be numeric.
#' #'   Can also be a column/element name or a vector of the same size than the input dataset.
#' #'   Set to NA if you don't want to control shape size.
#' @param col_var
#' #'   If numeric : Column/List index which control shape colors.
#' #'   Can also be a column/element name or a vector of the same size than the input dataset.
#' #'   Set to NA if you don't want to control shape color.
#' @param shape.scale Scale the size of the shapes, similar to cex.
#' @param cols.use 1 color or a vector containing multiple colors to color shapes.
#' @param dend_x_var A vector containing Column/List indexes or Column/List names to compute the x axis dendrogramm.
#' @param dist_method The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param hclust_method The agglomeration method to be used. This must be one of "single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid" or "median".
#' @param do_plot Boolean : whether to print the plot
#' @param figurenames
#' @param width
#' @param height
#' @param ... other parameters for plot_scpathway_contri_dot
#'
#' @return
#' @export
#'
#' @examples
plot_scpathway_dot<-function(Pagwas,
                             topn_path_celltype=10,
                             filter_p=0.05,
                             max_logp=10,
                             display_max_sizes=F,
                             size_var = "lm_beta",
                             col_var="logrankPvalue",
                             shape.scale = 8,
                             cols.use=c("lightgrey", "#E45826"),
                             dend_x_var = "lm_beta",
                             dist_method="euclidean",
                             hclust_method="ward.D",
                             do_plot = F,
                             figurenames = NULL,
                             width = 7,
                             height = 7,
                             ...){
  #library(tidyverse)
  #library("rhdf5")
  #Pathway_lm_score <- data.matrix(Pagwas$pca_cell_df) * t(data.matrix(Pagwas$Pathway_lm_results))
  #Pathway_lm_score<- as.data.frame(Pathway_lm_score)
  scPathrankP<- -log10(Pagwas$scPathways_rankPvalue+1e-20)

  top_function <-function(para_mat,n_path_to_keep=5){
    para_mat$path <- rownames(para_mat)
    #Only keep genes with a unique name and tidy data.
    para_mat <- para_mat %>% add_count(path) %>%
      filter(n==1) %>%
      select(-n) %>%
      gather(key = celltypes,value=paras,-path) %>%
      as_tibble()

    para_mat <- para_mat %>%
      group_by(celltypes) %>%
      mutate(paras_sum_mean=(paras*nrow(scPathrankP))/sum(paras))

    para_mat<- para_mat %>%
      group_by(path) %>%
      mutate(specificity=(paras_sum_mean*nrow(scPathrankP))/sum(paras_sum_mean)) %>%
      ungroup()

    d_spe<-para_mat %>% filter(paras_sum_mean>1.4)
    d_spe <- d_spe %>% group_by_("celltypes") %>% top_n(.,n_path_to_keep,specificity)
    return(d_spe)
  }


  spe2<-top_function(para_mat=scPathrankP,n_path_to_keep=topn_path_celltype)
  spe<-unique(spe2$path)

  #########

  #Pathway_lm<- as.data.frame(Pathway_lm_score[spe,])
  scPathrankP<-scPathrankP[spe,]

  scPathrankP[scPathrankP>max_logp]<-max_logp
  scPathrankP[scPathrankP< -log10(filter_p)]<-0

  #Pathway_lm$pathways<-rownames(Pathway_lm)
  scPathrankP$pathways<-rownames(scPathrankP)
  # Prep for ggplots
  #gg_lm = reshape2::melt(Pathway_lm,id.vars="pathways",variable.name="celltypes",value.name ="lm_beta")
  gg_rankp = reshape2::melt(scPathrankP,id.vars="pathways",variable.name="celltypes",value.name ="logrankPvalue")
  #gg_dot<-merge(gg_lm,gg_rankp)
  p<-plot_scpathway_contri_dot(data.to.plot = gg_rankp,
                           display_max_sizes=F,
                           size_var = "logrankPvalue",
                           col_var="logrankPvalue",
                           shape.scale = 8,
                           cols.use=c("lightgrey", "#E45826"),
                           dend_x_var = "logrankPvalue",
                           dist_method="euclidean",
                           hclust_method="ward.D",
                           ...
  )
  if (do_plot) print(p)
  if(!is.null(figurenames)){
    pdf(file = figurenames, width = width, height = height)
    print(plot_scpathway_contri_dot(data.to.plot = gg_rankp,
                          display_max_sizes=F,
                          size_var = "logrankPvalue",
                          col_var="logrankPvalue",
                          shape.scale = 8,
                          cols.use=c("lightgrey", "#E45826"),
                          dend_x_var = "logrankPvalue",
                          dist_method="euclidean",
                          hclust_method="ward.D",
                          ...
    ))
    dev.off()
  }

}


#' Dot-plot - Pacman-plot function
#'
#' Create dotplots to represent two discrete factors (x & y) described by several other factors. Each combination of the two discrete factors (x & y) can be described with : 1 continuous factor (setting shape size), 3 continuous or discrete factors (setting shape type, shape color and text on shape).
#'
#' @encoding UTF-8
#' @param data.to.plot Input data. Can be a list or a data.frame.
#'   If data.frame : Column 1 = x axis (Factor); Col2= y axis (Factor).
#'   If list : x and y axis are fixed by row and col names of list elements.
#' @param size_var
#'   If numeric : Column/List index which control shape sizes. This column/element has to be numeric.
#'   Can also be a column/element name or a vector of the same size than the input dataset.
#'   Set to NA if you don't want to control shape size.
#' @param col_var
#'   If numeric : Column/List index which control shape colors.
#'   Can also be a column/element name or a vector of the same size than the input dataset.
#'   Set to NA if you don't want to control shape color.
#' @param text_var
#'   If numeric : Column/List index which control text to add on shapes.
#'   Can also be a column/element name or a vector of the same size than the input dataset.
#'   Set to NA if you don't want to add text.
#' @param size_legend Custom name of shape legend.
#' @param col_legend Custom name of shape color.
#' @param cols.use 1 color or a vector containing multiple colors to color shapes.
#'   If coloring is continuous, default colors are taken from a "lightgrey" to "blue" gradient.
#'   If coloring is discrete, default colors are taken from the default ggplot2 palette.
#' @param shape.scale Scale the size of the shapes, similar to cex.
#' @param display_max_sizes Boolean : Display max shape size behind each shape ? (Default=TRUE)
#' @param scale.by Scale the size by size or radius.
#' @param scale.min Set lower limit for scaling, use NA for default values.
#' @param scale.max Set upper limit for scaling, use NA for default values.
#' @param plot.legend Plot the legends ?
#' @param do.return Return ggplot2 object ?
#' @param x.lab.pos Where to display x axis labels. This must be one of "bottom","top","both" or "none".
#' @param y.lab.pos Where to display y axis labels. This must be one of "left","right","both"or "none".
#' @param x.lab.size.factor Factor resizing x-axis labels (default=1)
#' @param y.lab.size.factor Factor resizing y-axis labels (default=1)
#' @param shape_var If numeric = Similar to pch : square=15; circle=16; triangle=17. Can also be a column/element name or a vector of the same size than the input dataset.
#' @param shape_use Shapes to uses (only when shape is controled by a discrete factor). Default shapes : \\u25A0 \\u25CF \\u25C6 \\u2BC8 \\u2BC7 \\u2BC6 \\u2BC5 \\u25D8 \\u25D9 \\u2726  \\u2605 \\u2736 \\u2737.
#' @param shape_legend Name of the shape legend if shape_var is a vector.
#' @param text.size Size of text to display on the shapes.
#' @param text.vjust Vertical justification of text to display on the shapes. Default value = 0, which mean no justification. Recommended value is between -0.5 and 0.5.
#' @param vertical_coloring Which color use to color the plot vertically ? (colors are repeated untill the end of the plot). Setting vertical and horizontal coloring at the same time is not recommended !
#' @param horizontal_coloring Which color use to color the plot horizontally ? (colors are repeated untill the end of the plot). Setting vertical and horizontal coloring at the same time is not recommended !
#' @param size.breaks.number Number of shapes with different size to display in the legend. Not used if size.breaks.values is not NA.
#' @param size.breaks.values Vector containing numerical labels for the size legend.
#' @param color.breaks.number Number of labels for the color gradient legend. Not used if color.breaks.values is not NA.
#' @param color.breaks.values Vector containing numerical labels for continuous color legend.
#' @param shape.breaks.number Number of shapes to display in the legend. Used when shape is controled by a continuous factor only. Not used if shape.breaks.values is not NA.
#' @param shape.breaks.values Vector containing numerical labels for continuous shape legend.
#' @param transpose Reverse x axis and y axis ?
#' @param dend_x_var A vector containing Column/List indexes or Column/List names to compute the x axis dendrogramm.
#' @param dend_y_var A vector containing Column/List indexes or Column/List names to compute the y axis dendrogramm.
#' @param dist_method The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param hclust_method The agglomeration method to be used. This must be one of "single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid" or "median".
#' @param do.plot Print the plot ? (default=TRUE)
#'
#' @import ggplot2
#' @importFrom grDevices colorRampPalette hcl
#' @importFrom stats as.dendrogram dist hclust na.omit
#' @importFrom FactoMineR PCA MCA FAMD
#' @importFrom scales rescale
#' @importFrom reshape2 dcast
#' @importFrom ggdendro segment dendro_data
#' @importFrom grImport2 readPicture symbolsGrob
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid textGrob grob gpar
#' @importFrom sisal dynTextGrob
#'
#' @return Print the plot (if do.plot=TRUE) and return a list containing input data, executed command, resulting dot plot and computed dendrograms (if do.return=TRUE)
#' @export
#' @examples
#' library(FlexDotPlot)
#' data(CBMC8K_example_data)
#' dotplot = dot_plot(data.to.plot=CBMC8K_example_data, size_var="RNA.avg.exp.scaled",
#' col_var="ADT.avg.exp.scaled", text_var="ADT.pct.exp.sup.cutoff",
#' shape_var="canonical_marker", shape_use = c("\u25CF","\u2737"),x.lab.pos="bottom",
#' y.lab.pos="left", cols.use=c("lightgrey","orange","red", "darkred"),
#' size_legend="RNA", col_legend="ADT", shape_legend="Canonical marker ?",
#' shape.scale =12, text.size=3, plot.legend = TRUE,
#' size.breaks.number=4, color.breaks.number=4, shape.breaks.number=5,
#' dend_x_var=c("RNA.avg.exp.scaled","ADT.avg.exp.scaled"),
#' dend_y_var=c("RNA.avg.exp.scaled","ADT.avg.exp.scaled"),
#' dist_method="euclidean",hclust_method="ward.D", do.return = TRUE)
#' @author Simon Leonard - simon.leonard@univ-rennes1.fr
plot_scpathway_contri_dot <- function(data.to.plot, size_var=NA,col_var=NA,
                                      text_var=NA, shape_var=16,
                                  size_legend="", col_legend="", shape_legend="",
                                  cols.use = "default",
                                  text.size=NA,  text.vjust=0,
                                  shape_use="default", shape.scale = 12,
                                  scale.by = "radius",
                                  scale.min = NA, scale.max = NA,
                                  plot.legend = TRUE,
                                  do.return = FALSE,
                                  x.lab.pos=c("both","top","bottom","none"),
                                  y.lab.pos=c("left","right","both","none"),
                                  x.lab.size.factor=1,
                                  y.lab.size.factor=1,
                                  vertical_coloring=NA,
                                  horizontal_coloring=NA,
                                  size.breaks.number=4,
                                  color.breaks.number=5,
                                  shape.breaks.number=5,
                                  size.breaks.values=NA, color.breaks.values=NA,
                                  shape.breaks.values=NA,
                                  display_max_sizes=TRUE,
                                  transpose=FALSE,
                                  dend_x_var=NULL, dend_y_var=NULL,
                                  dist_method=c("euclidean", "maximum", "manhattan", "canberra","binary", "minkowski"),
                                  hclust_method=c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2"),
                                  do.plot=TRUE
){

  ## For debuging
   #size_var=NA;col_var=NA; text_var=NA; shape_var=16;
   #size_legend=""; col_legend=""; shape_legend="";
   #cols.use = "default";
   #text.size=NA;  text.vjust=0;
   #shape_use="default"; shape.scale = 12;
   #scale.by = "radius"; scale.min = NA; scale.max = NA; plot.legend = TRUE; do.return = FALSE;
   #x.lab.rot = TRUE; x.lab.pos="both"; y.lab.pos="left";
   #x.lab.size.factor=1; y.lab.size.factor=1;
   #vertical_coloring=NA; horizontal_coloring=NA;
  #size.breaks.number=4; color.breaks.number=5; shape.breaks.number=5;
  #size.breaks.values=NA; color.breaks.values=NA; shape.breaks.values=NA;
  #display_max_sizes=TRUE;
  #transpose=FALSE;
  #dend_x_var=NULL; dend_y_var=NULL;
  #dist_method="euclidean";
  #hclust_method="ward.D2";
  #do.plot=TRUE

  x.lab.pos=match.arg(x.lab.pos)
  y.lab.pos=match.arg(y.lab.pos)


  # If cowplot library is loaded, we have to disable the cowplot default theme and set the ggplot2 default theme
  if("cowplot" %in% (.packages())){theme_set(theme_gray())}

  no_color_legend=FALSE # Is TRUE if col_var=NA => one legend to not plot
  no_size_legend=FALSE # Is TRUE if size_var=NA => one legend to not plot

  if (!(class(data.to.plot) %in% c("data.frame","list"))){stop("data.to.plot argument has to be a list or a data.frame")}
  if (all(is.na(c(size_var,col_var,text_var)))){
    if(((length(shape_var)!=nrow(data.to.plot)) )){
      if (length(shape_var)==1 & all(is.numeric(shape_var))){
        stop("Nothing to plot. Modify at least one of the following arguments : Shape size, shape color, text on shape, shape type")
      }
    }
  }


  ### 1: Data formatting ----
  ### 1.1: List to data.frame conversion ----
  if (inherits(data.to.plot, "list")){
    #Input=list -> Conversion to data.frame

    # Check - All the list elements are data.frame
    if (any(lapply(data.to.plot, class)!="data.frame")){stop("Provide a data.frame in each list element")}

    # Check - All the list elements have the same col names and row names (it means they all have the same dimension too)
    col_names=lapply(data.to.plot,colnames)
    if (all(sapply(col_names, identical, col_names[[1]]))==FALSE){stop("Provide data.frames with the same column names")}
    row_names=lapply(data.to.plot,row.names)
    if (all(sapply(row_names, identical, row_names[[1]]))==FALSE){stop("Provide data.frames with the same row names")}

    # Giving names to unnamed elements
    names(data.to.plot)=ifelse(names(data.to.plot)!="", names(data.to.plot), paste("Unnamed_column",1:length(data.to.plot), sep="_"))

    # d : Col 1 = row names; Col 2 = col names
    d=data.frame(row_names=rownames(data.to.plot[[1]])[row(data.to.plot[[1]])], col_names=colnames(data.to.plot[[1]])[col(data.to.plot[[1]])])
    # d2 = all the list elements converted in columns
    d2=do.call(data.frame,lapply(lapply(data.to.plot, c),unlist))

    data.to.plot=data.frame(d,d2)
    row.names(data.to.plot)=paste(data.to.plot$row_names,data.to.plot$col_names,sep="_")

  }

  ### 1.2: Determine/check parameters to plot ----
  #Input = dataframe
  if (ncol(data.to.plot)<3){
    stop("Provide a data.frame/list with at least 3 columns/elements")
  }

  save.data=data.to.plot
  data.to.plot=data.to.plot[,1:2]

  cat("Using : ")

  if(!length(size_var) %in% c(1,nrow(data.to.plot))){
    stop(paste("Length of size_var argument has to be equal to 1 or ",nrow(data.to.plot)," (the input dataset size)", sep=""))
  }
  if (length(size_var)==1){
    if (!is.na(size_var)){
      if (size_var %in% colnames(save.data)){
        if (is.numeric(save.data[,size_var])){
          cat(paste("\n -",size_var,"values to set shape size"))
          data.to.plot[,3]=save.data[,size_var]
          size_legend=ifelse(size_legend=="",size_var,size_legend)
        }else {stop(paste("size_var column (",size_var,") has to be numeric"))}
      } else if (is.numeric(size_var) & size_var %in% (1:ncol(save.data))){
        if (is.numeric(save.data[,size_var])){
          cat(paste("\n -",colnames(save.data)[size_var],"values to set shape size"))
          data.to.plot[,3]=save.data[,size_var]
          size_legend=ifelse(size_legend=="",colnames(save.data)[size_var], size_legend)
        }else {stop(paste("size_var column (",size_var,") has to be numeric"))}
      } else {stop("Size_var argument does not correspond to an element/column number or an element/column name")}
    }else {
      data.to.plot[,3]=1
      no_size_legend=TRUE
      cat(paste("\n - Nothing to set shape size"))
    }
  }else {
    if (is.numeric(size_var)){
      data.to.plot[,3]=size_var
      size_legend=size_legend
    }else{stop("Size var vector has to be numeric")}
  }

  if(!length(col_var) %in% c(1,nrow(data.to.plot))){
    stop(paste("Length of col_var argument has to be equal to 1 or ",nrow(data.to.plot)," (the input dataset size)", sep=""))
  }
  if (length(col_var)==1){
    if (!is.na(col_var)){
      if (col_var %in% colnames(save.data)){
        cat(paste("\n -",col_var,"values to set shape color"))
        data.to.plot[,4]=save.data[,col_var]
        col_legend=ifelse(col_legend=="",col_var,col_legend)
      } else if (is.numeric(col_var) & col_var %in% (1:ncol(save.data))){
        cat(paste("\n -",colnames(save.data)[col_var],"values to set shape color"))
        data.to.plot[,4]=save.data[,col_var]
        col_legend=ifelse(col_legend=="",colnames(save.data)[col_var], col_legend)
      } else {stop("Col_var argument does not correspond to an element/column number or an element/column name")}
    }else {
      data.to.plot[,4]="no_color"
      no_color_legend=TRUE
      cat(paste("\n - Nothing to set shape color"))
    }
  }else {
    data.to.plot[,4]=col_var
    col_legend=col_legend
  }

  if (length(cols.use)==1){
    if (cols.use=="default" & is.numeric(data.to.plot[,4])){
      cols.use=c("lightgrey", "blue")
    }
  }

  if(!length(text_var) %in% c(1,nrow(data.to.plot))){
    stop(paste("Length of text_var argument has to be equal to 1 or ",nrow(data.to.plot)," (the input dataset size)", sep=""))
  }
  if (length(text_var)==1){
    if (!is.na(text_var)){
      if (text_var %in% colnames(save.data)){
        cat(paste("\n -",text_var,"values to add text on shapes"))
        data.to.plot[,5]=save.data[,text_var]
      } else if (is.numeric(text_var) & text_var %in% (1:ncol(save.data))){
        cat(paste("\n -",colnames(save.data)[text_var],"values to add text on shapes"))
        data.to.plot[,5]=save.data[,text_var]
      } else {stop("Text_var argument does not correspond to an element/column number or an element/column name")}
    }else {
      data.to.plot[,5]=""
      cat(paste("\n - Nothing to add text on shapes"))
    }
  }else {
    data.to.plot[,5]=text_var
  }

  if(!length(shape_var) %in% c(1,nrow(data.to.plot))) {stop(paste("Length of shape_var argument has to be equal to 1 or ",nrow(data.to.plot)," (the input dataset size)", sep=""))}
  if(length(shape_var)==1){
    if (shape_var %in% colnames(save.data)){
      cat(paste("\n - ", shape_var," values to determine shapes (",
                length(unique(save.data[,shape_var]))," shapes detected)", sep=""))
      shape=save.data[,shape_var]
      shape_legend=ifelse(shape_legend=="",
                          shape_var,
                          shape_legend)
    }else if (is.numeric(shape_var)){
      shape=shape_var
    }else{stop("The shape_var argument is not numeric and does not correspond to a column name/list element name")
    }
  }else{shape=shape_var}

  data.to.plot[,1]=factor(data.to.plot[,1], levels=levels(as.factor(data.to.plot[,1])))
  data.to.plot[,2]=factor(data.to.plot[,2], levels=levels(as.factor(data.to.plot[,2])))
  data.to.plot[,5]=as.character(data.to.plot[,5])

  if (transpose){
    data.to.plot[,1:2]=data.to.plot[,2:1]
    save.data[,1:2]=save.data[,2:1]
  }



  ### 2 Calculate dendrograms ----
  check_FAMD_var=function(FAMD_var,type=c("x","y"), save.data){
    out=NULL
    if(!(is.vector(FAMD_var))){
      # Input is not a vector
      out$message=paste("FAMD_",type,"_var arguement has to be a vector",sep="")
      out$success=FALSE
    }else{
      if(is.vector(FAMD_var)){
        if(!is.numeric(FAMD_var)){
          if(all(FAMD_var %in% colnames(save.data))){
            FAMD_var=which(colnames(save.data) %in% FAMD_var)
          }else{
            out$message=paste("FAMD_",type,"_var names are not in element/column names",sep="")
            out$success=FALSE}
        }
        if(is.numeric(FAMD_var)){
          if(all(FAMD_var %in% 1:ncol(save.data))){
            if(2 %in% FAMD_var){}else{
              print(paste("In FAMD_",type,"_var : Adding y index",sep=""))
              FAMD_var=c(2,FAMD_var)
            }
            if(1 %in% FAMD_var){}else{
              print(paste("In FAMD_",type,"_var : Adding x index",sep=""))
              FAMD_var=c(1,FAMD_var)
            }

            out$famd_input=save.data[,FAMD_var]
            out$success=TRUE
          }else{
            out$message=paste("FAMD_",type,"_var indexes are not in element/column indexes",sep="")
            out$success=FALSE
          }
        }
      }
    }
    return(out)
  }

  format_FAMD_input=function(FAMD_input, type=c("x","y"), save.data){
    x_name=colnames(save.data)[1]
    y_name=colnames(save.data)[2]
    other_names=colnames(FAMD_input)[!colnames(FAMD_input) %in% c(x_name,y_name)]

    list=lapply(other_names, function(x){
      data=data.frame(dcast(FAMD_input, list(x_name,y_name), value.var=x))
      rownames(data)=data[,x_name]; data[,x_name]=NULL
      data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], as.factor)

      # return(data)
      if(type=="y"){
        return(data.frame(t(data)))
      }else{return(data)}

    })

    FAMD_finalinput=do.call(cbind,list)
    return(FAMD_finalinput)
  }

  run_FAMD_and_hclust=function(FAMD_final_input, metric, method){

    #Perform FAMD or PCA or MCA
    if( all(lapply(FAMD_final_input,class) %in% c("numeric","integer")) ){
      # Quantitate factors only -> PCA
      res.A=PCA(FAMD_final_input, graph=F, ncp=min(dim(FAMD_final_input))-1)
    } else if ( all(!lapply(FAMD_final_input,class) %in% c("numeric","integer")) ){
      # Qualitative factors only -> MCA
      res.A=MCA(FAMD_final_input, graph=F, ncp=min(dim(FAMD_final_input))-1)
    } else {
      # Mixed data -> FAMD
      res.A=FAMD(FAMD_final_input, graph=F, ncp=min(dim(FAMD_final_input))-1)
    }

    # Get coordinates, calculate distance and perform hierarchical clustering
    # Took from the HCPC function of FactoMineR package
    X = as.data.frame(res.A$ind$coord)

    do <- dist(X, method = metric)^2
    weight = rep(1, nrow(X))
    eff <- outer(weight, weight, FUN = function(x, y, n) {
      x * y/n/(x + y)
    }, n = sum(weight))
    dissi <- do * eff[lower.tri(eff)]

    hc <- hclust(dissi, method = method, members = weight)

    return(hc)
  }

  dist_method=match.arg(dist_method)
  hclust_method=match.arg(hclust_method)


  if(!is.null(dend_x_var)){
    # x dendrogramm

    # For the PCA/MCA/FAMD, we change the factor names in factor_1, factor2 ... to avoid syntax problem
    if(class(save.data[,1])!="factor"){save.data[,1]=as.factor(save.data[,1])}
    old_levels=levels(save.data[,1])
    save.data[,1]=paste("factor", as.numeric(save.data[,1]), sep="_")

    check_x=check_FAMD_var(dend_x_var, type="x",save.data=save.data)
    if (check_x$success){

      FAMD_x_input=format_FAMD_input(check_x$famd_input, type="x", save.data=save.data)
      hc_x_result=run_FAMD_and_hclust(FAMD_x_input, metric=dist_method, method=hclust_method)

    }else{
      #check_x$success==FALSE
      stop(check_x$message)
    }

    # Reorder x axis according to dendrogramms
    data.to.plot[,1]<- factor(data.to.plot[,1],
                              levels = old_levels[as.numeric(gsub("factor_","",hc_x_result$labels[hc_x_result$order]))])
    save.data[,1]=data.to.plot[,1]
  }else{hc_x_result=NULL}

  if(!is.null(dend_y_var)){
    #y dendrogram

    # For the PCA/MCA/FAMD, we change the factor names in factor_1, factor2 ... to avoid syntax problem
    if(class(save.data[,2])!="factor"){save.data[,2]=as.factor(save.data[,2])}
    old_levels=levels(save.data[,2])
    save.data[,2]=paste("factor",as.numeric(save.data[,2]),sep="_")

    check_y=check_FAMD_var(dend_y_var, type="y",save.data=save.data)
    if (check_y$success){

      FAMD_y_input=format_FAMD_input(check_y$famd_input, type="y", save.data=save.data)
      hc_y_result=run_FAMD_and_hclust(FAMD_y_input, metric=dist_method, method=hclust_method)

    }else{
      #check_y$success==FALSE
      stop(check_y$message)
    }

    # Reorder y axis according to dendrogramms
    data.to.plot[,2]<- factor(data.to.plot[,2],
                              levels = old_levels[as.numeric(gsub("factor_","",hc_y_result$labels[hc_y_result$order]))])
    save.data[,2]=data.to.plot[,2]

  }else{hc_y_result=NULL}



  ### 3: GGPlot object ----
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

  # Set theme function. Common function between classical dotplot and pacman plot
  # Transparent background, black border, no grid
  set_background=function(p, xlims, ylims, vertical_coloring, horizontal_coloring){
    p <- p + theme(
      panel.background = element_rect(fill="transparent",linetype="solid", color="black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA, fill = NA), axis.ticks = element_blank(),
      plot.margin = unit(c(0,0,0,0), "points")
    )
    p <- p + coord_cartesian(xlim=xlims,ylim=ylims,expand=FALSE, default=T)

    # Vertical coloring
    if(any(!is.na(vertical_coloring))){

      shading <- data.frame(min = c(0.5,seq(from = 1.5, to = max(as.numeric(as.factor(data.to.plot[,1]))), by = 1)),
                            max = c(seq(from = 1.5, to = max(as.numeric(as.factor(data.to.plot[,1]))) + 0.5, by = 1) ))
      # shading$col = rep_len(x=c(NA,"gray80"),length.out=length(unique(data.to.plot[,1]))))
      shading$col=rep(vertical_coloring, length.out=nrow(shading))

      p <- p + annotate(geom = "rect", xmin = shading$min, xmax = shading$max,
                        ymin = 0.5, ymax = max(as.numeric(data.to.plot[,2]))+0.5, fill = shading$col, col="black")

      # Adding black lines to delimitate shades
      # We use another annotate to not display first and last line; otherwise add colour="black" in the previous annotate
      p <- p + annotate(geom = "segment", x = c(shading$min[-1], shading$max[-length(shading$max)]),
                        xend = c(shading$min[-1], shading$max[-length(shading$max)]),
                        y = 0.5, yend = max(as.numeric(data.to.plot[,2]))+0.5,
                        colour = "black")
    }

    # Horizontal coloring
    if(any(!is.na(horizontal_coloring))){
      shading <- data.frame(min = seq(from = 0.5, to = max(as.numeric(as.factor(data.to.plot[,2]))), by = 1),
                            max = seq(from = 1.5, to = max(as.numeric(as.factor(data.to.plot[,2]))) + 0.5, by = 1))
      # shading$col = rep_len(x=c(NA,"gray80"),length.out=length(unique(data.to.plot[,2])))
      shading$col=rep(horizontal_coloring, length.out=nrow(shading))

      p <- p + annotate(geom = "rect", xmin=0.5,xmax=max(as.numeric(data.to.plot[,1]))+0.5, ymin = shading$min, ymax = shading$max,
                        fill = shading$col, col="black")

      # Adding black lines to delimitate shades
      # We use another annotate to not display first and last line; otherwise add colour="black" in the previous annotate
      p <- p + annotate(geom = "segment", x=0.5,xend=max(as.numeric(data.to.plot[,1]))+0.5, y = c(shading$min[-1], shading$max[-length(shading$max)]),
                        yend = c(shading$min[-1], shading$max[-length(shading$max)]),
                        colour = "black")
    }

    return(p)
  }

  get_shape_colors=function(data.to.plot, cols.use="default", color.breaks.values, color.breaks.number){
    # Functions wich define a color of each data.to.plot row according to the 4th column
    # Return a list containing :
    #   $palette : used color palette
    #   $colors : color of each row
    #   $labels : labels for legend (continuous color only)
    #   $breaks : legend breaks (continuous color only)


    shape_colors_labels=list()

    if (is.numeric(data.to.plot[,4])){
      if (length(x = cols.use) == 1) {
        shape_colors_labels$colors=rep(cols.use, nrow(data.to.plot)) # Monochrome
      }

      if (!(all(is.na(color.breaks.values)))){
        if(all(is.numeric(color.breaks.values))){
          cat("\n Putting color.breaks.values in color legend")
          shape_colors_labels$breaks=color.breaks.values
        }else{
          cat("\n Non numeric value in color.breaks.values, considering color.breaks.number instead")
          shape_colors_labels$breaks=(seq(min(na.omit(data.to.plot[,4])),max(na.omit(data.to.plot[,4])), length.out=color.breaks.number))
        }
      }else{
        shape_colors_labels$breaks=(seq(min(na.omit(data.to.plot[,4])),max(na.omit(data.to.plot[,4])), length.out=color.breaks.number))
      }

      color_breaks=shape_colors_labels$breaks
      shape_colors_labels$labels=ifelse((abs(color_breaks)<1e-2 | abs(color_breaks)>1e2) & color_breaks!=0,
                                        format(color_breaks, scientific=TRUE, digits=3),
                                        round(color_breaks,2)) # Values <1e-3 or >1e3 (excepted 0), are displayed with scientific notation

      map2color<-function(x,pal,limits=NULL){
        if(is.null(limits)) limits=range(x, na.rm = T)
        pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
      }

      shape_colors_labels$colors=map2color(data.to.plot[,4], pal=colorRampPalette(cols.use)(100))

      # if (length(x = cols.use) == 2){
      #   p <- p + scale_fill_gradient(low = cols.use[1], high = cols.use[2],
      #                                breaks=color_breaks, labels=color_labels) # Gradient from the first color to the second
      # } else {
      #   p <- p + scale_fill_gradientn(colours=cols.use, breaks=color_breaks, labels=color_labels)
      # }


    } else {
      # discrete colors
      if(all(cols.use !="default")){
        if (length(cols.use)>length(unique(data.to.plot[,4]))){
          cols.use=cols.use[1:length(unique(data.to.plot[,4]))]
          cat(paste("\n To much colors are supplied. Only the first",length(unique(data.to.plot[,4])),"are used"))
        } else if (length(cols.use)<length(unique(data.to.plot[,4]))){
          cols.use=rep_len(cols.use, length.out = length(unique(data.to.plot[,4])))
          cat(paste("\n The number of colors is lower than the modality number. Re-using colors"))
        }

        # p <- p + scale_fill_manual(values = cols.use)
      }
      else{
        #reproducing ggplot2 default discrete palette
        gg_color_hue <- function(n) {
          hues = seq(15, 375, length = n + 1)
          hcl(h = hues, l = 65, c = 100)[1:n]
        }

        cols.use=gg_color_hue(length(unique(data.to.plot[,4])))

      }

      shape_colors_labels$colors=cols.use[as.factor(data.to.plot[,4])]
      shape_colors_labels$labels=levels(as.factor(data.to.plot[,4]))

    }

    shape_colors_labels$palette=cols.use

    return(shape_colors_labels)
  }

  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  #Need to specify xlims and ylims of plot - this lets us calculate the normalised plot coordinates
  xlims <- c(0.5,length(unique(data.to.plot[,1]))+0.5)
  ylims <- c(0.5,length(unique(data.to.plot[,2]))+0.5)

  ### 3.1: PACMAN PLOT ----
  ### 3.2: DOT PLOT ----

    # ggplot object initilization
    p <- ggplot(mapping = aes(x = as.numeric(data.to.plot[,1]), y = as.numeric(data.to.plot[,2]), label=data.to.plot[,5]))

    # theme function
    p <- set_background(p, xlims, ylims, vertical_coloring, horizontal_coloring)

      # Displaying shapes
     # if(length(shape)==1){
        p <- p + geom_point(mapping = aes(size = data.to.plot[,3], color = data.to.plot[,4]), shape=shape)

        # Displaying maximal shapes area (only when shape %in% c(15,16,17,18))
        if (all(shape %in% c(15,16,17,18)) & display_max_sizes){
          background_shape=ifelse(shape %in% c(15,16,17),shape-15, shape-13)
          p <- p + geom_point(mapping = aes(size = max(na.omit(data.to.plot[,3]))), colour="black",shape=background_shape)
        }

      # Changing shape size + shape size legend
      if(!no_size_legend){
        if (!(all(is.na(size.breaks.values)))){
          if(all(is.numeric(size.breaks.values))){
            cat("Putting size.breaks.values in size legend")
            legend_breaks=size.breaks.values
          }else{
            cat("Non numeric value in size.breaks.values, considering size.breaks.number instead")
            legend_breaks=(seq(min(na.omit(data.to.plot[,3])),max(na.omit(data.to.plot[,3])), length.out=size.breaks.number))
          }
        }else{
          legend_breaks=(seq(min(na.omit(data.to.plot[,3])),max(na.omit(data.to.plot[,3])), length.out=size.breaks.number))
        }


        legend_labels=ifelse((abs(legend_breaks)<1e-2 | abs(legend_breaks)>1e2) & legend_breaks!=0,
                             format(legend_breaks, scientific=TRUE, digits=3),
                             round(legend_breaks,2)) # Values <1e-3 or >1e3 (excepted 0), are displayed with scientific notation
        p <- p + scale.func(range = c(0.1, shape.scale),limits = c(scale.min, scale.max), breaks = legend_breaks, labels = legend_labels)
        # /!\ It is important to no set the minimal range value to 0 (and >0.1) because its can induce shape size and position errors with specific shapes

      }


      # Coloring shapes
      # Assign color to each point
      plot_colors=get_shape_colors(data.to.plot, cols.use, color.breaks.values, color.breaks.number)
      # Display colors
      if (is.numeric(data.to.plot[,4])){
        p <- p + scale_color_gradientn(colors=plot_colors$palette, breaks=plot_colors$breaks, labels=plot_colors$labels, na.value = "transparent")
      }else{
        p <- p + scale_color_manual(values = plot_colors$palette, na.value="transparent")
      }


    # Displaying text
    p <- p + geom_text(aes(y = as.numeric(data.to.plot[,2]) + text.vjust), size=text.size)


  ### 3.3: Command lignes in common ----

  # Legend titles
  p$labels$size=size_legend
  p$labels$colour=col_legend
  p$labels$shape=shape_legend

  # Deleting axis labels (printed in an other grob)
  p <- p + scale_y_continuous(breaks = NULL,labels = NULL) + scale_x_continuous(breaks = NULL,labels = NULL)

  # Deleting axis titles
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

  # Remove panel border wich contain plot and legend
  p <- p + theme(panel.background = element_rect(fill="transparent",linetype=0))

  # Add panel border which contain plot only (legend outside) => /!\ It adds an additional margin in dotplot (not in pacman plot)
  p <- p + geom_rect(aes(xmin=0.5,xmax=max(as.numeric(data.to.plot[,1]))+0.5,ymin=0.5,ymax=max(as.numeric(data.to.plot[,2]))+0.5), alpha=0, colour="black")

  # Remove legend (printed in another grob)
  dot_plot_legend=get_legend(p)
  p <- p + theme(legend.position = "none")

  ### 4: Arrange plot with labels, dendrogramms and legends ----
  final.plot.list=list()

  # Plot number description :
  # 1 : Horizontal dendrogram
  # 2 : Top x labels
  # 3 : Vertical dendrogram
  # 4 : Left y labels
  # 5 : Dotplot
  # 6 : Right y labels
  # 7 : Size legend
  # 8 : Shape legend
  # 9 : Color legend
  # 10 : Bottom x labels

  layout=rbind(c(NA,NA,1,NA,NA,NA,NA),
               c(NA,NA,2,NA,NA,NA,NA),
               c(3,4,5,6,7,8,9),
               c(NA,NA,10,NA,NA,NA,NA))

  widths=c(1,3,10,3,3,3,3)
  x_lab_heights=3
  heigths=c(1,x_lab_heights,10,x_lab_heights)

  remove_geom <- function(ggplot2_object, geom_type) {
    # Delete layers that match the requested type.
    layers <- lapply(ggplot2_object$layers, function(x) {
      if (class(x$geom)[1] %in% geom_type) {
        NULL
      } else {
        x
      }
    })
    # Delete the unwanted layers.
    layers <- layers[!sapply(layers, is.null)]
    ggplot2_object$layers <- layers
    ggplot2_object
  }

  p_raw=remove_geom(p, c("GeomCustomAnn","GeomPoint","GeomRect","GeomText","GeomSegment"))


  ### 4.1 Axis Labels ----

  # x_coords=rescale(x = 1:length(levels(data.to.plot[,1])), to = c(0,1), from=c(0.5,length(levels(data.to.plot[,1]))+0.5))
  # final.plot.list[[2]]=dynTextGrob(levels(data.to.plot[,1]), x=x_coords, rot=ifelse(x.lab.rot, 90, 0),just="bottom", y=0.05, width=ifelse(x.lab.rot, 1, 1/length(levels(data.to.plot[,1]))* x.lab.size.factor))
  # final.plot.list[[10]]=dynTextGrob(levels(data.to.plot[,1]), x=x_coords, rot=ifelse(x.lab.rot, 90, 0),just="top", y=0.95, width=ifelse(x.lab.rot, 1, 1/length(levels(data.to.plot[,1]))* x.lab.size.factor))

  x_label_table=data.frame(xtext=levels(data.to.plot[,1]))
  final.plot.list[[2]]= p_raw + geom_text(data = x_label_table, mapping = aes_(label=~xtext), y=0.05, x=1:length(levels(data.to.plot[,1])), hjust=0, vjust=0.5, angle=90, size=3.88*x.lab.size.factor)
  final.plot.list[[2]]= final.plot.list[[2]]+
    coord_cartesian(ylim = c(0,1),xlim=c(0.5,length(unique(data.to.plot[,1]))+0.5),expand=F, default = T) +
    theme(plot.margin = unit(c(0,2,0,0), units = "points"), plot.background = element_rect(fill='transparent', color=NA))

  final.plot.list[[10]]= p_raw + geom_text(data = x_label_table, mapping = aes_(label=~xtext), y=0.95, x=1:length(levels(data.to.plot[,1])), hjust=1, vjust=0.5, angle=90, size=3.88*x.lab.size.factor)
  final.plot.list[[10]]= final.plot.list[[10]]+
    coord_cartesian(ylim = c(0,1),xlim=c(0.5,length(unique(data.to.plot[,1]))+0.5),expand=F, default = T) +
    theme(plot.margin = unit(c(2,0,0,0), units = "points"),plot.background = element_rect(fill='transparent', color=NA))

  if(x.lab.pos == "top"){
    heigths[4]=0
    final.plot.list[[10]]=grob()
  }else if (x.lab.pos=="bottom"){
    heigths[2]=0
    final.plot.list[[2]]=grob()
  }else if (x.lab.pos == "none"){
    heigths[c(2,4)]=0
    final.plot.list[[2]]=grob()
    final.plot.list[[10]]=grob()
  }

  # y_coords=rescale(x = seq(1, length(levels(data.to.plot[,2]))), to = c(0,1), from=c(0.5,length(levels(data.to.plot[,2]))+0.5))
  # final.plot.list[[4]]=dynTextGrob(levels(data.to.plot[,2]), x=0.95, y=y_coords, width=0.95,just="right", adjustJust = F)
  # final.plot.list[[4]]=textGrob(levels(data.to.plot[,2]),y=y_coords, x=0.95,gp = gpar(fontsize = 10), just="right")
  # final.plot.list[[6]]=dynTextGrob(levels(data.to.plot[,2]), x=0.05, y=y_coords,  width=0.95,just="left")

  y_label_table=data.frame(ytext=levels(data.to.plot[,2]))
  final.plot.list[[4]]= p_raw + geom_text(data = y_label_table, mapping = aes_(label=~ytext), x=1, y=1:length(levels(data.to.plot[,2])), hjust=1, vjust=0.5, size=3.88*y.lab.size.factor)
  final.plot.list[[4]]= final.plot.list[[4]]+coord_cartesian(xlim = c(0,1),ylim=c(0.5,length(unique(data.to.plot[,2]))+0.5),expand=F, default = T) +
    theme(plot.margin = unit(c(0,2,2,0), units = "points"),plot.background = element_rect(fill='transparent', color=NA))

  final.plot.list[[6]]= p_raw + geom_text(data = y_label_table, mapping = aes_(label=~ytext), x=0.05, y=1:length(levels(data.to.plot[,2])), hjust=0, vjust=0.5, size=3.88*y.lab.size.factor)
  final.plot.list[[6]]= final.plot.list[[6]]+coord_cartesian(xlim = c(0,1),ylim=c(0.5,length(unique(data.to.plot[,2]))+0.5),expand=F, default = T) +
    theme(plot.margin = unit(c(0,0,2,2), units = "points"), plot.background = element_rect(fill='transparent', color=NA))

  if(y.lab.pos == "left"){
    widths[4]=0
    final.plot.list[[6]]=grob()
  }else if (y.lab.pos=="right"){
    widths[2]=0
    final.plot.list[[4]]=grob()
  }else if (y.lab.pos == "none"){
    widths[c(2,4)]=0
    final.plot.list[[4]]=grob()
    final.plot.list[[6]]=grob()
  }



  ### 4.2 Dendrogramms ----
  # Arrange dotplot with dendrogramms

  if(!is.null(hc_x_result)){
    # Arrange x dendrogram
    ddata_x <- segment(dendro_data(hc_x_result))

    dendro_horizontal <- p_raw + geom_segment(data = ddata_x, mapping = aes_(x=~x,xend=~xend,
                                                                             y=~y,yend=~yend, label=NULL))

    dendro_horizontal <- dendro_horizontal+coord_cartesian(ylim = c(min(c(ddata_x$y,ddata_x$yend)),
                                                                    max(c(ddata_x$y,ddata_x$yend))),
                                                           xlim=c(0.5,length(unique(data.to.plot[,1]))+0.5),
                                                           expand=F, default = T)+ theme(plot.margin = unit(c(2,0,2,0), units = "points"))
  }

  if(!is.null(hc_y_result)){
    # Arrange y dendrogram
    ddata_y <- segment(dendro_data(hc_y_result))

    dendro_vertical <- p_raw + geom_segment(data = ddata_y, mapping = aes_(x=~length(unique(data.to.plot[,2]))+0.5-y,
                                                                           xend=~length(unique(data.to.plot[,2]))+0.5-yend,
                                                                           y=~x,yend=~xend, label=NULL))

    dendro_vertical <- dendro_vertical+coord_cartesian(xlim = range(c(length(unique(data.to.plot[,2]))+0.5-ddata_y$y,length(unique(data.to.plot[,2]))+0.5-ddata_y$yend)),
                                                       ylim=c(0.5,length(unique(data.to.plot[,2]))+0.5),
                                                       expand=F, default = T) + theme(plot.margin = unit(c(0,0,0,2), units = "points"))
  }


  if(!is.null(hc_x_result) & !is.null(hc_y_result)){

    final.plot.list[[1]]=dendro_horizontal
    final.plot.list[[3]]=dendro_vertical
    final.plot.list[[5]]=p


  } else if (!is.null(hc_x_result)){

    final.plot.list[[1]]=dendro_horizontal
    final.plot.list[[3]]=grob()
    final.plot.list[[5]]=p

    widths[1]=0

  } else if (!is.null(hc_y_result)){

    final.plot.list[[1]]=grob()
    final.plot.list[[3]]=dendro_vertical
    final.plot.list[[5]]=p

    heigths[1]=0

  } else {

    final.plot.list[[1]]=grob()
    final.plot.list[[3]]=grob()
    final.plot.list[[5]]=p

    heigths[1]=0
    widths[1]=0

  }


  ### 4.3 Legends ----
  # 7 : Size legend
  # 8 : Shape legend
  # 9 : Color legend

  if (plot.legend){

    if(is.numeric(shape) & length(shape)==nrow(data.to.plot) & plot.legend){
      #pacman plot

      # size legend
      if(!no_size_legend){
        final.plot.list[[7]]=size_leg
      }else{
        final.plot.list[[7]]=grob()
        widths[5]=0
      }

      # shape legend
      if(length(shape)==nrow(data.to.plot)){
        final.plot.list[[8]]=shape_leg
      }else{
        final.plot.list[[8]]=grob()
        widths[6]=0
      }

      # color legend
      if(!no_color_legend){
        final.plot.list[[9]]=col_leg
      }else{
        final.plot.list[[9]]=grob()
        widths[7]=0
      }

    } else {
      final.plot.list[[7]]=dot_plot_legend
      final.plot.list[[8]]=grob()
      final.plot.list[[9]]=grob()
      widths[6:7]=0
    }



  }else{
    final.plot.list[[7]]=grob()
    final.plot.list[[8]]=grob()
    final.plot.list[[9]]=grob()
    widths[5:7]=0
  }


  ### 4.4 Render plot ----
  # return(final.plot.list) # for debug

  final_plot=arrangeGrob(grobs = final.plot.list, layout_matrix = layout, widths = widths, heights = heigths)

  if(do.plot){
    grid.arrange(final_plot)
  }
  if (do.return) {
    final_output=list()
    final_output[["input_data"]]=save.data
    final_output[["command"]]=sys.calls()[[1]]
    final_output[["plot"]]=final_plot
    if(!is.null(hc_y_result)){final_output[["dendrogram_y"]]=as.dendrogram(hc_y_result)}else{final_output[["dendrogram_y"]]=NULL}
    if(!is.null(hc_x_result)){final_output[["dendrogram_x"]]=as.dendrogram(hc_x_result)}else{final_output[["dendrogram_x"]]=NULL}
    if(!is.null(hc_y_result) | !is.null(hc_x_result)){final_output[["raw_dend_ggplot"]]=p_raw}

    return(final_output)
  }

  ### End function ----
}
