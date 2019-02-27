# Visualize cms and groups


#' visHist
#'
#' Plot pvalue histograms of cms score distributions
#'
#' @param cms_res cms result matrix with one or two columns containing cms scores and/or smoothened cms scores.
#'
#' @details Plots cms score distribution similar to a pvalue histogram distribution.
#' Smoothened cms should be close to a normal distribution centered around 0.5 and cms scores should be approx. flat, if no dataset-specific bias are expected.
#'
#' @family visualize cms functions
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom ggplot2 ggplot aes_string geom_histogram ggtitle xlab theme_classic
#' @importFrom cowplot plot_grid
visHist <- function(cms_res){
  p <- do.call(plot_grid, c(lapply(colnames(cms_res), function(cms_name, ...){
    ggplot(as.data.frame(cms_res), aes_string(x=cms_name)) +
            geom_histogram(color="black", fill = col_hist[which(colnames(cms_res) %in% cms_name)]) +
            ggtitle(cms_name) + xlab(cms_name) + theme_classic()
  }), ncol = ifelse(length(colnames(cms_res)) > 1, 2, 1)))
  p
}


#' visOverview
#'
#' Plot cms, smoothened cms, group label and other Variable defined in colData in a reduced dimensional representation.
#'
#' @param cms_res cms result matrix with one or two columns containing cms scores and/or smoothened cms scores.
#' Rownames need to correspond to colnames of sce.
#' @param sce Combined \code{\link{SingleCellExperiment}} object containing either precalculated reduced dimension embeddings or counts or logcounts of all groups merged.
#' If precalculated dimension reduction embeddings are used, they need to be within the reducedDimensions slot.
#' @param group Character string specifying the name of variable used to define groups (batches). Should be have a corresonding element in colData(sce).
#' @param dim_red Character defining dimension reduction embeddings to plot. Default is TSNE as from \code{\link{runTSNE}}.
#' If none is provided and no dimension reduction slot named "TSNE" is within \code{\link{reducedDimNames}} tsne embeddings will be calculated by \code{\link{runTSNE}}.
#' @param smooth A logical value indicating if smoothened cms scores should be plotted.
#' If true cms_res needs to contain 2 colums with cms and cms_smoothened results. Default is set by the length of colnames(cms_res)
#' @param log10_val A logical indicating if -log10(cms) should be plotted to visualize differences of small values
#' @param other_Var Character string defining other variables to be plotted asided (e.g. some specified celltype label).
#' Need correspond to one of colData(sce).
#'
#' @details Plots reduced dimensions (tsne as default) of cells using the group variable and the cms score as colors.
#' Other color label as celltype label or smoothened cms scores can be plotted aside. Generates tsne embeddings, if none have been specified.
#' Embeddings from data integration methods (e.g. mnn.correct) can also been used as long as they are specified in \code{\link{reducedDimNames}}.
#'
#' @family visualize cms functions
#' @seealso \code{\link{visCms}}, \code{\link{visGroup}}
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom ggplot2 ggplot aes_string ylab xlab theme_void theme guide_legend guides element_blank element_line geom_point scale_color_manual
#' @importFrom cowplot plot_grid
#' @importFrom scater runTSNE
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDimNames reducedDim colData
#' @importFrom viridis scale_color_viridis
visOverview <- function(cms_res, sce, group, dim_red = "TSNE", smooth = ifelse(length(colnames(cms_res)) > 1, TRUE, FALSE), log10_val = FALSE, other_Var = NULL){
  cell_names <- colnames(sce)
  #Compare order sce and cms_restheme_classic
  stopifnot(rownames(cms_res) == cell_names)
  #custumized dimreduction
  if(!dim_red %in% "TSNE"){
    if(!dim_red %in% reducedDimNames(sce)){
      stop("Ambigous parameter 'dim_red', provide one of:
           * A dim_red method that is listed in reducedDimNames(sce).
           * Default('TSNE') will call runTSNE to calculate a subspace.")
    }
    red_dim <- as.data.frame(reducedDim(sce, dim_red))
  }else{
    #used tsne from scater package (check for availability first)
    if(is.null(reducedDim(sce, "TSNE"))){
      #use "logcounts" if availabe otherwise "counts"
      if(names(assays(sce)) %in% "logcounts"){
        sce <- runTSNE(sce)
      }else{
        sce <- runTSNE(sce, exprs_values = "counts")
      }
    }
    red_dim <- as.data.frame(reducedDim(sce, "TSNE"))
  }
    colnames(red_dim) <- c("red_dim1", "red_dim2") # ensure  the same colnames
    #data frame to plot
    df <- data.frame(sample_id = cell_names, group_var = as.factor(colData(sce)[, group]),
                     red_Dim1=red_dim$red_dim1, red_Dim2=red_dim$red_dim2, cms = cms_res[,"cms"])
    # use of -log10 values (for very low cms)
    if(log10_val == TRUE){
      df$cms <- -log10(df$cms)
    }

    t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
      xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
          theme_void() +
          theme(aspect.ratio=1, panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="grey", size=.3))


    t_group <- t + geom_point(size=1, alpha = 0.3, aes_string(color="group_var")) +
      guides(color=guide_legend(override.aes=list(size=1))) +
      scale_color_manual(values = col_group) +
      ggtitle(group)

    t_cms <- t + geom_point(size=1, alpha = 0.5, aes_string(color="cms")) +
      guides(color=guide_legend(override.aes=list(size=0.5))) +
      scale_color_viridis(option = "B") + ggtitle(paste0("cms : ", group))

    p <- plot_grid(t_group, t_cms)

    if(smooth == TRUE){
      df$cms_smooth <- cms_res[,"cms_smooth"]

      t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
        xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
        theme_void() +
        theme(aspect.ratio=1, panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color="grey", size=.3))

      t_smooth <- t + geom_point(size=1, alpha = 0.3, aes_string(color="cms_smooth")) +
        guides(color=guide_legend(override.aes=list(size=0.5))) +
        scale_color_viridis(option = "B") + ggtitle(paste0("smooth cms : ", group))
      p <- plot_grid(t_group, t_cms, t_smooth)

    }

    if(!is.null(other_Var)){
      df[,other_Var] <- colData(sce)[, other_Var]

      t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
        xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
        theme_void() +
        theme(aspect.ratio=1, panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color="grey", size=.3))

      t_other <- t + geom_point(size=1, alpha = 0.3, aes_string(color=other_Var)) +
        guides(color=guide_legend(override.aes=list(size=2))) + ggtitle(other_Var)

      if(is.numeric(df[,other_Var])){
        t_other <- t_other + scale_color_viridis(option = "C")
      }else{
        t_other <- t_other + scale_color_manual(values = col_hist)
      }

      if(smooth == TRUE){
      p <- plot_grid(t_group, t_cms, t_smooth, t_other)
      }else{
        p <- plot_grid(t_group, t_cms, t_other)
      }

    }

    #plot
    return(p)

}


### Single plots


#' visCms
#'
#' Plot cms scores in a reduced dimensional plot like tsne.
#'
#' @param cms_res cms result matrix with one or two columns containing cms scores and/or smoothened cms scores.
#' Rownames need to correspond to colnames of sce. Colnames are allowed to be set different, but need to be specified in cms_var.
#' @param sce Combined \code{\link{SingleCellExperiment}} object containing either precalculated reduced dimension embeddings or counts or logcounts of all groups merged.
#' If precalculated dimension reduction embeddings are used, they need to be within the reducedDimensions slot.
#' @param cms_var character string specifying the cms scores to use (usually "cms_smoothened" or "cms").
#' Needs to correspond to one of the colnames of cms_res. Default is "cms".
#' @param dim_red Character defining dimension reduction embeddings to plot. Default is TSNE as from \code{\link{runTSNE}}.
#' If none is provided and no dimension reduction slot named "TSNE" is within \code{\link{reducedDimNames}} tsne embeddings will be calculated by \code{\link{runTSNE}}.
#' @param log10_val A logical indicating if -log10(cms) should be plotted to visualize differences of small values
#' @param ... Additional arguments to pass to \code{\link{ggplot}}
#'
#' @details Plots a reduced dimension plot (tsne as default) colored by cms scores.
#' The dimesion reduction embedding can be specified, but only tsne embeddings will automatically be computed.
#' Embeddings from data integration methods (e.g. mnn.correct) can also been used as long as they are specified in \code{\link{reducedDimNames}} of sce.
#'
#' @seealso \code{\link{visOverview}}, \code{\link{visGroup}}
#' @family visualize cms functions
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom ggplot2 ggplot aes_string ylab xlab theme_void theme guide_legend guides element_blank element_line geom_point
#' @importFrom scater runTSNE
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDimNames reducedDim colData
#' @importFrom viridis scale_color_viridis
visCms <- function(cms_res, sce, cms_var = "cms", dim_red = "TSNE", log10_val = FALSE, ...){
  #Generate or specify dim reduction
  cell_names <- colnames(sce)
  #Compare order sce and cms_res
  stopifnot(rownames(cms_res) == cell_names)
  #custumized dimreduction
  if(!dim_red %in% "TSNE"){
    if(!dim_red %in% reducedDimNames(sce)){
      stop("Please provide a dim_red method to be used that is listed in the reducedDim slot of sce")
    }
    red_dim <- as.data.frame(reducedDim(sce, dim_red))
  }else{
    #used tsne from scater package (check for availability first)
    if(is.null(reducedDim(sce, "TSNE"))){
      #use "logcounts" if availabe otherwise "counts"
      if(names(assays(sce)) %in% "logcounts"){
        sce <- runTSNE(sce)
      }else{
        sce <- runTSNE(sce, exprs_values = "counts")
      }
    }
    red_dim <- as.data.frame(reducedDim(sce, "TSNE"))
  }
  colnames(red_dim) <- c("red_dim1", "red_dim2") # ensure  the same colnames
  #data frame to plot
  df <- data.frame(sample_id = cell_names, cms = cms_res[,cms_var],
                   red_Dim1=red_dim$red_dim1, red_Dim2=red_dim$red_dim2)
  # use of -log10 values (for very low cms)
  if(log10_val == TRUE){
    df$cms <- -log10(df$cms)
  }

  t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
    xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
    theme_void() +
    theme(aspect.ratio=1, panel.grid.minor=element_blank(),
          panel.grid.major=element_line(color="grey", size=.3))


  t_cms <- t + geom_point(size=1, alpha = 0.5, aes_string(color="cms")) +
    guides(color=guide_legend(override.aes=list(size=2), title = cms_var)) +
    scale_color_viridis(option = "B") + ggtitle(cms_var)

  t_cms
}


#' visGroup
#'
#' Plot group label in a reduced dimensional plot like tsne or integrated embeddings.
#'
#' @param sce Combined \code{\link{SingleCellExperiment}} object containing either precalculated reduced dimension embeddings or counts or logcounts of all groups merged.
#' If precalculated dimension reduction embeddings are used, they need to be within the reducedDimensions slot.
#' @param group character string specifying the group variable to label cells. Should be included into the colData slot of sce.
#' @param dim_red Character defining dimension reduction embeddings to plot. Default is TSNE as from \code{\link{runTSNE}}.
#' If none is provided and no dimension reduction slot named "TSNE" is within \code{\link{reducedDimNames}} tsne embeddings will be calculated by \code{\link{runTSNE}}.
#' @param ... Additional arguments to pass to \code{\link{ggplot}}
#'
#' @details Plots a reduced dimension plot (tsne as default) colored by group parameter.
#' The dimesion reduction embedding can be specified, but only tsne embeddings will automatically be computed.
#' Embeddings from data integration methods (e.g. mnn.correct) can also been used as long as they are specified in \code{\link{reducedDimNames}} of sce.
#'
#' @seealso \code{\link{visOverview}}, \code{\link{visCms}}
#' @family visualize cms functions
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom ggplot2 ggplot aes_string ylab xlab theme_void theme guide_legend guides element_blank element_line geom_point scale_color_manual
#' @importFrom scater runTSNE
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDimNames reducedDim colData
#' @importFrom viridis scale_color_viridis
visGroup <- function(sce, group, dim_red = "TSNE", ...){
  #Generate or specify dim reduction
  cell_names <- colnames(sce)
  #custumized dimreduction
  if(!dim_red %in% "TSNE"){
    if(!dim_red %in% reducedDimNames(sce)){
      stop("Please provide a dim_red method to be used that is listed in the reducedDim slot of sce")
    }
    red_dim <- as.data.frame(reducedDim(sce, dim_red))
  }else{
    #used tsne from scater package (check for availability first)
    if(is.null(reducedDim(sce, "TSNE"))){
      #use "logcounts" if availabe otherwise "counts"
      if(names(assays(sce)) %in% "logcounts"){
        sce <- runTSNE(sce)
      }else{
        sce <- runTSNE(sce, exprs_values = "counts")
      }
    }
    red_dim <- as.data.frame(reducedDim(sce, "TSNE"))
  }
  colnames(red_dim) <- c("red_dim1", "red_dim2") # ensure  the same colnames
  #data frame to plot
  df <- data.frame(sample_id = cell_names, group_var = colData(sce)[, group],
                   red_Dim1=red_dim$red_dim1, red_Dim2=red_dim$red_dim2)


  t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
    xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
    theme_void() +
    theme(aspect.ratio=1, panel.grid.minor=element_blank(),
          panel.grid.major=element_line(color="grey", size=.3))

  t_group <- t + geom_point(size=1, alpha = 0.5, aes_string(color="group_var")) +
    guides(color=guide_legend(override.aes=list(size=1), title = group)) +
    ggtitle(group)

  # fit colors to continous or discrete scale
  if(is.numeric(df$group_var)){
    t_group <- t_group + scale_color_viridis(option = "C")
  }else{
    t_group <- t_group + scale_color_manual(values = col_group)
  }

  t_group
}



