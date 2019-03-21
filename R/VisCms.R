# Visualize cms and groups


#' visHist
#'
#' Plot pvalue histograms of cms score distributions
#'
#' @param cms_res Matrix or data.frame. Cms scores to plot should be in columns and cells in rows.
#' @param ncol Numeric. Number of columns of the pval histogram.
#'
#'
#' @details Plots cms score distribution similar to a pvalue histogram distribution.
#' Without dataset-specific bias, cms scores should be approx. flat distributed.
#'
#' @family visualize cms functions
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' load(system.file("extdata/cms_sim30.rda", package = "CellMixS"))
#' visHist(cms_sim30)
#'
#' @importFrom ggplot2 ggplot aes_string geom_histogram ggtitle xlab theme_classic
#' @importFrom cowplot plot_grid
visHist <- function(cms_res, ncol = ifelse(length(colnames(cms_res)) > 1, 2, 1)){
  p <- do.call(plot_grid, c(lapply(colnames(cms_res), function(cms_name, ...){
    ggplot(as.data.frame(cms_res), aes_string(x=cms_name)) +
            geom_histogram(color="black", fill = col_hist[which(colnames(cms_res) %in% cms_name)]) +
            ggtitle(cms_name) + xlab(cms_name) + theme_classic()
  }), ncol = ncol))
  p
}


#' visOverview
#'
#' Plot an overview of cms, smoothened cms, group label and any colData variable in a reduced dimensional representation.
#'
#' @param cms_res Matrix or data.frame. Cms scores to plot should be in columns and cells in rows.
#' @param sce A \code{SingleCellExperiment} object with the combined data corresponding to 'cms_res'.
#' @param group Character. Name of group/batch variable. Needs to be one of \code{names(colData(sce))}.
#' @param dim_red Character. Name of embeddings to use as subspace for plotting. Default is "TSNE".
#' @param smooth Logical. Indicating if smoothened cms scores should be plotted.
#' @param log10_val Logical. Indicating if -log10(cms) should be plotted.
#' @param other_Var Character string. Name of other variables to be plotted asided.
#' Need correspond to one of \code{colData(sce)}.
#'
#' @details Plots reduced dimensions of cells colored by group variable and cms score.
#' If 'red_dim' is not defined in \code{reducedDimNames(sce)} a tsne is calculated using \code{runTSNE}.
#' Other color label as celltype label or smoothened cms scores can be plotted aside.
#' Embeddings from data integration methods (e.g. mnn.correct) can be used as long as they are specified in \code{reducedDimNames(sce)}.
#'
#' @family visualize cms functions
#' @seealso \code{\link{visCms}}, \code{\link{visGroup}}
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' load(system.file("extdata/cms_sim30.rda", package = "CellMixS"))
#' sce <- sim_30[[1]][, c(1:50, 500:550)]
#'
#' visOverview(cms_sim30, sce, "batch")
#'
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

    if(log10_val == TRUE){
      t_cms <- t_cms + guides(color=guide_legend(title="-log10(cms)"))
    }

    p <- plot_grid(t_group, t_cms)

    if(smooth == TRUE){
      df$cms_smooth <- cms_res[,"cms_smooth"]

      if(log10_val == TRUE){
        df$cms_smooth <- -log10(df$cms_smooth)
      }

      t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
        xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
        theme_void() +
        theme(aspect.ratio=1, panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color="grey", size=.3))

      t_smooth <- t + geom_point(size=1, alpha = 0.3, aes_string(color="cms_smooth")) +
        guides(color=guide_legend(override.aes=list(size=0.5))) +
        scale_color_viridis(option = "B") + ggtitle(paste0("smooth cms : ", group))

      if(log10_val == TRUE){
        t_smooth <- t_smooth + guides(color=guide_legend(title="-log10(cms_smooth)"))
      }

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
#' Plot cms scores in a reduced dimensional plot.
#'
#' @param cms_res Matrix or data.frame. Cms scores to plot should be in columns and cells in rows.
#' @param sce  A \code{SingleCellExperiment} object with the combined data corresponding to 'cms_res'.
#' @param cms_var Character Name of the cms scores to use (usually "cms_smoothened" or "cms").
#' @param dim_red Character. Name of embeddings to use as subspace for plotting. Default is "TSNE".
#' @param log10_val Logical. Indicating if -log10(cms) should be plotted.
#' @param ... Additional arguments to pass to \code{\link{ggplot}}.
#'
#' @details Plots a reduced dimension plot colored by cms scores.
#' The dimesion reduction embedding can be specified, but only tsne embeddings will automatically be computed using \code{runTSNE}.
#' Embeddings from data integration methods (e.g. mnn.correct) can be used as long as they are specified in \code{reducedDimNames(sce)}.
#'
#' @seealso \code{\link{visOverview}}, \code{\link{visGroup}}
#' @family visualize cms functions
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' load(system.file("extdata/cms_sim30.rda", package = "CellMixS"))
#' sce <- sim_30[[1]][, c(1:50, 500:550)]
#'
#' visCms(cms_res = cms_sim30, sce = sce)
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

  if(log10_val == TRUE){
    t_cms <- t_cms + guides(color=guide_legend(title=paste0("-log10(", cms_var, ")")))
  }

  t_cms
}


#' visGroup
#'
#' Plot group label in a reduced dimensional plot.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param group Character. Name of group/batch variable. Needs to be one of \code{names(colData(sce))}.
#' @param dim_red Character. Name of embeddings to use as subspace for plotting. Default is "TSNE".
#' @param ... Additional arguments to pass to \code{\link{ggplot}}.
#'
#' @details Plots a reduced dimension plot colored by group parameter.
#' The dimesion reduction embedding can be specified, but only tsne embeddings will automatically be computed by \code{runTSNE}.
#'  Embeddings from data integration methods (e.g. mnn.correct) can be used as long as they are specified in \code{reducedDimNames(sce)}.
#'
#' @seealso \code{\link{visOverview}}, \code{\link{visCms}}
#' @family visualize cms functions
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' sce <- sim_30[[1]][, c(1:50, 500:550)]
#'
#' visGroup(sce, "batch")
#'
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



