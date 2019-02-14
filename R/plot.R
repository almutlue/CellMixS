#Plot functions


#' plot.hist
#'
#' Plot pvalue histograms of cms score distributions
#'
#' @param cms cms result matrix with one or two columns containing cms scores and/or scaled cms scores.
#' @param ... Additional arguments to pass to \code{\link{ggplot2}}.
#'
#' @details Plots cms score distribution similar to a pvalue histogram distribution.
#' Scaled cms should be close to a normal distribution centered around 0.5 and unscaled should be approx. flat, if no dataset-specific bias are expected.
#'
#' @family plot cms functions
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom ggplot2 ggplot aes_string geom_histogram ggtitle xlab theme_classic
#' @importFrom gridExtra grid.arrange
plot.hist <- function(cms, ...){
  p <- do.call(grid.arrange, c(lapply(colnames(cms), function(cms_name, ...){
    ggplot(as.data.frame(cms), aes_string(x=cms_name)) +
            geom_histogram(color="black", fill = col_hist[which(colnames(cms) %in% cms_name)]) +
            ggtitle(cms_name) + xlab(cms_name) + theme_classic()
  }), ncol = ifelse(length(colnames(cms)) > 1, 2, 1)))
  p
}


#' plot.overview
#'
#' Plot cms, scaled cms, group label and other Variable defined in colData in a reduced dimensional representation.
#'
#' @param cms_res cms result matrix with one or two columns containing cms scores and/or scaled cms scores.
#' Rownames need to correspond to colnames of sce.
#' @param sce Combined \code{\link{SingleCellExperiment}} object containing either precalculated reduced dimension embeddings or counts or logcounts of all groups merged.
#' If precalculated dimension reduction embeddings are used, they need to be within the reducedDimensions slot.
#' @param group Character string specifying the name of variable used to define groups (batches). Should be have a corresonding element in colData(sce).
#' @param dim_red Character defining dimension reduction embeddings to plot. Default is TSNE as from \code{\link{runTSNE}}.
#' If none is provided and no dimension reduction slot named "TSNE" is within \code{\link{reducedDimNames}} tsne embeddings will be calculated by \code{\link{runTSNE}}.
#' @param scaled A logical value indicating if scaled cms scores should be plotted.
#' If true cms_res needs to contain 2 colums with cms and cms_scaled results. Default is set by the length of colnames(cms_res)
#' @param log10_val A logical indicating if -log10(cms) should be plotted to visualize differences of small values
#' @param other_Var Character string defining other variables to be plotted asided (e.g. some specified celltype label).
#' Need to be specified in colData(sce).
#' @param ... Additional arguments to pass to \code{\link{runTSNE}} (CHECK!!!)
#'
#' @details Plots reduced dimensions (tsne as default) of cells using the group variable and the cms score as colors.
#' Other color label as celltype label or scaled cms scores can be plotted aside. Generates tsne embeddings, if none have been specified.
#' Embeddings from data integration methods (e.g. mnn.correct) can also been used as long as they are specified in \code{\link{reducedDimNames}}.
#'
#' @family plot cms functions
#' @seealso \code{\link{plot.cms}}, \code{\link{plot.group}}
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
plot.overview <- function(cms_res, sce, group, dim_red = "TSNE", scaled = ifelse(length(colnames(cms_res)) > 1, TRUE, FALSE), log10_val = FALSE, other_Var = NULL, ...){
  cell_names <- colnames(sce)
  #Compare order sce and cms_restheme_classic
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

    if(scaled == TRUE){
      df$cms_scaled <- cms_res[,"cms_scaled"]

      t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
        xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
        theme_void() +
        theme(aspect.ratio=1, panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color="grey", size=.3))

      t_scaled <- t + geom_point(size=1, alpha = 0.3, aes_string(color="cms_scaled")) +
        guides(color=guide_legend(override.aes=list(size=0.5))) +
        scale_color_viridis(option = "B") + ggtitle(paste0("scaled cms : ", group))
      p <- plot_grid(t_group, t_cms, t_scaled)

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

      if(scaled == TRUE){
      p <- plot_grid(t_group, t_cms, t_scaled, t_other)
      }else{
        p <- plot_grid(t_group, t_cms, t_other)
      }

    }

    #plot
    return(p)

}


### Single plots


#' plot.cms
#'
#' Plot cms scores in a reduced dimensional plot like tsne.
#'
#' @param cms_res cms result matrix with one or two columns containing cms scores and/or scaled cms scores.
#' Rownames need to correspond to colnames of sce. Colnames are allowed to be set different, but need to be specified in cms_var.
#' @param sce Combined \code{\link{SingleCellExperiment}} object containing either precalculated reduced dimension embeddings or counts or logcounts of all groups merged.
#' If precalculated dimension reduction embeddings are used, they need to be within the reducedDimensions slot.
#' @param cms_var character string specifying the cms scores to use (usually "cms_scaled" or "cms").
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
#' @seealso \code{\link{plot.overview}}, \code{\link{plot.group}}
#' @family plot cms functions
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
plot.cms <- function(cms_res, sce, cms_var = "cms", dim_red = "TSNE", log10_val = FALSE, ...){
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


#' plot.group
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
#' @seealso \code{\link{plot.overview}}, \code{\link{plot.cms}}
#' @family plot cms functions
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
plot.group <- function(sce, group, dim_red = "TSNE", ...){
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

## Plot summariezed cms function

#' compare.integration
#'
#' Creates a summary boxplot of cms scores (for different integration methods).
#'
#' @param cms_res data frame or list of cms scores to be summarized. Each column of the dataframe or each element of the list should correspond to a set of cms score that shall be summarized.
#' List elements need to have the same length and only one set of cms scores per list element is allowed.
#'
#' @details Plots summarized cms scores from an input list or dataframe. This function is intended to visualize and compare different methods and views of the same dataset, not to compare different datasets.
#'
#' @seealso \code{\link{compare.cluster}}
#' @family plot cms functions
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom ggplot2 ggplot aes ylab xlab scale_color_manual theme_classic geom_boxplot labs
#' @importFrom tidyr gather_
compare.integration <- function(cms_res){
    if(is.list(cms_res)){
      average_table <- do.call(cbind.data.frame, cms_res)
      colnames(average_table) <- names(cms_res)
    }else{
      average_table <- as.data.frame(cms_res)
    }


    #change to long format
    #long format
    keycol <- "alignment"
    valuecol <- "average_metric"
    gathercols <- colnames(average_table)
    average_long <- gather_(average_table, keycol, valuecol, gathercols, factor_key=TRUE)


    #plot
    summarized_cms <- ggplot(average_long, aes(x=alignment, y=average_metric, color=alignment)) +
      geom_boxplot()  +
      labs(title="Summarized cms",x="alignment", y = "cms") + scale_color_manual(values= col_hist) + theme_classic()

    summarized_cms
  }


#' compare.cluster
#'
#' Creates summary boxplots of cms scores for different groups/cluster.
#'
#' @param cms_res data frame of cms scores to be summarized. The cms score should be in the first coulmn. It should either contain a second column with corresponding levels for cluster_var or  sce need to be specified.
#' @param cluster_var character string specifying the name of the factor level variable to summarize cms scores on.
#' @param sce SingleCellexperiment object should only be specified if cluster-var is not already provided in cms_res. Should contain 'cluster_var' within colData.
#' If sce is specified, rownames of cms_res need to correspond to colnames of sce (should contain the same cells).
#'
#' @details Plots summarized cms scores for a specified list. This function is intended to visualize and compare cms scores among clusters or other dataset variables.
#'
#' @seealso \code{\link{compare.integration}}
#' @family plot cms functions
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom ggplot2 ggplot aes ylab xlab scale_fill_manual theme_classic geom_boxplot labs
#' @importFrom tidyr gather_
#' @importFrom SingleCellExperiment colData
compare.cluster <- function(cms_res, cluster_var, cms_var = "cms", sce = NULL){
  if(!is.null(sce)){
    cms_res <- cms_res[colnames(sce),]
    cms_table <- data.frame(cms = cms_res[,cms_var], cluster = as.factor(colData(sce)[,cluster_var]))
  }else{
    cms_table <- cms_res
    if(ncol(cms_table)!= 2){
     stop("Input error: Please provide a data frame with two columns one with the cms score and one with the coresponding levels of cluster_var")
    }
    colnames(cms_table) <- c("cms", "cluster")
    cms_table$cluster <- as.factor(cms_table$cluster)
  }

 #plot
  summarized_cms <- ggplot(cms_table, aes(x=cluster, y=cms, fill=cluster)) +
    geom_boxplot()  +
    labs(title="Summarized cms", x=cluster_var, y = "cms") + scale_fill_manual(values = col_hist) + theme_classic()

  summarized_cms
}

