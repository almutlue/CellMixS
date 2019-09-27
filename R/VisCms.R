# Visualize metrics and groups


#' visHist
#'
#' Plot pvalue histograms of metric score distributions
#'
#' @param res_object \code{SingleCellExperiment} object, matrix or data.frame.
#' The SingleCellExperiment object should contain the result scores (e.g. cms)
#' to plot in \code{colData(res_object)}.
#' Matrix or data frame should have result scores in columns and cells in rows.
#' @param metric Character vector. Specify names of \code{colData(sce)} to be
#' plotted. Applys only if `res_object` is a \code{SingleCellExperiment} object.
#' Default is 'cms'. If prefix is TRUE all columns starting with `metric` will
#' be plotted.
#' @param prefix Boolean. Is `metric` used to specify column's prefix(true) or
#' complete column names (False).
#' @param n_col Numeric. Number of columns of the pval histogram.
#' @param metric_prefix Former parameter to define prefix of the metric to be
#' plotted. Will stop and ask for the new syntax.
#'
#'
#' @details Plots metric score distribution similar to a pvalue histogram
#' distribution. Without dataset-specific bias, cms scores should be approx.
#' flat distributed. If `res_object` is a matrix or data.frame,
#' it will create a histogram for each column. If `res_object` is a
#' \code{SingleCellExperiment} object, it will create a histogram of all
#' \code{colData(res_object)} that start with or are specified in `metric`.
#'
#' @family visualize metric functions
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:50)]
#' sce_cms <- cms(sce, "batch", k = 20, n_dim = 2)
#' visHist(sce_cms)
#'
#'
#' @importFrom ggplot2 ggplot aes_string geom_histogram ggtitle xlab
#' theme_classic
#' @importFrom cowplot plot_grid
#' @importFrom dplyr as_tibble select starts_with
#' @importFrom SingleCellExperiment colData
#' @importFrom methods is
visHist <- function(res_object, metric = "cms", prefix = TRUE, n_col = 1,
                    metric_prefix = NULL){
    ## Check input structure and select data to plot
    if( !is.null(metric_prefix) ){
        stop("'metric_prefix' has been replaced by the parameter 'metric'.
             Please change it's name and check the man page.")
    }
    if( is(res_object, "SingleCellExperiment") ){
        #select columns to plot
        if( prefix ){
            cms_res <- as_tibble(colData(res_object)) %>%
                select(starts_with(metric))
        }else{
            cms_res <- as_tibble(colData(res_object)) %>%
                select(metric)
            }
    }else{
        cms_res <- res_object
    }

    #Check input/presence of score
    if( ncol(cms_res) == 0 ){
        stop("Error: 'res_object' does not contain any metric results.
             Please continue by one of:
             * Run `cms` on your SingleCellExperiment object before plotting.
             * Specify colData(res_object) column to plot by `metric`.
             * Specify a matrix with results to plot as `res_object`.")
    }

    #plot function
    p <- do.call(plot_grid, c(lapply(colnames(cms_res), function(cms_name){
        ggplot(as.data.frame(cms_res), aes_string(x=cms_name)) +
            geom_histogram(color="black",
                           fill=col_hist[which(colnames(cms_res) %in% cms_name)],
                           breaks=seq(0, 1, by=0.05)) + xlab(cms_name) +
            theme_classic()
    }), ncol = n_col))
    p
    }


#' visOverview
#'
#' Plot an overview of metric results, group label and any colData variable in
#' a reduced dimensional representation.
#'
#' @param sce_cms A \code{SingleCellExperiment} object with the result scores
#' (e.g. cms) to plot in \code{colData(sce_cms)}.
#' @param group Character. Name of group/batch variable. Needs to be one of
#' \code{names(colData(sce))}.
#' @param metric Character vector. Specify names of \code{colData(sce)} to be
#' plotted. Applys only if `res_object` is a \code{SingleCellExperiment} object.
#' Default is 'cms'. If prefix is TRUE all columns starting with `metric` will
#' be plotted.
#' @param prefix Boolean. Is `metric` used to specify column's prefix(true) or
#' complete column names (False).
#' @param dim_red Character. Name of embeddings to use as subspace for plotting.
#'  Default is "TSNE".
#' @param log10_val Logical. Indicating if -log10(metric) should be plotted.
#' @param other_var Character string. Name(s) of other variables to be plotted
#' asided. Need correspond to one of \code{colData(sce)}.
#' @param metric_prefix Former parameter to define prefix of the metric to be
#' plotted. Will stop and ask for the new syntax.
#'
#' @details Plots reduced dimensions of cells colored by group variable and
#' metric score. If 'red_dim' is not defined in \code{reducedDimNames(sce)} a
#' tsne is calculated using \code{runTSNE}. Other color label as celltype label
#' or smoothened scores can be plotted aside. Embeddings from data integration
#' methods (e.g. mnn.correct) can be used if they are specified in
#' \code{reducedDimNames(sce)}.
#'
#' @family visualize metric functions
#' @seealso \code{\link{visMetric}}, \code{\link{visGroup}}
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:30, 300:330)]
#' sce_cms <- cms(sce, "batch", k = 20, n_dim = 2)
#'
#' visOverview(sce_cms, "batch", other_var = "batch")
#'
#'
#' @importFrom ggplot2 ggplot aes_string ylab xlab theme_void theme guide_legend
#' guides element_blank element_line geom_point scale_color_manual
#' @importFrom cowplot plot_grid
#' @importFrom scater runTSNE
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDimNames reducedDim colData
#' @importFrom viridis scale_color_viridis
#' @importFrom purrr map negate
#' @importFrom dplyr bind_rows as_tibble select starts_with select_if bind_cols
#' @importFrom magrittr %>% set_colnames
#' @importFrom methods is
visOverview <- function(sce_cms, group, metric = "cms", prefix = TRUE,
                        dim_red = "TSNE", log10_val = FALSE, other_var = NULL,
                        metric_prefix = NULL){
    ## Check input structure and select data to plot
    if( !is.null(metric_prefix) ){
        stop("'metric_prefix' has been replaced by the parameter 'metric'.
             Please change it's name and check the man page.")
    }
    if( !is(sce_cms, "SingleCellExperiment") ){
        stop("Error:'sce_cms' must be a 'SingleCellExperiment' object.")
    }
    if( !group %in% names(colData(sce_cms)) ){
        stop("Error: 'group' variable must be in 'colData(sce_cms)'")
    }
    #select columns to plot
    if( prefix ){
        cms_res <- as_tibble(colData(sce_cms)) %>%
            select(starts_with(metric))
    }else{
        cms_res <- as_tibble(colData(sce_cms)) %>%
            select(metric)
    }

    #Check the presence of results
    if( ncol(cms_res) == 0 ){
        stop("Error: 'sce_cms' does not contain any metric results.
             Please continue by one of:
             * Run `cms` on your SingleCellExperiment object before plotting.
             * Specify a colData(res_object) column in `metric`.")
    }

    cell_names <- colnames(sce_cms)

    ### ------ Start Check and prepare dim reduction slot --------####

    if( !dim_red %in% "TSNE" ){
        if( !dim_red %in% reducedDimNames(sce_cms) ){
            stop("Ambigous parameter 'dim_red', provide one of:
                 * A dim_red method that is listed in reducedDimNames(sce_cms).
                 * Default('TSNE') will call runTSNE to calculate a subspace.")
        }
        red_dim <- as.data.frame(reducedDim(sce_cms, dim_red))
        }else{
            #used tsne from scater package (check for availability first)
            if( !"TSNE" %in% reducedDimNames(sce_cms) ){
                #use "logcounts" if availabe otherwise "counts"
                if( names(assays(sce_cms)) %in% "logcounts" ){
                    sce_cms <- runTSNE(sce_cms)
                }else{
                    sce_cms <- runTSNE(sce_cms, exprs_values = "counts")
                }
            }
            red_dim <- as.data.frame(reducedDim(sce_cms, "TSNE"))
        }
    colnames(red_dim) <- c("red_dim1", "red_dim2")

    ### ------ Finish Check and prepare dim reduction slot --------####
    other_var_tib <- as_tibble(colData(sce_cms)[,other_var]) %>%
        set_colnames(other_var)

    #data frame to plot
    df <- data.frame(sample_id = cell_names,
                     group_var = as.factor(colData(sce_cms)[, group]),
                     red_Dim1=red_dim$red_dim1,
                     red_Dim2=red_dim$red_dim2) %>%
        bind_cols(cms_res, other_var_tib)


    # use of -log10 values (for very low cms)
    if( isTRUE(log10_val )){
        mi_log10 <- function(x)(-log10(x))
        df <- df %>% mutate_at(colnames(cms_res), mi_log10)
    }

    t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
        xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
        theme_void() +
        theme(aspect.ratio=1, panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color="grey", size=.3))

    #plot function for discrete input
    t_discr <- function(group_var){
        t_group <- t +
            geom_point(size=1, alpha = 0.3, aes_string(color=group_var)) +
            guides(color=guide_legend(override.aes=list(size=1))) +
            scale_color_manual(values = col_group) +
            ggtitle(group_var)
    }

    #plot function for continous input
    t_cont <- function(cms_res_var){
        t_cms <- t +
            geom_point(size=1, alpha = 0.5, aes_string(color=cms_res_var)) +
            guides(color=guide_legend(override.aes=list(size=0.5))) +
            scale_color_viridis(option = "B") +
            ggtitle(paste0(cms_res_var, " : ", group))
        if( isTRUE(log10_val) & cms_res_var %in% colnames(cms_res) ){
            t_cms <- t_cms +
                guides(color=guide_legend(
                    title=paste0("-log10(",cms_res_var, ")")))
        }
        t_cms
    }

    t_group   <- "group_var" %>% map(t_discr)

    t_cms_res <- colnames(cms_res) %>% map(t_cont)

    if( !is.null(other_var) ){
        t_other_var_cont <- as_tibble(colData(sce_cms)[,other_var]) %>%
            set_colnames(other_var) %>%
            select_if(is.numeric) %>%
            colnames() %>%
            map(t_cont)

        t_other_var_dis <- as_tibble(colData(sce_cms)[,other_var]) %>%
            set_colnames(other_var) %>%
            select_if(negate(is.numeric)) %>%
            colnames() %>%
            map(t_discr)

        t_all <- c(t_group, t_cms_res, t_other_var_cont, t_other_var_dis)
    }else{
        t_all <- c(t_group, t_cms_res)
    }


    p <- plot_grid(plotlist =  t_all)

    return(p)

    }


### Single plots


#' visMetric
#'
#' Plot metric scores in a reduced dimensional plot.
#'
#' @param sce_cms A \code{SingleCellExperiment} object with the result scores
#' (e.g. cms) to plot within \code{colData(res_object)}.
#' @param metric_var Character Name of the metric scores to use.
#' Default is "cms".
#' @param dim_red Character. Name of embeddings to use as subspace for plotting.
#'  Default is "TSNE".
#' @param log10_val Logical. Indicating if -log10(metric) should be plotted.
#'
#' @details Plots a reduced dimension plot colored by metric scores.
#' The dimension reduction embedding can be specified, but only tsne embeddings
#' will automatically be computed using \code{runTSNE}. Embeddings from data
#' integration methods (e.g. mnn.correct) can be used as long as they are
#' present in \code{reducedDimNames(sce)}.
#'
#' @seealso \code{\link{visOverview}}, \code{\link{visGroup}}
#' @family visualize metric functions
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:30, 300:320)]
#' sce_cms <- cms(sce, "batch", k = 20, n_dim = 2)
#'
#' visMetric(sce_cms)
#'
#' @importFrom ggplot2 ggplot aes_string ylab xlab theme_void theme guide_legend
#' guides element_blank element_line geom_point
#' @importFrom scater runTSNE
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDimNames reducedDim colData
#' @importFrom viridis scale_color_viridis
#' @importFrom dplyr mutate_at
#' @importFrom magrittr %>%
#' @importFrom methods is
visMetric<- function(sce_cms, metric_var = "cms", dim_red = "TSNE",
                     log10_val = FALSE){

    ## Check input structure
    if( !is(sce_cms, "SingleCellExperiment") ){
        stop("Error:'sce_cms' must be a 'SingleCellExperiment' object.")
    }
    if( !metric_var %in% names(colData(sce_cms)) ){
        stop("Error: 'metric_var' variable must be in 'colData(sce_cms)'")
    }

    ## Check reduced dimensions
    cell_names <- colnames(sce_cms)
    if(!dim_red %in% "TSNE"){
        if(!dim_red %in% reducedDimNames(sce_cms)){
            stop("Ambigous parameter 'dim_red', provide one of:
                 * A dim_red method that is listed in reducedDimNames(sce_cms).
                 * Default('TSNE') will call runTSNE to calculate a subspace.")
        }
        red_dim <- as.data.frame(reducedDim(sce_cms, dim_red))
        }else{
            if(is.null(reducedDim(sce_cms, "TSNE"))){
                if(names(assays(sce_cms)) %in% "logcounts"){
                    sce_cms <- runTSNE(sce_cms)
                }else{
                    sce_cms <- runTSNE(sce_cms, exprs_values = "counts")
                }
            }
            red_dim <- as.data.frame(reducedDim(sce_cms, "TSNE"))
        }
    colnames(red_dim) <- c("red_dim1", "red_dim2")

    #data frame to plot
    df <- data.frame(sample_id = cell_names,
                     metric = colData(sce_cms)[,metric_var],
                     red_Dim1=red_dim$red_dim1, red_Dim2=red_dim$red_dim2)

    # use of -log10 values (for very low metric scores)
    if( isTRUE(log10_val )){
        mi_log10 <- function(x)(-log10(x))
        df <- df %>% mutate_at("metric", mi_log10)
    }
    t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
        xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
        theme_void() +
        theme(aspect.ratio=1, panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color="grey", size=.3))

    t_metric <- t +
        geom_point(size=1, alpha = 0.5, aes_string(color="metric")) +
        guides(color= guide_legend(override.aes=list(size=2),
                                title = metric_var)) +
        scale_color_viridis(option = "B") + ggtitle(metric_var)

    if(isTRUE(log10_val)){
        t_metric <- t_metric +
            guides(color=guide_legend(title=paste0("-log10(", metric_var, ")")))
    }
    t_metric
    }


#' visGroup
#'
#' Plot group label in a reduced dimensional plot.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}.
#' @param dim_red Character. Name of embeddings to use as subspace for plotting.
#' Default is "TSNE".
#'
#' @details Plots a reduced dimension plot colored by group parameter.
#' The dimesion reduction embedding can be specified, but only tsne embeddings
#' will automatically be computed by \code{runTSNE}. Embeddings from data
#' integration methods (e.g. mnn.correct) can be used as long as they are
#' specified in \code{reducedDimNames(sce)}.
#'
#' @seealso \code{\link{visOverview}}, \code{\link{visMetric}}
#' @family visualize functions
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:50, 300:350)]
#'
#' visGroup(sce, "batch")
#'
#'
#' @importFrom ggplot2 ggplot aes_string ylab xlab theme_void theme guide_legend
#'  guides element_blank element_line geom_point scale_color_manual
#' @importFrom scater runTSNE
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDimNames reducedDim colData
#' @importFrom viridis scale_color_viridis
#' @importFrom methods is
visGroup <- function(sce, group, dim_red = "TSNE"){

    ## Check input structure
    if( !is(sce, "SingleCellExperiment") ){
        stop("Error:'sce' must be a 'SingleCellExperiment' object.")
    }
    if( !group %in% names(colData(sce)) ){
        stop("Error: 'group' variable must be in 'colData(sce)'")
    }

    #Generate or specify dim reduction
    cell_names <- colnames(sce)

    ### ------ Start Check and prepare dim reduction slot --------####
    if(!dim_red %in% "TSNE"){
        if(!dim_red %in% reducedDimNames(sce)){
            stop("Please provide a dim_red method listed in reducedDims of sce")
        }
        red_dim <- as.data.frame(reducedDim(sce, dim_red))
    }else{
        #use tsne from scater package (check for availability first)
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
    ### ------ Finish Check and prepare dim reduction slot --------####

    #data frame to plot
    df <- data.frame(sample_id = cell_names,
                     group_var = colData(sce)[, group],
                     red_Dim1=red_dim$red_dim1,
                     red_Dim2=red_dim$red_dim2)

    t <- ggplot(df, aes_string(x="red_Dim1", y="red_Dim2")) +
        xlab(paste0(dim_red,"_1")) + ylab(paste0(dim_red,"_2")) +
        theme_void() +
        theme(aspect.ratio=1, panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color="grey", size=.3))

    t_group <- t +
        geom_point(size=1, alpha = 0.5, aes_string(color="group_var")) +
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



