# Functions to plot summarized metrics

### summarize function

#' visIntegration
#'
#' Creates a summary plot of metric scores (for different integration methods).
#'
#' @param res_object \code{SingleCellExperiment} object, list, matrix or
#' data.frame. The SingleCellExperiment object should contain the result scores
#' (cms) to compare within \code{colData(res_object)}. List, matrix or data
#' frame should have result scores in list elements resp. columns.
#' @param metric_prefix Character. Prefix to specify names of
#' \code{colData(sce)} to be compared. Applys only if `res_object` is a
#' \code{SingleCellExperiment} object. Default is 'cms'.
#' @param violin A logical. If true violin plots are plotted,
#' while the default (FALSE) will plot ridge plots.
#'
#' @details Plots summarized cms scores from an \code{SingleCellExperiment}
#' object, list or dataframe. This function is intended to visualize and
#' compare different methods and views of the same dataset, not to compare
#' different datasets.
#'
#' @seealso \code{\link{visCluster}}, \code{ggridges}
#' @family visualize  functions
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#'
#' sce <- sim_list[["batch20"]][, c(1:30,300:320)]
#' sce_mnn <- cms(sce,"batch", k = 20, dim_red = "MNN", res_name = "MNN",
#' n_dim = 2)
#'
#' visIntegration(sce_mnn, metric_prefix = "cms.", violin = TRUE)
#'
#' @importFrom ggplot2 ggplot aes ylab xlab scale_color_manual theme_classic
#' labs geom_violin
#' @importFrom ggridges geom_density_ridges
#' @importFrom tidyr gather
#' @importFrom dplyr as_tibble select starts_with as_data_frame
#' @importFrom magrittr %>%
#' @importFrom methods is
visIntegration <- function(res_object, metric_prefix = "cms", violin = FALSE){
    if( is.list(res_object) ){
        average_table <- res_object %>% cbind.data.frame()
    }else if( is(res_object, "SingleCellExperiment") ){
        average_table <- as_tibble(colData(res_object)) %>%
            select(starts_with(metric_prefix))
    }else{
        average_table <- as_data_frame(res_object)
    }

    #change to long format
    #long format
    gathercols <- colnames(average_table)
    average_long <- gather(average_table, "keycol", "valuecol", gathercols,
                           factor_key=TRUE)


    #plot
    if( isTRUE(violin) ){
        summarized_metric <- ggplot(average_long, aes_string(x="keycol",
                                                             y="valuecol",
                                                             fill="keycol")) +
            geom_violin()  +
            labs(title="Summarized metric", x="integration",
                 y = paste0("average ", metric_prefix)) +
            scale_fill_manual(values = col_hist) + theme_classic()
    }else{
        summarized_metric <- ggplot(average_long, aes_string(y="keycol",
                                                             x="valuecol",
                                                             fill="keycol")) +
            geom_density_ridges(scale = 1)  +
            labs(title="Summarized metric",y="integration",
                 x = paste0("average ", metric_prefix)) +
            scale_fill_manual(values= col_hist) + theme_classic()
    }
    summarized_metric
}


#' visCluster
#'
#' Creates summary plots of metric scores for different groups/cluster.
#'
#' @param sce_cms A \code{SingleCellExperiment} object with the result scores
#' (e.g. cms) to plot within \code{colData(res_object)}.
#' @param cluster_var Character. Name of the factor level variable to summarize
#' metric scores on.
#' @param metric_var Character Name of the metric scores to use.
#' Default is "cms".
#' @param violin A logical. If true violin plots are plotted, while the default
#' (FALSE) will plot ridge plots.
#'
#' @details Plots summarized metric scores.
#' This function is intended to visualize and compare metric scores among
#' clusters or other dataset variables spcified in `cluster_var`.
#'
#' @seealso \code{\link{visIntegration}}
#' @family visualize functions
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:30,300:320)]
#' sce_cms <- cms(sce, "batch", k = 20, n_dim = 2)
#'
#' visCluster(sce_cms, "batch")
#'
#' @importFrom ggplot2 ggplot aes ylab xlab scale_fill_manual theme_classic labs
#'  geom_violin
#' @importFrom SingleCellExperiment colData
#' @importFrom ggridges geom_density_ridges
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate_at as_tibble
visCluster <- function(sce_cms, cluster_var, metric_var = "cms",
                       violin = FALSE){

    metric_table <- as_tibble(colData(sce_cms)) %>%
        select(metric_var, cluster_var) %>%
        mutate_at(cluster_var, as.factor)

    #plot
    if(violin == TRUE){
        summarized_metric <- ggplot(metric_table,
                                    aes_string(x=cluster_var, y=metric_var,
                                               fill=cluster_var)) +
            geom_violin() +
            labs(title=paste0("Summarized ", metric_var), x=cluster_var,
                 y = metric_var) +
            scale_fill_manual(values = col_hist) + theme_classic()

    }else{
        summarized_metric <- ggplot(metric_table,
                                    aes_string(y=cluster_var, x=metric_var,
                                               fill=cluster_var)) +
            geom_density_ridges(scale = 1)  +
            labs(title=paste0("Summarized ", metric_var), y=cluster_var,
                 x = metric_var) +
            scale_fill_manual(values = col_hist) + theme_classic()
    }
    summarized_metric
}
