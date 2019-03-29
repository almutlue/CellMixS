# Functions to summarize cms

### summarize function

#' cmsSummary
#'
#' Summarizes cms scores, if specified based on groups.
#'
#' @param cms_res data frame of cms scores to be summarized.
#' @param sum_var variable to group summaries. Need to correspond to one of \code{colData(sce)}.
#' @param sce A \code{SingleCellExperiment} object corresponding to 'cms_res'.
#'
#' @details Summarises cms scores by mean to make different conditions/methods/.. comparable.
#' In default the mean of each cloumn of cms_res is returned.
#' Groups to summarize can be specified by sum_var and sce.
#'
#' @family cms functions
#' @seealso \code{\link{compareIntegration}}, \code{\link{compareCluster}}
#'
#' @return Data frame with mean cms scores (for each group).
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' load(system.file("extdata/cms_sim30.rda", package = "CellMixS"))
#' sce <- sim_30[[1]][, c(1:50,500:550)]
#' cmsSummary(cms_sim30)
#' cmsSummary(cms_sim30, sum_var = "batch", sce = sce)
#'
#' @importFrom dplyr as_tibble group_by_ summarize_all funs
#' @importFrom SingleCellExperiment colData
#' @importFrom magrittr %>%
cmsSummary <- function(cms_res, sum_var = NULL, sce = NULL){
  if(is.null(sum_var)){
    cms_summarized <- as.data.frame(t(colMeans(as.data.frame(cms_res))))
  }else{
    if(is.null(sce)){
      stop("Missing variable: Please provide a 'sce' object.")
    }
    if(!sum_var %in% names(colData(sce))){
      stop("Missing variable: Could not find 'sum_var'.
           Please specify one of names(colData(sce)).")
    }
    #to prevent type conversion if cms has 1 column only
    cms_res_sorted <- as.data.frame(cms_res[colnames(sce),])
    colnames(cms_res_sorted) <- colnames(cms_res)
    rownames(cms_res_sorted) <- colnames(sce)
    cms_merged <- cbind.data.frame(cms_res_sorted, colData(sce)[,sum_var])
    colnames(cms_merged) <- c(colnames(cms_res_sorted), sum_var)

    cms_summarized <- as_tibble(cms_merged) %>% group_by_(sum_var) %>% summarize_all(funs(mean))
    cms_summarized <- as.data.frame(cms_summarized)
  }
  cms_summarized
}


## Plot summarized cms function

#' compareIntegration
#'
#' Creates a summary plot of cms scores (for different integration methods).
#'
#' @param cms_res data frame, matrix or list of cms scores to be summarized.
#' Each column of the dataframe or each element of the list corresponds to a set of cms score that shall be summarized.
#' List elements need to have the same length and only one set of cms scores per list element.
#' @param scale Scale param for \code{geom_density_ridges} from \code{ggridges}.
#' @param violin A logical. If true violin plots are plotted, while the default (FALSE) will plot ridge plots.
#'
#' @details Plots summarized cms scores from an input list or dataframe.
#' This function is intended to visualize and compare different methods and views of the same dataset, not to compare different datasets.
#'
#' @seealso \code{\link{compareCluster}}, \code{ggridges}
#' @family visualize cms functions
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' load(system.file("extdata/cms_sim30.rda", package = "CellMixS"))
#'
#' sce <- sim_30[[1]][, c(1:50,500:550)]
#' cms_mnn <- cms(sce, k = 30, group = "batch", dim_red = "MNN")
#' cms_list <- list("raw"= cms_sim30[,"cms"], "mnn" = cms_mnn[,"cms"])
#'
#' compareIntegration(cms_list)
#'
#' @importFrom ggplot2 ggplot aes ylab xlab scale_color_manual theme_classic labs geom_violin
#' @importFrom ggridges geom_density_ridges
#' @importFrom tidyr gather_
compareIntegration <- function(cms_res, scale = 1, violin = FALSE){
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
  if(violin == TRUE){
    summarized_cms <- ggplot(average_long, aes_string(x="alignment", y="average_metric", fill="alignment")) +
      geom_violin()  +
      labs(title="Summarized metric", x="alignment", y = "cms") +
      scale_fill_manual(values = col_hist) + theme_classic()
  }else{
    summarized_cms <- ggplot(average_long, aes_string(y="alignment", x="average_metric", fill="alignment")) +
      geom_density_ridges(scale = scale)  +
      labs(title="Summarized metric",y="alignment", x = "cms") +
      scale_fill_manual(values= col_hist) + theme_classic()
  }
  summarized_cms
}


#' compareCluster
#'
#' Creates summary plots of cms scores for different groups/cluster.
#'
#' @param cms_res data frame of cms scores to be summarized. Cms scores should be in the first column.
#' Either a second column with corresponding levels of \code{cluster_var} or \code{sce} with levels in \code{colData(sce)} need to be specified.
#' @param cluster_var Character. Name of the factor level variable to summarize cms scores on.
#' @param cms_var Character. Name of the cms_res to use. Default is "cms".
#' @param sce \code{SingleCellexperiment} object. Should only be specified if \code{cluster_var} is not already provided in \code{cms_res}.
#' @param violin A logical. If true violin plots are plotted, while the default (FALSE) will plot ridge plots.
#'
#' @details Plots summarized cms scores for a specified list.
#' This function is intended to visualize and compare cms scores among clusters or other dataset variables.
#'
#' @seealso \code{\link{compareIntegration}}
#' @family visualize cms functions
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' load(system.file("extdata/cms_sim30.rda", package = "CellMixS"))
#'
#' sce <- sim_30[[1]][, c(1:50,500:550)]
#' colData(sce)$group <- sample(c("A", "B", "C"), ncol(sce), replace = TRUE)
#' compareCluster(cms_sim30, "group", cms_var = "cms", sce = sce, violin = TRUE)
#'
#' @importFrom ggplot2 ggplot aes ylab xlab scale_fill_manual theme_classic labs geom_violin
#' @importFrom tidyr gather_
#' @importFrom SingleCellExperiment colData
#' @importFrom ggridges geom_density_ridges
compareCluster <- function(cms_res, cluster_var, cms_var = "cms", sce = NULL, violin = FALSE){
  if(!is.null(sce)){
    cms_res_sorted <- as.data.frame(cms_res[colnames(sce),])
    colnames(cms_res_sorted) <- colnames(cms_res)
    rownames(cms_res_sorted) <- colnames(sce)
    cms_table <- data.frame(cms = cms_res_sorted[,cms_var], cluster = as.factor(colData(sce)[,cluster_var]))
    colnames(cms_table) <- c(cms_var, cluster_var)
  }else{
    cms_table <- cms_res[,c(cms_var, cluster_var)]
    colnames(cms_table) <- c(cms_var, cluster_var)
    cms_table[,cluster_var] <- as.factor(cms_table[,cluster_var])
  }

  #plot
  if(violin == TRUE){
    summarized_cms <- ggplot(cms_table, aes_string(x=cluster_var, y=cms_var, fill=cluster_var)) +
      geom_violin()  +
      labs(title="Summarized cms", x=cluster_var, y = cms_var) +
      scale_fill_manual(values = col_hist) + theme_classic()

  }else{
    summarized_cms <- ggplot(cms_table, aes_string(y=cluster_var, x=cms_var, fill=cluster_var)) +
      geom_density_ridges(scale = 1)  +
      labs(title="Summarized cms", y=cluster_var, x = cms_var) +
      scale_fill_manual(values = col_hist) + theme_classic()
  }
  summarized_cms
}
