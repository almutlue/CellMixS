#Other metrics


#' entropy
#'
#' @param sce \code{SingleCellExperiment} object, with the integrated data.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}.
#' @param k Numeric. Number of k-nearest neighbours (Knn) to use.
#' @param dim_red Character. Name of embeddings to use as subspace for distance
#' distributions. Default is "PCA".
#' @param assay_name Character. Name of the assay to use for PCA.
#' Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce))}. Default is "logcounts".
#' @param n_dim Numeric. Number of dimensions to include to define the subspace.
#' @param res_name Character. Appendix of the result score's name
#' (e.g. method used to combine batches).
#'
#' @details The entropy function calculates the shannon entropy of the group
#' variable within each cell's k-nearest neighbourhood.
#' For balanced batches the shannon entropy close to 1 indicates high randomness
#' and mixing. For unbalanced batches entropy should be interpreted with
#' caution, but could work as a relative meassure in a comparative setting.
#'
#'
#' @return A \code{SingleCellExperiment} with the entropy score within colData.
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:15, 400:420, 16:30)]
#'
#' sce <- entropy(sce, "batch", k = 20)
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom magrittr %>% set_rownames set_colnames
#' @importFrom purrr %>% map
#' @importFrom BiocGenerics table
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData<-
entropy <- function(sce, group, k, dim_red = "PCA", assay_name = "logcounts",
                     n_dim = 10, res_name = NULL){
    #------------------Check input parameter ---------------------------------#
    if(!is(sce, "SingleCellExperiment")){
        stop("Error: 'sce' must be a 'SingleCellExperiment' object.")
    }
    if(!group %in% names(colData(sce))){
        stop("Error: 'group' variable must be in 'colData(sce)'")
    }
    if(!is(colData(sce)[,group], "factor")){
        sce[[group]] <- as.factor(colData(sce)[, group])
    }
    #----------------- determine knn matrix ----------------------------------#
    cell_names <- colnames(sce)
    names(cell_names) <- cell_names

    subspace <- .defineSubspace(sce, assay_name, dim_red, n_dim)
    #determine knn
    knn <- findKNN(subspace, k = k) %>% map(set_rownames, cell_names)

    # group assignment of knn cells for each cell
    knn[[group]] <- matrix(
        colData(sce)[as.numeric(as.character(knn$index)), group],
        nrow = nrow(knn$index)) %>%
        set_rownames(cell_names)

    #----------------- calculate entropy  ----------------------------------#
    entropy <- apply(knn[[group]], 1, function(x){
        freq_batch = table(x)/length(x)
        freq_batch_positive = freq_batch[freq_batch > 0]
        shannon <- -sum(freq_batch_positive * log(freq_batch_positive))/
            log(length(levels(as.factor(x))))
    })

    #Add to colData of sce
    entropy <- as.data.frame(entropy) %>% set_colnames("entropy")
    if( !is.null(res_name) ){
        entropy <- as.data.frame(entropy) %>% set_colnames(res_name)
    }

    colData(sce) <- cbind(colData(sce), entropy)
    sce
}
