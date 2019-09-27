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
#' @details The entropy function calculates the Shannon entropy of the group
#' variable within each cell's k-nearest neighbourhood.
#' For balanced batches a Shannon entropy close to 1 indicates high randomness
#' and mixing. For unbalanced batches entropy should be interpreted with
#' caution, but could work as a relative measure in a comparative setting.
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
    if( !is(sce, "SingleCellExperiment") ){
        stop("Error: 'sce' must be a 'SingleCellExperiment' object.")
    }
    if( !group %in% names(colData(sce)) ){
        stop("Error: 'group' variable must be in 'colData(sce)'")
    }
    if( !is(colData(sce)[,group], "factor") ){
        sce[[group]] <- as.factor(colData(sce)[, group])
    }
    colData(sce)[, group] <- droplevels(colData(sce)[, group])
    if( k >= ncol(sce) ){
        warning("'k' exceeds number of cells. Is set to max (all cells).")
        k <- ncol(sce) - 1
    }
    #----------------- determine knn matrix ----------------------------------#
    cell_names <- colnames(sce)
    names(cell_names) <- cell_names

    subspace <- .defineSubspace(sce, assay_name, dim_red, n_dim)
    #determine knn
    knn <- findKNN(subspace, k = k) %>% map(set_rownames, cell_names)

    #group assignment of knn cells for each cell
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

    stopifnot(all(rownames(entropy) == colnames(sce)))

    colData(sce) <- cbind(colData(sce), entropy)
    sce
}


#' isi
#'
#' @param sce \code{SingleCellExperiment} object, with the integrated data.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}.
#' @param k Numeric. Number of k-nearest neighbours (Knn) to use.
#' @param dim_red Character. Name of embeddings to use as subspace for distance
#' distributions. Default is "PCA".
#' @param assay_name Character. Name of the assay to use for PCA.
#' Only relevant if no existing 'dim_red' is provided.
#' @param n_dim Numeric. Number of dimensions to include to define the subspace.
#' @param weight Boolean. If TRUE, batch probabilities to calculate the isi
#'  score are weighted by the mean distance of their cells towards the cell
#'  of interest. Relevant for metrics: 'isi'.
#' @param res_name Character. Appendix of the result score's name
#' (e.g. method used to combine batches).
#'
#'@details The isi function calculates the inverse Simpson index of the group
#' variable within each cell's k-nearest neighbourhood.
#' The Simpson index describes the probability that two entities are taken at
#' random from the dataset and its inverse represent the effective number of
#' batches in a neighbourhood. The inverse Simpson index has been proposed as a
#' diversity score for batch mixing in single cell RNAseq by Korunsky et al.
#' They provide a distance-based neighborhood weightening in their Lisi package.
#' Here, we provide a simplified way of weightening probabilitities, if the 
#' \code{weight} argument is enabled.
#' @return A \code{SingleCellExperiment} with the entropy score within colData.
#' @export
#'
#' @references
#' Korsunsky I Fan J Slowikowski K Zhang F Wei K et. al. (2018).
#' Fast, sensitive, and accurate integration of single cell data with Harmony.
#' bioRxiv (preprint)
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:15, 400:420, 16:30)]
#'
#' sce <- isi(sce, "batch", k = 20)
#' @importFrom SingleCellExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom magrittr %>% set_rownames set_colnames set_names
#' @importFrom purrr %>% map
#' @importFrom BiocGenerics table
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData<-
#' @importFrom dplyr bind_rows
isi <- function(sce, group, k, dim_red = "PCA", assay_name = "logcounts",
                 n_dim = 10, weight = TRUE, res_name = NULL){
    #------------------Check input parameter ---------------------------------#
    if( !is(sce, "SingleCellExperiment") ){
        stop("Error: 'sce' must be a 'SingleCellExperiment' object.")
    }
    if( !group %in% names(colData(sce)) ){
        stop("Error: 'group' variable must be in 'colData(sce)'")
    }
    if( !is(colData(sce)[,group], "factor") ){
        sce[[group]] <- as.factor(colData(sce)[, group])
    }
    colData(sce)[, group] <- droplevels(colData(sce)[, group])
    if( k >= ncol(sce) ){
        warning("'k' exceeds number of cells. Is set to max (all cells).")
        k <- ncol(sce) - 1
    }
    #----------------- determine knn matrix ----------------------------------#
    cell_names <- colnames(sce)
    names(cell_names) <- cell_names
    batch_ids <- levels(colData(sce)[,group])

    subspace <- .defineSubspace(sce, assay_name, dim_red, n_dim)
    #determine knn
    knn <- findKNN(subspace, k = k) %>% map(set_rownames, cell_names)

    # group assignment of knn cells for each cell
    knn[[group]] <- matrix(
        colData(sce)[as.numeric(as.character(knn$index)), group],
        nrow = nrow(knn$index)) %>%
        set_rownames(cell_names)

    #----------------- calculate Inverse simpson index -----------------------#
    #We simplify the LISI score implemented in harmony by using a
    #fixed neighborhood of knn.
    knn[["p"]] <- lapply(batch_ids, function(batch){
        idx <- apply(knn[[group]], 1, function(cell){
            id <- which(cell %in% batch)
        })
        if( length(idx) == 0 ){
            p <- rep(0, length(cell_names)) %>% set_names(cell_names)
            weights <- rep(0, length(cell_names)) %>% set_names(cell_names)
            summary <- cbind(p, weights)
        }else{
            p <- lapply(cell_names, function(cell){
                if( length(idx[[cell]]) == 0 ){
                    summary <- c(0, 0)
                }else{
                    p_1 <- length(idx[[cell]])/k
                    p_cell <- p_1
                    weights <- 1/(knn[[2]][cell, idx[[cell]]] + 1)
                    summary <- c(p_cell, mean(weights))
                }
            }) %>% bind_cols() %>% t() %>% set_colnames(c("p", "weights"))
        }
    }) %>% set_names(batch_ids)

    #Get distance based weights
    knn[["weights"]] <- lapply(batch_ids, function(batch){
        weights <- knn[["p"]][[batch]][,"weights"]
    }) %>% set_names(batch_ids) %>% bind_rows
    knn[["weights"]] <- knn[["weights"]]/rowSums(knn[["weights"]]) *
        length(batch_ids)

    #Calculate p^2
    knn[["p2"]] <- lapply(batch_ids, function(batch){
        p <- knn[["p"]][[batch]][,"p"]
    }) %>% set_names(batch_ids) %>% bind_rows
    if( weight ){
        knn[["p2"]] <-(data.frame(knn[["p2"]]) * data.frame(knn[["weights"]]))^2
    }else{
        knn[["p2"]] <- knn[["p2"]]^2
    }
    rownames(knn[["p2"]]) <- cell_names

    #Inverse Simpson
    knn[["in_simpson"]] <- 1/rowSums(knn[["p2"]])

    #Add to colData of sce
    in_simpson <- as.data.frame(knn[["in_simpson"]]) %>% set_colnames("isi")
    if( !is.null(res_name) ){
        in_simpson <- in_simpson %>% set_colnames(res_name)
    }
    stopifnot(all(rownames(in_simpson) == colnames(sce)))
    colData(sce) <- cbind(colData(sce), in_simpson)
    sce
}
