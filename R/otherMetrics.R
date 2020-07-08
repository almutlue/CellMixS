#Other metrics


#' entropy
#'
#' @param sce \code{SingleCellExperiment} object, with the integrated data.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}.
#' @param k Numeric. Number of k-nearest neighbours (knn) to use.
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

    n_batch <- length(levels(colData(sce)[,group]))

    #----------------- calculate entropy  ----------------------------------#
    entropy <- apply(knn[[group]], 1, function(x){
        freq_batch = table(x)/length(x)
        freq_batch_positive = freq_batch[freq_batch > 0]
        shannon <- -sum(freq_batch_positive * log(freq_batch_positive))/
            log(n_batch)
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
#' @param k Numeric. Number of k-nearest neighbours (knn) to use.
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
#' They provide a distance-based neighbourhood weightening in their Lisi package.
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
    #fixed neighbourhood of knn.
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
    knn[["weights"]] <- do.call("cbind", lapply(batch_ids, function(batch){
        weights <- knn[["p"]][[batch]][,"weights"]
    })) %>% set_colnames(batch_ids)
    knn[["weights"]] <- knn[["weights"]]/rowSums(knn[["weights"]]) *
        length(batch_ids)

    #Calculate p^2
    knn[["p2"]] <- do.call("cbind", lapply(batch_ids, function(batch){
        p <- knn[["p"]][[batch]][,"p"]
    })) %>% set_colnames(batch_ids)
    if( weight ){
        knn[["p2"]] <-(knn[["p2"]] * knn[["weights"]])^2
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



#' mixMetric
#'
#' @param sce \code{SingleCellExperiment} object, with the integrated data.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}.
#' @param k Numeric. Number of k-nearest neighbours (knn) to use.
#' @param dim_red Character. Name of embeddings to use as subspace for distance
#' distributions. Default is "PCA".
#' @param assay_name Character. Name of the assay to use for PCA.
#' Only relevant if no existing 'dim_red' is provided.
#' @param n_dim Numeric. Number of dimensions to include to define the subspace.
#' @param k_pos Position of the cell, which rank to use for scoring,
#' defaults to 5.
#' @param res_name Character. Appendix of the result score's name
#' (e.g. method used to combine batches).
#'
#'@details The mixMetric function implements the mixingMetric function from
#'Seurat (See \code{\link[Seurat]{MixingMetric}}. It takes the median rank of
#'the '__k_pos__ neighbour from each batch as estimation for the data's entropy
#'according to the batch variable. The same result can be assesed using the
#' \code{\link[Seurat]{MixingMetric}} function and a seurat object from the
#' __Seurat__ package.
#'
#' @return A \code{SingleCellExperiment} with the mixing metric within colData.
#' @export
#'
#' @references
#' Stuart T Butler A Hoffman P Hafemeister C Papalexi E et. al. (2019)
#' Comprehensive Integration of Single-Cell Data.
#' Cell.
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:15, 400:420, 16:30)]
#'
#' sce <- mixMetric(sce, "batch", k = 20)
#' @importFrom SingleCellExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom magrittr %>% set_rownames set_colnames set_names
#' @importFrom stats median
#' @importFrom purrr %>% map
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData<-
#' @importFrom dplyr bind_rows
mixMetric <- function(sce, group, k = 300, dim_red = "PCA",
                      assay_name = "logcounts", n_dim = 10, k_pos = 5,
                      res_name = NULL){
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

    #----------------- calculate mixing metric -----------------------#
    #We implement seurat's mixing metric to prevent dependencies on seurat.
    #Get the position of the 5th neighbour within each batch
    #Add cell itself to knn list (to reproduce seurat's results -
    #reflects not how mm is actually described)
    knn[["group_mod"]] <- cbind(colData(sce)[,group], knn[[group]])
    knn[["pos"]] <- lapply(batch_ids, function(batch){
        idx <- apply(knn[["group_mod"]], 1, function(cell){
            id <- which(cell %in% batch)
        })
        if( !is.list(idx) )
            idx <- split(idx, col(idx))
            names(idx) <- cell_names
        if( length(idx) == 0 ){
            pos <- rep(k, length(cell_names)) %>% set_names(cell_names)
        }else{
            pos <- lapply(cell_names, function(cell){
                if( length(idx[[cell]]) < k_pos ){
                    pos <- k
                }else{
                    pos <-idx[[cell]][k_pos]
                }
            })
        }
    }) %>% bind_rows() %>% t() %>%  set_colnames(batch_ids)

    #Get median
    knn[["mm"]] <- lapply(cell_names, function(cell){
        cell_med <- median(knn[["pos"]][cell,])
    }) %>% bind_rows() %>% t() %>% set_rownames(cell_names) %>%
        set_colnames("median")


    #Add to colData of sce
    mm <- as.data.frame(knn[["mm"]]) %>% set_colnames("mm")
    if( !is.null(res_name) ){
        mm <- mm %>% set_colnames(res_name)
    }
    stopifnot(all(rownames(mm) == colnames(sce)))
    colData(sce) <- cbind(colData(sce), mm)
    sce
}


#' locStructure
#'
#' @param sce \code{SingleCellExperiment} object, with the integrated data.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}.
#' @param dim_combined Charactyer. Name of the reduced dimensional
#' representation of the integrated data.
#' Needs to be one of \code{reducedDimNames(sce))}.
#' @param k Numeric. Number of k-nearest neighbours (knn) to use.
#' @param dim_red Character. Name of embeddings to calculate neighbourhoods
#' before integration. Default is "PCA".
#' @param assay_name Character. Name of the assay to use for PCA of the original
#'  (not integrated) data. Should not refer to "corrected" counts.
#' @param n_dim Numeric. Number of dimensions to include for the original data.
#' @param n_combined Numeric. Number of dimensions to include for the
#' integrated data.
#' @param res_name Character. Appendix of the result score's name
#' (e.g. method used to combine batches).
#'
#'@details The locStructure function implements the localStructure function
#'from Seurat (See \code{\link[Seurat]{LocalStruct}}. For each group it
#'calculates the k nearest neighbour within PCA space before integration and
#'compares it to the knn within the reduced dimensional representation after
#'integration. The score represents the proportion of overlapping neighbours.
#'The \code{\link[Seurat]{LocalStruct}} function is based on the
#'\code{\link[Seurat]{RunPCA}} function, while here \code{\link[scater]{runPCA}}
#'is used. This can cause small deviance from the
#'\code{\link[Seurat]{LocalStruct}} function, but overall these functions are
#' equivalent.
#'
#' @return A \code{SingleCellExperiment} with the mixing metric within colData.
#' @export
#'
#' @references
#' Stuart T Butler A Hoffman P Hafemeister C Papalexi E et. al. (2019)
#' Comprehensive Integration of Single-Cell Data.
#' Cell.
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[["batch20"]][, c(1:50, 300:350)]
#'
#' sce <- locStructure(sce, "batch", "MNN", k = 20, assay_name = "counts")
#' @importFrom SingleCellExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom magrittr %>% set_rownames set_colnames set_names
#' @importFrom purrr %>% map
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData<-
#' @importFrom dplyr bind_rows
#' @importFrom scater runPCA
locStructure <- function(sce, group, dim_combined, k = 100, dim_red = "PCA",
                      assay_name = "logcounts", n_dim = 10, n_combined = 10,
                      res_name = NULL){
    #------------------Check input parameter ---------------------------------#
    if( !is(sce, "SingleCellExperiment") ){
        stop("Error: 'sce' must be a 'SingleCellExperiment' object.")
    }
    if( !group %in% names(colData(sce)) ){
        stop("Error: 'group' variable must be in 'colData(sce)'")
    }
    if( k > min(table(droplevels(colData(sce)[, group]))) ){
        stop("'k' exceeds number of cells/batch.")
    }
    if( !assay_name %in% names(assays(sce)) ){
        stop("Error: 'assay_name' not found.")
    }
    if( !dim_combined %in% reducedDimNames(sce) ){
        stop("Error: 'dim_combined' not found.")
    }
    if( ncol(reducedDims(sce)[[dim_combined]]) < n_combined ){
        stop("Error: 'n_combined' exceeds the dimensions of 'dim_combined'.")
    }
    if( !is(colData(sce)[,group], "factor") ){
        sce[[group]] <- as.factor(colData(sce)[, group])
    }
    colData(sce)[, group] <- droplevels(colData(sce)[, group])

    #----------------- determine subspace ----------------------------------#
    cell_names <- colnames(sce)
    names(cell_names) <- cell_names
    batch_ids <- levels(colData(sce)[,group])

    subspace_com <- reducedDims(sce)[[dim_combined]][, seq_len(n_combined)]


    #----------------- calculate local Structure -----------------------#
    #We implement seurat's mixing metric to prevent dependencies on seurat.
    #Overlap between knn_combined and knn within each group
    gr_overlap <- lapply(batch_ids, function(batch){
        cell_group <- colnames(sce[,colData(sce)[,group] == batch])
        #calculate PCA on these batch cells only
        sce_sub <- sce[,cell_group]
        sce_sub <- runPCA(sce_sub, ncomponents = n_dim, scale = TRUE,
                          exprs_values = assay_name, ntop = 2000)
        #filter subspace by group
        subspace_b <- reducedDims(sce_sub)[["PCA"]]
        subspace_com_b <- subspace_com[cell_group, ]
        #determine knn
        knn <- findKNN(subspace_b, k = k) %>%
            map(set_rownames, rownames(subspace_b))
        knn_combined <- findKNN(subspace_com_b, k = k) %>%
            map(set_rownames, rownames(subspace_com_b))
        #determine overlap
        n_overlap <- vapply(cell_group, function(cell){
            length(intersect(knn[["index"]][cell,],
                             knn_combined[["index"]][cell,]))/k
        }, 1)
    }) %>% unlist() %>% data.frame() %>% set_colnames("overlap")


    #Add to colData of sce
    gr_overlap <- as.data.frame(gr_overlap[cell_names,]) %>%
        set_colnames("overlap") %>% set_rownames(cell_names)
    if( !is.null(res_name) ){
        gr_overlap <- gr_overlap %>% set_colnames(res_name)
    }
    stopifnot(all(rownames(gr_overlap) == colnames(sce)))
    colData(sce) <- cbind(colData(sce), gr_overlap)
    sce
}


