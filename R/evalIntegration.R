# Wrapper to run different evaluation metrics

#' evalIntegration
#'
#' Function to evaluate sc data integration providing a framework for different
#' metrics. Metrics to evaluate mixing and preservance of the local/individual
#' structure are provided.
#'
#' @param metrics Character vector. Name of the metrics to apply. Must be one to
#' all of 'cms', 'ldfDiff', 'isi', 'mixing_metric', 'local_structure',
#' 'entropy'.
#' @param sce \code{SingleCellExperiment} object, with the integrated data.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}.
#' @param dim_red Character. Name of embedding to use as subspace for distance
#' distributions. Default is "PCA".
#' @param assay_name Character. Name of the assay to use for PCA.
#' Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce))}. Default is "logcounts".
#' @param n_dim Numeric. Number of dimensions to include to define the subspace.
#' @param res_name Character vector. Appendix of the result score's name
#' (e.g. method used to combine batches).
#' Needs to have the same length as metrics or NULL.
#' @param k  Numeric. Number of k-nearest neighbours (Knn) to use.
#' @param k_min Numeric. Minimum number of Knn to include
#' (see \code{\link{cms}}). Relevant for metrics: 'cms'.
#' @param smooth Logical. Indicating if cms results should be smoothened within
#' each neighbourhood using the weigthed mean. Relevant for metric: 'cms'.
#' @param cell_min Numeric. Minimum number of cells from each group to be
#' included into the AD test. Should be > 4. Relevant for metric: 'cms'.
#' @param batch_min Numeric. Minimum number of cells per batch to include in to
#' the AD test. If set, neighbours will be included until batch_min cells from
#' each batch are present. Relevant for metrics: 'cms'.
#' @param unbalanced Boolean. If TRUE, neighbourhoods with only one batch present
#' will be set to NA. This way they are not included into any summaries or
#' smoothening. Relevant for metrics: 'cms'.
#' @param weight Boolean. If TRUE, batch probabilities to calculate the isi
#' score are weighted by the mean distance of their cells towards the cell
#' of interest. Relevant for metrics: 'isi'.
#' @param k_pos Numeric. Position of cell to be used as reference within mixing
#' metric. See \code{\link[Seurat]{MixingMetric}} for details.
#' Relevant for metric: 'mixing_metric'
#' @param sce_pre_list A list of \code{SingleCellExperiment} objects with single
#' datasets before integration. Names should correspond to levels in
#' \code{colData(sce_combined)[,group]}. Relevant for metric: 'ldfDiff'
#' @param dim_combined Character. Name of embeddings to use as subspace to
#' calculate LDF after integration. Default is \code{dim_red}.
#' Relevant for metric 'ldfDiff'.
#' @param assay_pre Character. Name of the assay to use for PCA.
#' Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce_pre))}. Default is "logcounts".
#' Relevant for metric 'ldfDiff'.
#' @param n_dim_orig Number of PCs to use in original space.
#' See \code{\link[Seurat]{LocalStruct}} for details.
#' Relevant for metric 'local_structure'.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#' cms scores shall be calculated in parallel. Relevant for metric: 'cms'.
#'
#' @details evalIntegration is a wrapper function for different metrics to
#' understand results of integrated single cell data sets.
#' In general there are metrics evaluationg the *mixing* of datasets, that is,
#' metrics that show whether there still is a bias for different datasets
#' after integration. Furthermore there are metrics to evaluate how well the
#' dataset internal structure has been retained, that is, metrics that show whether
#' there has been (potentially biological) signal removed or noise added by
#' integration.
#'
#' @section Metrics:
#' Here we provide the following metrics:
#' \describe{
#'   \item{cms}{Cellspecific Mixing Score. Metric that tests the hypothesis
#'   that group-specific distance distributions of knn cells have the same
#'   underlying unspecified distribution. The score can be interpreted as the
#'   probability of having unbiased data according to the batch variable
#'   (see \code{\link{cms}}).}
#'   \item{isi}{Inverse Simpson Index. Metric that uses the Inverse Simpson’s
#'   Index to calculate the diversification within a specified
#'   neighbourhood as the probability of sampling a cell from the same batch
#'   twice. For 2 batches a score close to 2 indicates high randomness/mixing
#'   The inverse Simpson index has been proposed as a diversity score for batch
#'   mixing in single cell RNAseq by Korunsky et al. They provide a
#'   distance-based neighborhood weightening in their Lisi package.}
#'   \item{mixing_metric}{Mixing Metric. Metric using the median position of the
#'    kth cell from each batch within its knn as a score. The lower the better
#'    mixed is the neighbourhood (See \code{\link[Seurat]{MixingMetric}}.)}
#'    \item{entropy}{Shannon entropy. Metric calculating the Shannon entropy of
#'    the batch/group variable within each cell's k-nearest neigbours.
#'    For balanced batches the entropy is closer to 1 the higher the variables
#'    randomness. For unbalanced batches entropy should only be used as a
#'    relative metric in a comparative setting (See \code{\link{entropy}}.)}
#'    \item{ldfDiff}{Local density factor differences. Metric that determines
#'    cell-specific changes in the Local Density Factor before and after data
#'    integration. A metric/difference close to 0 indicates no distortion of
#'    the previous structure (see \code{\link{ldfDiff}}).}
#'    \item{local_structure}{Local structure. Metric that compares the
#'    intersection of knn from the same batch before and after integration
#'    returning the average between all groups. The higher the more neighbours
#'    were reproduced after integration (See \code{\link[Seurat]{LocalStruct}}.
#'    )}
#' }
#'
#' @return A \code{SingleCellExperiment} with the chosen metric's score within
#' colData.
#' @export
#'
#' @references
#' Korsunsky I Fan J Slowikowski K Zhang F Wei K et. al. (2018).
#' Fast, sensitive, and accurate integration of single cell data with Harmony.
#' bioRxiv (preprint).
#'
#' Stuart T Butler A Hoffman P Hafemeister C Papalexi E et. al. (2019)
#' Comprehensive Integration of Single-Cell Data.
#' Cell.
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:15, 300:320, 16:30)]
#' sce_batch1 <- sce[,colData(sce)$batch == "1"]
#' sce_batch2 <- sce[,colData(sce)$batch == "2"]
#' pre <- list("1" = sce_batch1, "2" = sce_batch2)
#'
#' sce <- evalIntegration(metrics = c("cms", "mixing_metric", "isi", "entropy"), sce, "batch", k = 20)
#' sce <- evalIntegration("ldfDiff", sce, "batch", k = 20, sce_pre_list = pre)
#'
#' @importFrom scater runPCA
#' @importFrom SingleCellExperiment reducedDims colData counts reducedDims<-
#' @importFrom Seurat as.Seurat LocalStruct MixingMetric
#' @importFrom magrittr %>% set_names
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData<- assays assay assay<-
evalIntegration <- function(metrics, sce, group, dim_red = "PCA",
                            assay_name = "logcounts", n_dim = 10,
                            res_name = NULL, k = NULL, k_min = NULL,
                            smooth = NULL, cell_min = NULL, batch_min = NULL,
                            unbalanced = FALSE, weight  = TRUE, k_pos = NULL,
                            sce_pre_list = NULL, dim_combined = NULL,
                            assay_pre = NULL, n_dim_orig = NULL,
                            BPPARAM=SerialParam()){
    #------------------- Check input parameter----------------------
    metric_params <- c("cms", "ldfDiff", "isi", "mixing_metric",
                       "local_structure", "entropy")
    if( !all(metrics %in% metric_params) ){
        stop("Error: 'metrics' is unknown. Please define one or more of 'cms', 'isi',
             'ldfDiff', 'mixing_metric', 'local_structure', 'entropy'")
    }
    if( !is(sce, "SingleCellExperiment" )){
        stop("Error: 'sce' must be a 'SingleCellExperiment' object.")
    }
    if( !group %in% names(colData(sce)) ){
        stop("Error: 'group' variable must be in 'colData(sce)'")
    }
    if( !is(colData(sce)[, group], "factor") ){
        sce[[group]] <- as.factor(colData(sce)[, group])
    }

    #------------------overall parameter settings ------------------------
    if( !is.null(res_name) ){
        if( length(res_name) != length(metrics) ){
            stop("Error: Define 'res_name' for all metrics to calculate'")
        }
        names(res_name) <- metrics
    }
    default <- ifelse(is.null(k), TRUE, FALSE)

    #------------------ run metrics -------------------------------
    #------------------ mixing-metrics ----------------------------
    if( "cms" %in% metrics ){
        #Set default parameter
        if( is.null(k_min) ){
            k_min <- NA
        }
        if( is.null(cell_min) ){
            cell_min <- 10
        }
        if( is.null(smooth) ){
            smooth <- TRUE
        }
        if( is.null(res_name) ){
            res_name[["cms"]] <- NULL
        }
        #Check parameter
        if( is.null(k) ){
            stop("Please specify 'k', the number of nearest neigbours to check
                 for equal mixing, e.g. median of cells/celltype.")
        }
        #run cms
        sce <- cms(sce, k = k, group = group, dim_red = dim_red,
                   assay_name = assay_name, res_name = res_name[["cms"]],
                   k_min = k_min, batch_min = batch_min,
                   unbalanced = unbalanced, smooth = smooth, n_dim = n_dim,
                   cell_min = cell_min, BPPARAM = BPPARAM)
        }

    if( "isi" %in% metrics ){
        #Set default parameter
        if( is.null(k) ){
            stop("Please specify 'k', the number of nearest neigbours to check
                 for equal mixing, e.g. median of cells/celltype.")
        }
        if( is.null(res_name) ){
            res_name[["isi"]] <- NULL
        }

        #run isi
        sce <- isi(sce, k = k, group = group, dim_red = dim_red,
                   assay_name = assay_name, weight = weight,
                   res_name = res_name[["isi"]], n_dim = n_dim)
    }

    if( "mixing_metric" %in% metrics ){
        #Set default parameter
        if( is.null(k) | default ){
            k <- 300
        }
        if( is.null(k_pos) ){
            k_pos <- 5
        }
        #Check parameter
        if( k >= ncol(sce) ){
            warning("'k' exceeds number of cells. Is set to max (all cells).")
            k <- ncol(sce) - 1
        }
        if( n_dim  > ncol(reducedDims(sce)[[dim_red]]) ){
            warning("'n_dim' exceeds number of provided reduced dimensions.
                    Is set to max (all dims).")
            n_dim <- ncol(reducedDims(sce)[[dim_red]])
        }
        #Run mixing metric
        reducedDims(sce)[[dim_red]] <- .defineSubspace(sce, assay_name, dim_red,
                                                       n_dim)
        #Make sure it can be converted into seurat
        if( !"logcounts" %in% names(assays(sce)) ){
            assay(sce, "logcounts") <- log2(counts(sce) + 1)
        }
        if( is.null(rownames(sce)) ){
            rownames(sce) <- paste0("gene", seq_len(nrow(sce)))
        }

        seurat <- as.Seurat(sce)
        mix_dist <- MixingMetric(seurat, grouping.var = group,
                                 reduction = dim_red,
                                 dims = seq_len(n_dim),
                                 k = k_pos, max.k = k)
        if( !is.null(res_name) ){
            colData(sce)[, paste0("mm.", res_name[["mixing_metric"]])] <- mix_dist
        }else{
            colData(sce)[, "mm"] <- mix_dist
        }
    }

    if( "entropy" %in% metrics ){
        #Check parameter
        if( is.null(k) | default ){
            stop("Please specify 'k', the number of nearest neigbours to check
                 for equal mixing, e.g. median of cells/celltype.")
        }
        if( is.null(res_name) ){
            res_name[["entropy"]] <- NULL
        }
        #run entropy
        sce <- entropy(sce, group = group, k = k, dim_red = dim_red,
                       assay_name = assay_name, n_dim = n_dim,
                       res_name = res_name[["entropy"]])
    }

    #------------------ structure metrics ----------------------------
    if( "ldfDiff" %in% metrics ){
        #Set default parameter
        if( is.null(k) | default ){
            stop("Please specify 'k', the number of nearest neigbours to check
                 for structual changes")
        }
        if( is.null(dim_combined) ){
            dim_combined <- dim_red
        }
        if( is.null(assay_name) ){
            assay_name <- "logcounts"
        }
        if( is.null(assay_pre) ){
            assay_pre <- "logcounts"
        }
        if( is.null(res_name) ){
            res_name[["ldfDiff"]] <- NULL
        }
        #run ldfDiff
        sce <- ldfDiff(sce_pre_list = sce_pre_list, sce_combined = sce, group,
                       k = k, dim_red = dim_red, dim_combined = dim_combined,
                       assay_pre = assay_pre, assay_combined = assay_name,
                       n_dim = n_dim, res_name = res_name[["ldfDiff"]])
        }

    if( "local_structure" %in% metrics ){
        #Set default parameter
        if( is.null(k) | default ){
            k <- 100
        }
        if( is.null(n_dim_orig) ){
            n_dim_orig <- 10
        }
        #Check parameter
        if( k > min(table(droplevels(colData(sce)[, group]))) ){
            warning("'k' exceeds number of cells/batch.
                    Is set to min(cells/batch).")
            k <- min(table(droplevels(colData(sce)[, group])))
        }
        if( n_dim  > ncol(reducedDims(sce)[[dim_red]]) ){
            warning("'n_dim' exceeds number of provided reduced dimensions.
                    Is set to max (all dims).")
            n_dim <- ncol(reducedDims(sce)[[dim_red]])
        }
        #Run localStructure
        reducedDims(sce)[[dim_red]] <- .defineSubspace(sce, assay_name, dim_red,
                                                       n_dim)
        #Make sure it can be converted into seurat
        if( !"logcounts" %in% names(assays(sce)) ){
            assay(sce, "logcounts") <- log2(counts(sce) + 1)
        }
        if( is.null(rownames(sce)) ){
            rownames(sce) <- paste0("gene", seq_len(nrow(sce)))
        }
        seurat <- as.Seurat(sce)
        local_struct <- LocalStruct(seurat, grouping.var = group, neighbors = k,
                                    reduction = dim_red,
                                    reduced.dims = seq_len(n_dim),
                                    orig.dims = seq_len(n_dim_orig))
        #Transform batch-sorted list into vector with cells ordered by sce
        cells_by_batch <- colnames(sce)[order(colData(sce)[,group])]
        local_s_batch <- local_struct %>% unlist %>%
            set_names(cells_by_batch)
        local_s_sce <- local_s_batch[colnames(sce)]

        if( !is.null(res_name) ){
            colData(sce)[,paste0("localStruct.", res_name[["local_structure"]])] <-
                local_s_sce
        }else{
            colData(sce)[,"localStruct"] <- local_s_sce
        }
    }

    return(sce)
}

