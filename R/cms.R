# Function to calculate cms

## Main function of the package!

#' cms
#'
#' Calculates cell-specific mixing scores based on euclidean distances within a
#' subspace of integrated data.
#'
#' @param sce A \code{SingleCellExperiment} object with the combined data.
#' @param k Numeric. Number of k-nearest neighbours (Knn) to use.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}
#' @param dim_red Character. Name of embeddings to use as subspace for distance
#' distributions. Default is "PCA".
#' @param assay_name Character. Name of the assay to use for PCA.
#' Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce))}. Default is "logcounts".
#' @param res_name Character. Appendix of the result score's name
#' (e.g. method used to combine batches).
#' @param k_min Numeric. Minimum number of Knn to include.
#' Default is NA (see Details).
#' @param smooth Logical. Indicating if cms results should be smoothened within
#' each neighbourhood using the weigthed mean.
#' @param n_dim Numeric. Number of dimensions to include to define the subspace.
#' @param cell_min Numeric. Minimum number of cells from each group to be
#' included into the AD test.
#' Should be > 4 to make the ad.test function working.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#' cms scores shall be calculated in parallel.
#'
#' @details The cms function tests the hypothesis, that group-specific distance
#' distributions of knn cells have the same underlying unspecified distribution.
#' It performs Anderson-Darling tests as implemented in the
#' \code{kSamples package}.
#' In default the function uses all distances and group label defined in knn.
#' If \code{k_min} is specified, the first local minimum of the overall distance
#' distribution with at least k_min cells is used. This can be used to adapt to
#' the local structure of the datatset e.g. prevent cells from a
#' different cluster to be included. If 'dim_red' is not defined or default cms
#' will calculate a PCA using \code{runPCA}. Results will be appended to
#' \code{colData(sce)}. Names can be specified using \code{res_name}.
#' If multiple cores are available cms scores can be calculated in parallel
#' (does not work on Windows). Parallelization can be specified using BPPARAM.
#'
#' @family cms functions
#' @seealso \code{\link{.cmsCell}}, \code{\link{.smoothCms}}.
#'
#' @return A \code{SingleCellExperiment} with cms (and cms_smooth) within
#' colData.
#'
#' @references
#' Scholz, F. W. and Stephens, M. A. (1987).
#' K-Sample Anderson-Darling Tests.
#' J. Am. Stat. Assoc.
#'
#'
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
#' sce <- sim_list[[1]][, c(1:50)]
#'
#' sce_cms <- cms(sce, k = 20, group = "batch", n_dim = 2)
#'
#' @importFrom scater runPCA
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom SummarizedExperiment assays
#' @importFrom BiocNeighbors findKNN
#' @importFrom magrittr %>% set_rownames set_colnames
#' @importFrom purrr %>% map
#' @importFrom listarrays set_dimnames
#' @importFrom dplyr bind_rows .data
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData<-
#' @importFrom BiocParallel bplapply SerialParam
cms <- function(sce, k, group, dim_red = "PCA", assay_name = "logcounts",
                res_name = NULL, k_min = NA, smooth = TRUE, n_dim = 20,
                cell_min = 4, BPPARAM=SerialParam()){
    #------------------Check input parameter ---------------------------------#
    if(cell_min < 4){
        stop("Error: 'cell_min' is < 4. Must be > 4 to estimate cms.")
    }
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
    knn <- findKNN(subspace, k=k) %>% map(set_rownames, cell_names)

    # group assignment of knn cells for each cell
    knn[[group]] <- matrix(
        colData(sce)[as.numeric(as.character(knn$index)), group],
        nrow = nrow(knn$index)) %>%
        set_rownames(cell_names)

    #----------------- calculate cms score  ----------------------------------#
    cms_raw <- bplapply(cell_names,.cmsCell, group = group, knn = knn,
                        k_min=k_min, cell_min = cell_min, BPPARAM = BPPARAM) %>%
        bind_rows() %>% t() %>%
        set_colnames("cms")

    if(isTRUE(smooth)){
        res_cms <- .smoothCms(knn, cms_raw, cell_names, k_min, k)
    }else{
        res_cms <- cms_raw
    }

    #Add to colData of sce
    if(!is.null(res_name)){
        res_cms <- res_cms %>% set_colnames(paste0(colnames(.), ".", res_name))
    }
    colData(sce) <- cbind(colData(sce), res_cms)
    sce
}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("."))
