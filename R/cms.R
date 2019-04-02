# Function to calculate cms

## Main function of the package!

#' cms
#'
#' Calculates cell-specific mixing scores based on euclidean distances within a subspace of integrated data.
#'
#' @param sce A \code{SingleCellExperiment} object with the combined data.
#' @param k Numeric. Number of k-nearest neighbours (Knn) to use.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(colData(sce))}
#' @param dim_red Character. Name of embeddings to use as subspace for distance distributions.
#' Default is "PCA".
#' @param assay_name Character. Name of the assay to use for PCA.
#' Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce))}. Default is "logcounts".
#' @param res_name Character. Appendix of the result score's name (e.g. method used to combine batches).
#' @param kmin Numeric. Minimum number of Knn to include.
#' Default is NA (see Details).
#' @param smooth Logical. Indicating if cms results should be smoothened within each neighbourhood using the weigthed mean.
#' @param n_dim Numeric. Number of dimensions to include to define the subspace.
#' @param cell_min Numeric. Minimum number of cells from each group to be included into the AD test.
#' Should be > 10 to make the ad.test function working.
#'
#' @details The cms function tests the hypothesis, that group-specific distance distributions of knn cells have the same underlying unspecified distribution.
#' It performs Anderson-Darling tests as implemented in the \code{kSamples package}.
#' In default the function uses all distances and group label defined in knn.
#' If \code{kmin} is specified, the first local minimum of the overall distance distribution with at least kmin cells is used.
#' This can be used to adapt to the local structure of the datatset e.g. prevent cells from a distinct different cluster to be included.
#' If 'dim_red' is not defined or default cms will calculate a PCA using \code{runPCA}.
#' Results will be appended to \code{colData(sce)}.
#' Names can be specified using \code{res_name}.
#'
#' @family cms functions
#' @seealso \code{\link{.cmsCell}}, \code{\link{.smoothCms}}.
#'
#' @return A matrix with cells as rows and cms (and cms_smooth) as columns.
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
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' sce <- sim_30[[1]][, c(1:50)]
#'
#' cms_smooth <- cms(sce, k = 20, group = "batch")
#' cms_raw <- cms(sce, k = 20, group = "batch", smooth = FALSE)
#'
#' @importFrom scater runPCA
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom SummarizedExperiment assays
#' @importFrom BiocNeighbors findKNN
#' @importFrom magrittr %>% set_rownames set_colnames
#' @importFrom purrr %>% map
#' @importFrom listarrays set_dimnames
#' @importFrom dplyr bind_rows
cms <- function(sce, k, group, dim_red = "PCA", assay_name = "logcounts",
                res_name = NULL, kmin = NA, smooth = TRUE, n_dim = 20,
                cell_min = 10){

    #------------------Check input parameter ---------------------------------#
    if(cell_min < 10){
        stop("Error: 'cell_min' is < 10. Must be > 10 to estimate cms.")
    }
    if(!is(sce, "SingleCellExperiment")){
        stop("Input error: class('sce') must be 'SingleCellExperiment'.")
    }

    #check group variable class
    if(!is(colData(sce)[,group], "factor")){
        sce[[group]] <- as.factor(colData(sce)[, group])
    }

    cell_names <- colnames(sce)
    names(cell_names) <- cell_names

    #---------------------------------------------------------------------------#

    #----------------- determine knn matrix -----------------------------------#
    subspace <- .defineSubspace(sce, assay_name, dim_red, n_dim)
    #determine knn
    knn <- findKNN(subspace, k=k) %>% map(set_rownames, cell_names)

    # group assignment of knn cells for each cell
    knn[[group]] <- matrix(
        colData(sce)[as.numeric(as.character(knn$index)), group],
        nrow = nrow(knn$index)) %>%
        set_rownames(cell_names)
    #---------------------------------------------------------------------------#

    #----------------- calculate cms score  -----------------------------------#

    cms_raw <- cell_names %>%
        map(.cmsCell, group = group, knn = knn, kmin=kmin, cell_min = cell_min) %>%
        bind_rows() %>% t() %>%
        set_colnames("cms")

    if(isTRUE(smooth)){
        res_cms <- .smoothCms(knn, cms_raw, cell_names, kmin, k)
    }else{
        res_cms <- cms_raw
    }

    #Add to colData of sce
    if(!is.null(res_name)){
        res_cms <- res_cms %>% set_colnames(paste0(colnames(.), "_", res_name))
    }
    colData(sce) <- cbind(colData(sce), res_cms)
    sce
}


