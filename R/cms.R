# Function to calculate cms

## Main function of the package!

#' cms
#'
#' Calculates cell-specific mixing scores based on euclidean distances within subspace of integrated data.
#'
#' @param sce SingleCellExperiment \code{\link{SingleCellExperiment}} object containing all cells of interest.
#' @param k Numeric, defining the number of k-nearest neighbours to use for testing the "mixing" within the neighbourhood.
#' @param group Character string specifying the name of variable used to define groups (batches). Should be have a corresonding element in knn.
#' @param dim_red Character string specifying the embeddings to use as subspace to calculate distance distributions.
#' The default PCA performs a principle component analysis and uses the n_pc first PCs to define the subspace.
#' @param assay_name Character string defining the count assay to use. If not otherwise specified by 'dim_red', a PCA on this assay is performed to calculate a new PCA subspace. Must be one of 'names(assays(sce))'.
#' Default is 'logcounts' and should not be changed, if other embeddings than PCA should be used.
#' @param kmin Numeric giving the minimum number of k nearest neighbours to use. cells to include.
#' If k_min is defined the overall distance density distribution of knn will be determined and only cells up to the first local minima that includes more than kmin cells will be used.
#' @param smooth A logical indicating if cms results should be smoothened within each neighbourhood.
#' @param n_pc Numeric giving the number of PCs to include to define the subspace to calculate cell distances in.
#' @param cell_min Minimum number of cells from each group to be included into the AD test for calculation of cms.
#' Should be > 10 to make the ad.test function working.
#'
#' @details The cms tests the hypothesis, that group-specific distance distributions of the k-nearest neighbouring cells have the same underlying unspecified distribution using the Anderson-Darling test as implemented in \code{\link{cmsCell}}.
#' In default the function uses all distances and group label defined in knn. f *kmin* is specified, the first local minimum of the overall distance distribution with at least kmin cells is used.
#' This can be used to adapt to the local structure of the datatset e.g. prevent cells from a distinct different cluster to be included, while keeping the overall number of k nearest neighbours reasonable large.
#' cms.pca uses a *SingleCellExperiment* \code{\link{SingleCellExperiment}} object as input and calculates euclidean distances within the subspace defined by the n_pc argument (default = 30).
#' For *smoothening* weigthed means of cms scores within each cell's k-nearest neighbourhood are calculated.
#' Reciprocal distances are used as weights.
#'
#' @family cms functions
#' @seealso \code{\link{cmsCell}}, \code{\link{smoothCms}}.
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom scater runPCA
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom SummarizedExperiment assay
#' @importFrom FNN get.knn
cms <- function(sce, k, group, dim_red = "PCA", assay_name = "logcounts", kmin = NA, smooth = TRUE, n_pc = 30, cell_min = 10){

  #------------------Check input parameter ---------------------------------#
  if(cell_min < 10){
    stop("Error: 'cell_min' is < 10. Must be > 10 to estimate cms.")
  }

  #check group variable class
  colData(sce)[,group] <- as.factor(colData(sce)[,group])
  cell_names <- colnames(sce)

  # Check assay and dim_red
  if(!assay_name %in% c("logcounts", "counts") & !dim_red %in% c("pca", "PCA", "Pca")){
    stop("Ambigious parameter: Please specify parameter for distance calculations.
         * If precalculated embeddings shall be used, keep 'assay_name' as default.
         * If a PCA based on 'assay_name' shall be used, keep 'dim_red' as default.")
  }
  if(is.null(reducedDim(sce, dim_red))){
    if(!assay_name %in% names(assays(sce))){
      stop("Parameter 'assay_name' not found: Please provide a valid value.")
    }

    warning("'dim_red' not found: PCA subspace is used to calculate distances.")

    #run PCA, if PCs do not exist
    sce <- runPCA(sce, ncomponents = n_pc, exprs_values = assay_name)
    dim_red <- "PCA"

  }else if(!assay_name %in% c("logcounts", "counts")){
    sce <- runPCA(sce, ncomponents = n_pc, exprs_values = assay_name)
    dim_red <- "PCA"
  }
  #---------------------------------------------------------------------------#

  #----------------- determine knn matrix -----------------------------------#
  subspace <- reducedDim(sce, dim_red)
  #determine knn
  knn <- get.knn(subspace, k=k, algorithm = 'cover_tree')
  rownames(knn[[1]]) <- cell_names  #indices of knn cells per cell
  rownames(knn[[2]]) <- cell_names  #euclidean dist of knn cells per cell
  names(knn) <- c("indices", "distance")

  # group assignment of knn cells for each cell
  knn[[group]] <- do.call(rbind, lapply(cell_names, function(cell_id){
    colData(sce)[knn[["indices"]][cell_id,], group]
    }))
  rownames(knn[[group]]) <- cell_names
  #---------------------------------------------------------------------------#

  #----------------- calculate cms score  -----------------------------------#

  cms_raw <- do.call(rbind, lapply(cell_names, cmsCell, group = group, knn = knn, kmin=kmin, cell_min = cell_min))
  rownames(cms_raw) <- cell_names
  colnames(cms_raw) <- "cms"

  if(smooth == TRUE){
    res_cms <- smoothCms(knn, cms_raw, cell_names, kmin, k)
  }else{
    res_cms <- cms_raw
  }
  res_cms
}


