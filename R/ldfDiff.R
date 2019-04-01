# Function to calculate LDF change

## Function to evaluate data integration

#' ldfDiff
#'
#' Determines cell-specific changes in the Local Density Factor before and after data integration.
#'
#' @param sce_pre_list A list of \code{SingleCellExperiment} objects with single datasets before integration.
#' Names should correspond to levels in \code{colData(sce_combined)$group}
#' @param sce_combined A \code{SingleCellExperiment} object with the combined data.
#' @param group Character. Name of group/batch variable that separates elements of \code{sce_pre_list}.
#' Needs to be one of \code{names(colData(sce_combined))}.
#' @param k Numeric. Number of k-nearest neighbours (Knn) to use.
#' @param dim_red Character. Name of embeddings to use as subspace to calculate LDF before integration. Default is "PCA".
#' @param dim_combined Character. Name of embeddings to use as subspace to calculate LDF after integration. Default is \code{dim_red}.
#' @param assay_pre Character. Name of the assay to use for PCA. Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce_pre))}. Default is "logcounts".
#' @param assay_combined Character. Name of the assay to use for PCA. Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce_combined))}. Default is "logcounts".
#' @param n_dim Numeric. Number of PCs to include to define subspaces.
#'
#' @details The ldfDiff function calculates differences in LDF for each element in \code{sce_pre_list} and their corresponding cells in \code{sce_combined} using \code{\link{ldfSce}}.
#' If 'dim_red' is not defined a PCA will be calculated using \code{runPCA}.
#' In this case 'assay_pre' need to refer to the data slot that shall define the subspace.
#' Similar refer 'dim-combined' and 'assay_combined' to the integrated subspace or to the resp. "corrected" count data slot.
#' 'k' can be used to define the level of local structure that is tested.
#' The smaller 'k' the more focus is on detailed structures, while a large k will tets overall changes.
#'
#'
#' @family ldf functions
#' @seealso \code{\link{ldfSce}}, \code{\link{.ldfKnn}}.
#'
#' @return A data.frame with cells as rows and difference in LDF as column named "diff".
#'
#' @references
#' Latecki, Longin Jan and Lazarevic, Aleksandar and Pokrajac, Dragoljub (2007).
#' Outlier Detection with Kernel Density Functions.
#' Mach. Learn. Data Min. Pattern Recognit..
#' Springer Berlin Heidelberg.
#'
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' sce <- sim_30[[4]][, c(1:50, 500:550)]
#' sce_batch1 <- sce[,colData(sce)$batch == "1"]
#' sce_batch2 <- sce[,colData(sce)$batch == "2"]
#' sce_pre_list <- list("1" = sce_batch1, "2" = sce_batch2)
#'
#' ldf_diff <- ldfDiff(sce_pre_list, sce, k = 10, group = "batch", dim_combined = "MNN", n_dim = 5)
#'
#'
ldfDiff <- function(sce_pre_list, sce_combined, group, k = 75, dim_red = "PCA", dim_combined = dim_red, assay_pre = "logcounts", assay_combined = "logcounts", n_dim = 20){
  diff <- do.call(rbind, lapply(names(sce_pre_list),
                                ldfSce, sce_pre_list = sce_pre_list,
                                sce_combined = sce_combined,
                                group = group,
                                k = k,
                                dim_red = dim_red,
                                dim_combined = dim_combined,
                                assay_pre = assay_pre,
                                assay_combined = assay_combined,
                                n_dim = n_dim))
}



#' ldfSce
#'
#' Determines cell-specific changes in the Local Density Factor before and after data integration for one specific group.
#'
#' @param sce_name Character. Name of the element in \code{sce_pre_list} to calculate LDF differences in.
#' @param sce_pre_list A list of \code{SingleCellExperiment} objects with single datasets before integration.
#' Names should correspond to levels in \code{colData(sce_combined)$group} and \code{sce_name}.
#' @param sce_combined A \code{SingleCellExperiment} object with the combined data.
#' @param group Character. Name of group/batch variable that separates elements of \code{sce_pre_list}.
#' Needs to be one of \code{names(colData(sce_combined))}.
#' @param k Numeric. Number of k-nearest neighbours (Knn) to use.
#' @param dim_red Character. Name of embeddings to use as subspace to calculate LDF before integration. Default is "PCA".
#' @param dim_combined Character. Name of embeddings to use as subspace to calculate LDF after integration. Default is \code{dim_red}.
#' @param assay_pre Character. Name of the assay to use for PCA. Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce_pre))}. Default is "logcounts".
#' @param assay_combined Character. Name of the assay to use for PCA. Only relevant if no existing 'dim_red' is provided.
#' Must be one of \code{names(assays(sce_combined))}. Default is "logcounts".
#' @param n_dim Numeric. Number of PCs to include to define subspaces.
#'
#' @details The ldfSce function calculates differences in LDF for one specified element in \code{sce_pre_list} and their corresponding cells in \code{sce_combined}.
#' If 'dim_red' is not defined a PCA will be calculated using \code{runPCA}.
#' In this case 'assay_pre' need to refer to the data slot that shall define the subspace.
#' Similar refer 'dim-combined' and 'assay_combined' to the integrated subspace or to the resp. "corrected" count data slot.
#' 'k' can be used to define the level of local structure that is tested.
#' The smaller 'k' the more focus is on detailed structures, while a large k will tets overall changes.
#' K-nearest neighbours (KNN) are determined in the subspaces before integration defined by 'dim_red'.
#' The same set of KNN are used to determine LDF before and after integration.
#'
#'
#' @family ldf functions
#' @seealso \code{\link{ldfDiff}}, \code{\link{.ldfKnn}}.
#'
#' @return A data.frame with cells as rows and difference in LDF as column named "diff".
#'
#' @references
#' Latecki, Longin Jan and Lazarevic, Aleksandar and Pokrajac, Dragoljub (2007).
#' Outlier Detection with Kernel Density Functions.
#' Mach. Learn. Data Min. Pattern Recognit..
#' Springer Berlin Heidelberg.
#'
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' load(system.file("extdata/sim30.rda", package = "CellMixS"))
#' sce <- sim_30[[4]][, c(1:50, 500:550)]
#' sce_batch1 <- sce[,colData(sce)$batch == "1"]
#' sce_pre_list <- list("1" = sce_batch1)
#'
#' ldf_1 <- ldfSce("1", sce_pre_list, sce, k = 10, group = "batch", dim_combined = "MNN", n_dim = 5)
#'
#' @importFrom stats dist
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom SummarizedExperiment assays
#' @importFrom  BiocNeighbors findKNN
ldfSce <-function(sce_name, sce_pre_list, sce_combined, group, k = 75, dim_red = "PCA", dim_combined = dim_red, assay_pre = "logcounts", assay_combined = "logcounts", n_dim = 20){
  #sce before integration
  sce_pre <- sce_pre_list[[sce_name]]
  cell_names <- colnames(sce_pre)

  #sub sce after integration
  sce_post <- sce_combined[, colData(sce_combined)[,group] == sce_name]

  #determine subspace
  subspace <- .defineSubspace(sce_pre, assay_pre, dim_red, n_dim)

  #----------------- determine ldf pre-----------------------------------#
  # test size of k has
  if(k >= nrow(subspace)){
    stop("Parameter 'k' is greater than dataset size: Please provide a valid value.")
  }
  knn <- findKNN(subspace, k=k)
  rownames(knn[[1]]) <- rownames(subspace)  #index of knn cells per cell
  rownames(knn[[2]]) <- rownames(subspace) #euclidean dist of knn cells per cell

  #assign names to cell index
  knn[["cell_name"]] <- do.call(rbind, lapply(rownames(subspace), function(cell_id){
    rownames(subspace)[knn[["index"]][cell_id,]]
  }))

  rownames(knn[["cell_name"]]) <- rownames(subspace)
  ldf_pre <- .ldfKnn(subspace, knn_object = knn, k = k, c = 0.5)

  #---------------------------------------------------------------------------#

  #----------------- determine ldf post-----------------------------------#
  # euclidean distances in integrated subspace
  knn_int <- knn
  subspace_int <- .defineSubspace(sce_combined, assay_combined, dim_combined, n_dim)
  subspace_int <- subspace_int[rownames(knn_int$index),]
  # calculate new distances (keeping neighbours)
  knn_int$distance <- do.call(rbind, lapply(rownames(knn_int$index), function(cell){
    knn_cells <- knn_int[["cell_name"]][cell,]
    subspace_sub <- subspace_int[c(cell,knn_cells),]
    distance <- as.matrix(dist(subspace_sub))[cell,-1]
  }))

  rownames(knn_int[["distance"]]) <- rownames(knn_int$index)

  ldf_post <- .ldfKnn(subspace_int, knn_object = knn_int, k = k, c = 0.5)

  #----------------------------------------------------------------------------#
  #calculate differences in LDF
  diff <- data.frame("diff" = ldf_post[["LDF"]][,"LDF"] - ldf_pre[["LDF"]][,"LDF"],
                     row.names = rownames(ldf_pre[["LDF"]]))
  }
