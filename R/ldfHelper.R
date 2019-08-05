# LDF helper functions


### LDF function modified from the DDoutlier package

#' .ldfKnn
#'
#' Calculates the Local Density Factor as implemented in the \code{DDoutlier}
#' package with a predefined KNN neighbourhood.
#'
#' @param dataset Matrix with cell embeddings with cells as rows and reduced
#' dimensions as cloumns. Subspace to determine LDF in.
#' @param knn_object List with k-nearest neighbours (KNN) as provided by
#' \code{get.knn} from the \code{FNN} package.
#' First element named "indices" contains indices of KNN in \code{dataset}.
#' Second element named "distance" contains distances of KNN in \code{dataset}.
#' Third element named "cell_name" contains rownames of KNN in \code{dataset}.
#' @param k Numeric. Number of KNN used. Should correspond to \code{knn_object}.
#' @param h Numeric. Bandwidth  for  kernel  functions.
#' The  greater  the  bandwidth, the smoother kernels and lesser weight are put
#' on outliers. Default is 1
#' @param c Scaling constant for comparison of LDE to neighboring observations.
#' Default is 1.
#'
#' @details LDF fuction modified from the \code{DDoutlier} package.
#' Calculates a Local Density Estimate (LDE) and Local Density Factor (LDF) with
#'  a gaussian kernel. Modified to use a predefined KNN neighbourhood.
#' For \code{\link{ldfSce}} this is essential to determine LDF after data
#' integration on the same set of cells.
#'
#' @return List with two elements "LDE" and "LDF".
#'
#' @family helper functions
#' @seealso \code{\link{ldfSce}}
#'
.ldfKnn <- function(dataset, knn_object, k=k, h=1, c=1){

  dim <- ncol(dataset)

  LDE <- do.call(rbind, lapply(rownames(dataset), function(cell_nam){
    k_dist <- knn_object[["distance"]][knn_object[["cell_name"]][cell_nam,], k]
    dist_po <- knn_object[["distance"]][cell_nam,]
    reach_dist <- apply(cbind(k_dist, dist_po), 1, max)

    kernel <- 1/((((2*pi)^(dim/2)))*((h*k_dist)^dim)) *
        exp(-((reach_dist^2)/(2*((h*k_dist)^2))))

    LDE <- (1/k)*sum(kernel)
  }))

  rownames(LDE) <- rownames(dataset)
  colnames(LDE) <- "LDE"

  LDF <- do.call(rbind, lapply(rownames(dataset), function(cell_nam){
    knn_nam <- knn_object[["cell_name"]][cell_nam,]
    LDF <- sum(LDE[knn_nam,]/k)/(LDE[cell_nam,]+(c*sum(LDE[knn_nam,]/k)))
    LDF <- round(LDF, digits = 3)
  }))

  rownames(LDF) <- rownames(dataset)
  colnames(LDF) <- "LDF"

  return_list <- list(LDE=LDE, LDF=LDF)
  return(return_list)

}

## Define subspace

#' .defineSubspace
#'
#' Helper function for ldfSce and cms to define or recalculate the subspace for
#' analysis.
#'
#' @param sce A \code{SingleCellExperiment} object with the data to define the
#' subspace.
#' @param assay_name Character. Name of the assay to use for PCA. Only relevant
#' if no existing 'dim_red' is provided. Must be one of
#' \code{names(assays(sce))}.
#' @param dim_red Character. Name of embeddings to use as subspace.
#' @param n_dim Numeric. Number of subspace elements to include to define
#' subspace.
#'
#' @details Function to determine the subspace for \code{ldfDiff} and
#' \code{cms}. Checks whether the defined 'dim_red' is present.
#' Only if no subspace is defined or present it will perform a PCA using
#' \code{runPCA}. To calculate PCA counts defined in 'assay_name' are used.
#'
#' @family helper functions
#' @seealso \code{\link{ldfSce}}, \code{\link{cms}}.
#'
#' @return A matrix of cell embeddings with reduced dimensions as columns.
#'
#' @importFrom scater runPCA
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom SummarizedExperiment assays
.defineSubspace <- function(sce, assay_name, dim_red, n_dim){

  # Check assay and dim_red
  if( !assay_name %in% c("logcounts", "counts") &
      !dim_red %in% c("pca", "PCA", "Pca") ){
    stop("Ambigious parameter: Specify subspace parameter.
         * For precalculated embeddings keep 'assay_name' as default.
         * For PCA based on 'assay_name' keep 'dim_red' as default.")
  }
  if(!assay_name %in% c("logcounts", "counts") & dim_red %in% c("PCA")){
    sce <- runPCA(sce, ncomponents = n_dim, exprs_values = assay_name)
    dim_red <- "PCA"
  }
  if(!dim_red %in% reducedDimNames(sce)){
    if(!assay_name %in% names(assays(sce))){
      stop("Parameter 'assay_name' not found: Provide a valid value.")
    }

    warning("'dim_red' not found:
            PCA subspace is used to calculate distances.")

    #run PCA, if PCs do not exist
    sce <- runPCA(sce, ncomponents = n_dim, exprs_values = assay_name)
    dim_red <- "PCA"

  }else if(!assay_name %in% c("logcounts", "counts")){
    sce_pre <- runPCA(sce, ncomponents = n_dim, exprs_values = assay_name)
    dim_red <- "PCA"
  }
  # check number dimension to use
  if(n_dim > ncol(reducedDim(sce, dim_red))){
    stop("Parameter 'n_dim' is greater than reduced dimensional space:
         Please provide a valid value.")
  }
  subspace <- reducedDim(sce, dim_red)[,seq_len(n_dim)]

}
