# Functions to calculate cms

#' cms_cell
#'
#' Function to calculate a cellspecific mixing score (cms) of groups/batches.
#'
#' @param cell Character string defining the name of the cell to calculate cms for. Should be one of the rownames of knn.
#' @param group Character string specifying the variable name used to define groups (batches). Should have a corresonding element in knn.
#' @param knn List containing distances and group information of cells. Should have one slot type double named "distance" with distances towards the k nearest neighbours of the cell and one slot type double named after the group variable.
#' Both slots need to be aligned and the specified cell's name need to correspond to one of the rownames, while the columns represent the ordered knn.
#' @param kmin Numeric giving the minimum number of k nearest neighbours to use. cells to include.
#' If k_min is defined the overall distance density distribution of knn will be determined and only cells up to the first local minima that includes more than kmin cells will be used.
#' @param cell_min Minimum number of cells from each group to be included into the AD test for calculation of cms.
#' Should be > 10 to make the \code{\link{ad.test}} working.
#'
#' @details The cms tests the hypothesis, that group-specific distance distributions of the k-nearest neighbouring cells have the same underlying unspecified distribution using the Anderson-Darling test.
#' In default the function uses all distances and group label defined in knn. If k_min is specified, the overall distance distribution is checked for the first local minimum, that includes at least kmin cells.
#' This can be used to adapt to the local structure of the datatset e.g. prevent cells from a distinct different cluster to be included, while keeping the overall number of k nearest neighbours reasonable large.
#' This function is mainly intended as core function for a set of functions to calculate cms for a single cell containing object: \code{\link{cms}}.
#'
#' @seealso \code{\link{ad.test}} for Anderson-Darling test on k-samples, \code{\link{cms}}, \code{\link{smooth_cms}}
#' @family cms functions
#' @return
#' @export
#'
#' @examples
#' @importFrom kSamples ad.test
cms_cell <- function(cell, group, knn, kmin = NA, cell_min = 10){
  #get knn distances and group assignments
  knn_cell <- cbind(knn[["distance"]][cell,], knn[[group]][cell,])
  knn_cell <- as.data.frame(knn_cell)
  colnames(knn_cell) <- c("distance", group)
  knn_cell[,group] <- as.factor(knn_cell[,group])

  #Filter cells within a distinct different density distribution (only if local_min = TRUE (default))
  if(!is.na(kmin)){
    knn_cell <- filter_locmin(knn_cell, kmin)
  }

  #filter groups with to few cells (cell_min, default 10)
  groups_included <- levels(knn_cell[,group])[which(table(knn_cell[,group]) > cell_min)]

  #do not perform AD test if only one group with enough cells is in neighbourhood
  if(length(groups_included) <= 1){
    p <- 0
  }else{
    dist_included <- lapply(groups_included, function(group_level){
      dist <- knn_cell$distance[which(knn_cell[,group] %in% group_level)]
    })
    names(dist_included) <- groups_included
    #perform AD test with remaining cells
    k_samp <- ad.test(dist_included)
    p <- mean(k_samp[["ad"]][," asympt. P-value"])
  }
  p
}

### Helper for cms.cell

#' filter.locmin
#'
#' Function to filter number of knn to use for cms calculation by overall distance density distribution.
#'
#' @param knn_cell Data frame with distances and group information of cells. Should have one column named "distance" with distances towards the k nearest neighbours and one column named after the group variable.
#' Rows correspond to the knn cells and do not need rownames.
#' @param kmin Numeric giving the minimum number of k nearest neighbours to use. cells to include.
#' If k_min is defined the overall distance density distribution of knn will be determined and only cells up to the first local minima that includes more than kmin cells will be used.
#'
#' @details  Internal function to filter cells used for cms testing to come from a continous overall density distribution function (similar to cluster definitions).
#' filter.locmin is only applied, if k-min is specified as parameter in \code{\link{cms_cell}} or \code{\link{cms}}.
#'
#' @seealso \code{\link{cms_cell}} for cms calculation
#' @family cms functions
#' @return
#'
#' @examples
#'
#' @importFrom stats density
filter_locmin <- function(knn_cell, kmin){
  distances <- density(knn_cell$distance)$x
  dist_density <- density(knn_cell$distance)$y
  #Find local minima
  loc_min <- distances[which(diff(sign(diff(dist_density)))== 2)+1]
  # Find first local minimum (that is larger than a minimal threshold of cells (kmin))
  all_locmin <- unlist(lapply(loc_min, function(minim){ind <- length(which(knn_cell$distance <= minim))}))
  loc_th <- suppressWarnings(min(which(all_locmin > kmin)))                #indice of first local minimum (>kmin)
  loc_dist_th <- loc_min[loc_th]                          #up to which distance are cells included

  #filter cells before first local minimum (>kmin)
  if(length(loc_min) > 0 & !is.na(loc_dist_th)){
    knn_cell <- knn_cell[which(knn_cell$distance <= loc_dist_th),]
  }
  knn_cell
}

### Calculate cms

#' cms
#'
#' Calculates cell-specific mixing scores based on euclidean distances within subspace of integrated data.
#'
#' @param sce SingleCellExperiment \code{\link{SingleCellExperiment}} object containing all cells of interest.
#' @param k Numeric, defining the number of k-nearest neighbours to use for testing the "mixing" within the neighbourhood.
#' @param group Character string specifying the name of variable used to define groups (batches). Should be have a corresonding element in knn.
#' @param embedding Character string specifying the embeddings to use as subspace to calculate distance distributions.
#' The default PCA performs a principle component analysis and uses the n_pc first PCs to define the subspace.
#' @param assay_name Character string defining the count assay to use. If not otherwise specified by 'embeddings', a PCA on this assay is performed to calculate a new PCA subspace. Must be one of 'names(assays(sce))'.
#' Default is 'logcounts' and should not be changed, if other embeddings than PCA should be used.
#' @param kmin Numeric giving the minimum number of k nearest neighbours to use. cells to include.
#' If k_min is defined the overall distance density distribution of knn will be determined and only cells up to the first local minima that includes more than kmin cells will be used.
#' @param smooth A logical indicating if cms results should be smoothened within each neighbourhood.
#' @param n_pc Numeric giving the number of PCs to include to define the subspace to calculate cell distances in.
#' @param cell_min Minimum number of cells from each group to be included into the AD test for calculation of cms.
#' Should be > 10 to make the ad.test function working.
#'
#' @details The cms tests the hypothesis, that group-specific distance distributions of the k-nearest neighbouring cells have the same underlying unspecified distribution using the Anderson-Darling test as implemented in \code{\link{cms.cell}}.
#' In default the function uses all distances and group label defined in knn. f *kmin* is specified, the first local minimum of the overall distance distribution with at least kmin cells is used.
#' This can be used to adapt to the local structure of the datatset e.g. prevent cells from a distinct different cluster to be included, while keeping the overall number of k nearest neighbours reasonable large.
#' cms.pca uses a *SingleCellExperiment* \code{\link{SingleCellExperiment}} object as input and calculates euclidean distances within the subspace defined by the n_pc argument (default = 30).
#' For *smoothening* weigthed means of cms scores within each cell's k-nearest neighbourhood are calculated.
#' Reciprocal distances are used as weights.
#'
#' @family cms functions
#' @seealso \code{\link{cms_cell}}, \code{\link{smooth_cms}}.
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
cms <- function(sce, k, group, embedding = "PCA", assay_name = "logcounts", kmin = NA, smooth = TRUE, n_pc = 30, cell_min = 10){

  #------------------Check input parameter ---------------------------------#
  if(cell_min < 10){
    stop("Error: 'cell_min' is < 10. Must be > 10 to estimate cms.")
  }

  #check group variable class
  colData(sce)[,group] <- as.factor(colData(sce)[,group])
  cell_names <- colnames(sce)

  # Check assay and embedding
  if(!assay_name %in% c("logcounts", "counts") & !embedding %in% c("pca", "PCA", "Pca")){
    stop("Ambigious parameter: Please specify parameter for distance calculations.
         * If precalculated embeddings shall be used, keep 'assay_name' as default.
         * If a PCA based on 'assay_name' shall be used, keep 'embedding' as default.")
  }
  if(is.null(reducedDim(sce, embedding))){
    if(!assay_name %in% names(assays(sce))){
      stop("Parameter 'assay_name' not found: Please provide a valid value.")
    }

    warning("Embedding not found: PCA subspace is used to calculate distances.")

    #run PCA, if PCs do not exist
    sce <- runPCA(sce, ncomponents = n_pc, exprs_values = assay_name)
    embedding <- "PCA"

  }else if(!assay_name %in% c("logcounts", "counts")){
    sce <- runPCA(sce, ncomponents = n_pc, exprs_values = assay_name)
    embedding <- "PCA"
  }
  #---------------------------------------------------------------------------#

  #----------------- determine knn matrix -----------------------------------#
  subspace <- reducedDim(sce, embedding)
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

  cms_raw <- do.call(rbind, lapply(cell_names, cms_cell, group = group, knn = knn, kmin=kmin, cell_min = cell_min))
  rownames(cms_raw) <- cell_names
  colnames(cms_raw) <- "cms"

  if(smooth == TRUE){
    res_cms <- smooth_cms(knn, cms_raw, cell_names, kmin, k)
  }else{
    res_cms <- cms_raw
  }
  res_cms
}



## Helpers for cms function


#' smooth_cms
#'
#' Performs weighted smoothening of cms scores
#'
#' @param knn List containing distances and batch information of cells. Should have one slot type double named "distance" with distances towards the k nearest neighbours of the cell and one slot type double named after the group variable specified in 'batch'.
#' Both slots need to be aligned and the specified cell's name need to correspond to one of the rownames, while the columns represent the ordered knn.
#' @param cms_raw Matrix containing raw cms scores for all cells specified in cell_names and knn. Rownames need to correspond to cellnames and colnames should be "cms".
#' @param cell_names Character vector with cell names corresponding to the rownames of the list elements in the knn object and the cms_raw object.
#' @param kmin Numeric giving the minimum number of k nearest neighbours to use.
#' If k_min is defined the overall distance density distribution of knn will be determined and only cells up to the first local minima that includes more than kmin cells will be used.
#' @param k Numeric, defining the number of k-nearest neighbours to use for testing the "mixing" within the neighbourhood.
#'
#' @details Internal function to smoothes cms scores. In a complete random setting cms scores are uniform distributed.
#' To reduce the resulting random variance and enable visualization of local tendencies and pattern cms scores can be smoothened assuming that within one region mixing is uniform.
#' Generates smoothened cms scores using weigthed means of cms scores within the k-nearest neighbourhood.
#' Reciprocal distances are used as weights. Returns a matrix with "smoothened_cms" and original "cms" values.
#'
#' @family cms functions
#' @seealso \code{\link{cms_cell}}, \code{\link{cms}}
#'
#' @return
#'
#' @examples
#'
#' @importFrom stats weighted.mean
smooth_cms <- function(knn, cms_raw, cell_names, kmin, k){

  # cms assignment of knn cells for each cell
  knn[["cms"]] <- do.call(rbind, lapply(cell_names, function(cell_id){
    cms_raw[knn[["indices"]][cell_id,], "cms"]}))
  rownames(knn[["cms"]]) <- cell_names

  # how many cells to smooth over
  k_smooth <- ifelse(!is.na(kmin), kmin, k)

  # calculate weigthed mean
  cms_smooth <- do.call(rbind, lapply(cell_names, function(cell_id){
    #add 1 to ensure that distances < 1 are less weighted than the original cell
    weights <- c(1,(1/(knn[[2]][cell_id,c(1:k_smooth)] + 1)))
    knn_cms <- c(cms_raw[cell_id,"cms"],knn[["cms"]][cell_id, c(1:k_smooth)])
    cms_new <- weighted.mean(knn_cms, weights)
  }))

  res_cms <- cbind(cms_smooth, cms_raw)
  rownames(res_cms) <- cell_names
  colnames(res_cms) <- c("cms_smooth", "cms")
  res_cms
}


### summarize function

#' cms_summary
#'
#' Summarizes cms scores, if specified based on groups
#'
#' @param cms_res data frame of cms scores to be summarized. Columns should correspond to cms scores, rows to cells.
#' @param sum_var variable to group summaries on. Need to correspond to a colData variable
#' @param sce SingleCellExperiment object describing the cells specified in cms_res.
#'
#' @details Returns the mean of cms score to make different conditions/methods/.. comparable.
#' In default the mean of each cloumn of cms_res is returned. Groups to summarize can be specified by sum_var and sce.
#'
#' @family cms functions
#' @seealso \code{\link{compare_integration}}, \code{\link{compare_cluster}}
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom dplyr as_tibble group_by_ summarize_all funs
#' @importFrom SingleCellExperiment colData
#' @importFrom magrittr %>%
cms_summary <- function(cms_res, sum_var = NULL, sce = NULL){
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
