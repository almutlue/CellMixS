# Helper and internal functions for cms

#' .cmsCell
#'
#' Function to calculate a cellspecific mixing score (cms) of groups/batches.
#'
#' @param cell Character. Name of the cell to calculate cms for. Needs to be one of \code{rownames(knn)}.
#' @param group Character. Name of group/batch variable. Needs to be one of \code{names(knn)}.
#' @param knn List with three elements. First "index" with indices of KNN cells.
#' Second "distance" with distances to KNN cells. Third a slot named by \code{group} variable with group level of KNN cells.
#' @param kmin Numeric. Minimum number of Knn to include. Default is NA (see Details).
#' @param cell_min Numeric. Minimum number of cells from each group to be included into the AD test.
#' Should be > 10 to make the ad.test function working.
#'
#' @details The cms function tests the hypothesis, that group-specific distance distributions of knn cells have the same underlying unspecified distribution.
#' It performs Anderson-Darling tests as implemented in the \code{kSamples package}.
#' In default the function uses all distances and group label defined in knn.
#' If \code{kmin} is specified, the first local minimum of the overall distance distribution with at least kmin cells is used.
#' This can be used to adapt to the local structure of the datatset e.g. prevent cells from a distinct different cluster to be included.
#'
#' @seealso \code{\link{ad.test}} for Anderson-Darling test on k-samples, \code{\link{cms}}, \code{\link{.smoothCms}}
#'
#' @family helper functions
#'
#' @return A p.value as resulting from the ad.test.
#'
#' @importFrom kSamples ad.test
#' @importFrom magrittr set_colnames %>% extract2
#' @importFrom dplyr mutate mutate_at group_by_at filter summarize
.cmsCell <- function(cell, group, knn, kmin = NA, cell_min = 10){
    #get knn distances and group assignments
    knn_cell <- cbind(knn[["distance"]][cell, ], knn[[group]][cell, ]) %>%
        as.data.frame %>%
        set_colnames(c("distance", group)) %>%
        mutate(distance = as.numeric(as.character(distance))) %>%
        mutate_at(group, factor)

    #Filter cells within a distinct different density distribution (only if local_min = TRUE (default))
    if( !is.na(kmin) ){
        knn_cell <- .filterLocMin(knn_cell, kmin)
    }

    #filter groups with to few cells (cell_min, default 10)
    groups_included <- knn_cell %>% group_by_at(group) %>%
        summarize("n_group" = n()) %>%
        filter(n_group > cell_min) %>%
        extract2(group) %>% levels()

    #do not perform AD test if only one group with enough cells is in neighbourhood
    if (length(groups_included) <= 1) {
        p <- 0
    } else {
        dist_included <- lapply(groups_included, function(group_level){
            dist <- knn_cell$distance[which(knn_cell[, group] %in% group_level)]
        })
        names(dist_included) <- groups_included
        #perform AD test with remaining cells
        k_samp <- ad.test(dist_included)
        p <- mean(k_samp[["ad"]][, " asympt. P-value"])
    }
    p
}


#' .filterLocMin
#'
#' Function to filter knn by overall distance density distribution.
#'
#' @param knn_cell Data frame with one column "distance" and one column named by the group variable.
#' Rows correspond to the knn cells and do not need rownames.
#' @param kmin Numeric. Minimum number of Knn to include.
#'
#' @details  Internal function to filter cells used for cms testing to come from a continous overall density distribution function (similar to cluster definitions).
#' 'filterLocMin' is only applied, if k-min is specified as parameter in \code{\link{.cmsCell}} or \code{\link{cms}}.
#'
#' @seealso \code{\link{.cmsCell}}
#' @family helper functions
#' @return data.frame with two columns (index, distance) for filtered knn cells.
#'
#' @importFrom stats density
.filterLocMin <- function(knn_cell, kmin){
    distances <- density(knn_cell$distance)$x
    dist_density <- density(knn_cell$distance)$y
    #Find local minima
    loc_min <- distances[which(diff(sign(diff(dist_density)))== 2)+1]
    # Find first local minimum (that is larger than a minimal threshold of cells (kmin))
    all_locmin <- unlist(lapply(loc_min, function(minim){
        ind <- length(which(knn_cell$distance <= minim))
    }))
    loc_th <- suppressWarnings(min(which(all_locmin > kmin)))
    loc_dist_th <- loc_min[loc_th]

    #filter cells before first local minimum (>kmin)
    if(length(loc_min) > 0 & !is.na(loc_dist_th)){
        knn_cell <- knn_cell[which(knn_cell$distance <= loc_dist_th),]
    }
    knn_cell
}

#' .smoothCms
#'
#' Performs weighted smoothening of cms scores
#'
#' @param knn List with three elements. First "index" with indices of KNN cells.
#' Second "distance" with distances to KNN cells. Third a slot named by \code{group} variable with group level of KNN cells.
#' @param cms_raw Matrix with raw cms scores for all cells specified in \code{cell_names} and \code{knn}. Colnames need to be "cms.
#' @param cell_names Character vector with cell names corresponding to the rownames of the list elements in \code{knn} and \code{rownames(cms_raw)}.
#' @param kmin Numeric. Minimum number of Knn to include. Default is NA (see Details).
#' @param k Numeric. Number of k-nearest neighbours (Knn) to use.
#'
#' @details Internal function to smooth cms scores. In a complete random setting cms scores are uniform distributed.
#' To reduce the resulting random variance and enable visualization of local pattern cms scores can be smoothened assuming that within one region mixing is uniform.
#' Generates smoothened cms scores using weigthed means of cms scores within the k-nearest neighbourhood.
#' Reciprocal distances are used as weights.
#'
#' @family helper functions
#' @seealso \code{\link{.cmsCell}}, \code{\link{cms}}
#'
#' @return matrix with two columns ("cms_smooth", "cms").
#'
#'
#' @importFrom stats weighted.mean
#' @importFrom magrittr %>% set_colnames
#' @importFrom purrr %>% map
#' @importFrom dplyr bind_rows
.smoothCms <- function(knn, cms_raw, cell_names, kmin, k){

    # cms assignment of knn cells for each cell
    knn[["cms"]] <- cell_names %>%
        map(function(cell_id) cms_raw[knn[["index"]][cell_id, ], "cms"]) %>%
        bind_rows() %>% t()


    # how many cells to smooth over
    k_smooth <- ifelse(!is.na(kmin), kmin, k)

    # calculate weigthed mean
    cms_smooth <- cell_names %>%
        map(function(cell_id){
            #add 1 to ensure that distances < 1 are less weighted than the original cell
            weights <- c(1, (1/(knn[[2]][cell_id, seq_len(k_smooth)] + 1)))
            knn_cms <- c(cms_raw[cell_id, "cms"], knn[["cms"]][cell_id, seq_len(k_smooth)])
            cms_new <- weighted.mean(knn_cms, weights)
        }) %>% bind_rows() %>% t() %>%
        set_colnames("cms_smooth")

    res_cms <- cbind(cms_smooth, cms_raw)
}

## DefineSubspace
### See ldf_helper.R
