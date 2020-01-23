# Helper and internal functions for cms

#' .cmsCell
#'
#' Function to calculate a cellspecific mixing score (cms) of groups/batches.
#'
#' @param cell Character. Name of the cell to calculate cms for.
#' Needs to be one of \code{rownames(knn)}.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(knn)}.
#' @param knn List with three elements. First "index" with indices of knn cells.
#' Second "distance" with distances to knn cells. Third a slot named by
#' \code{group} variable with group level of knn cells.
#' @param k_min Numeric. Minimum number of knn to include.
#' Default is NA (see Details).
#' @param cell_min Numeric. Minimum number of cells from each group to be
#' included into the AD test. Should be > 4 to make 'ad.test' working.
#' @param batch_min Numeric. Minimum number of cells per batch to include in to
#' the AD test. If set neighbours will be included until batch_min cells from
#' each batch are present.
#' @param unbalanced Boolean. If True neighbourhoods with only one batch present
#' will be set to NA. This way they are not included into any summaries or
#' smoothening.
#' @param sce A \code{SingleCellExperiment} object with the combined data.
#'
#' @details The cms function tests the hypothesis, that group-specific distance
#' distributions of knn cells have the same underlying unspecified distribution.
#' It performs Anderson-Darling tests as implemented in the
#' \code{kSamples package}.
#' In default the function uses all distances and group label defined in knn.
#' If \code{k_min} is specified, the first local minimum of the overall distance
#' distribution with at least kmin cells is used. This can be used to adapt to
#' the local structure of the datatset e.g. prevent cells from a distinct
#' different cluster to be included.
#'
#' @seealso \code{\link{ad.test}}, \code{\link{cms}}, \code{\link{.smoothCms}}
#'
#' @family helper functions
#'
#' @return A p.value as resulting from the ad.test.
#'
#' @importFrom kSamples ad.test
#' @importFrom magrittr set_colnames %>% extract2
#' @importFrom dplyr mutate mutate_at group_by_at filter summarize n
.cmsCell <- function(cell, group, knn, k_min = NA, batch_min = NULL,
                     cell_min = 4, unbalanced = FALSE, sce){
    #get knn distances and group assignments
    knn_cell <- cbind(knn[["distance"]][cell, ], knn[[group]][cell, ]) %>%
        as.data.frame %>%
        set_colnames(c("distance", group)) %>%
        mutate(distance = as.numeric(as.character(.data$distance))) %>%
        mutate_at(group, factor)

    #Filter cells within a distinct different density distribution
    if( !is.na(k_min) ){
        knn_cell <- .filterLocMin(knn_cell, k_min)
    }

    if( !is.null(batch_min) ){
        knn_cell <- .filterKnn(knn_cell, batch_min, group = group, sce = sce)
    }


    #filter groups with to few cells (cell_min, default 4)
    groups_included <- knn_cell %>% group_by_at(group) %>%
        summarize("n_group" = n()) %>%
        filter(.data$n_group > cell_min) %>%
        extract2(group) %>% droplevels() %>% levels()

    #do not perform AD test if only one group with enough cells is in knn.
    if( length(groups_included) <= 1 ){
        p <- ifelse(unbalanced, NA, 0)
    }else{
        dist_included <- lapply(groups_included, function(group_level){
            dist <- knn_cell$distance[which(knn_cell[, group] %in% group_level)]
        })
        names(dist_included) <- groups_included
        #filter cells with the same representation
        if( any(map(dist_included, sum) == 0) ){
            warning("Distances equal to 0 - cells with identical
                    representations detected. NA assigned!")
            p <- NA
        }else{
            #perform AD test with remaining cell
            k_samp <- ad.test(dist_included)
            p <- mean(k_samp[["ad"]][, " asympt. P-value"])
        }
    }
    p
}


#' .filterLocMin
#'
#' Function to filter knn by overall distance density distribution.
#'
#' @param knn_cell Data frame with one column "distance" and one column named
#' by the group variable. Rows correspond to the knn cells and do not need
#' rownames.
#' @param k_min Numeric. Minimum number of Knn to include.
#'
#' @details  Internal function to filter cells used for cms testing to come
#' from a continous overall density distribution function
#' (similar to cluster definitions). 'filterLocMin' is only applied, if k-min
#' is specified as parameter in \code{\link{.cmsCell}} or \code{\link{cms}}.
#'
#' @seealso \code{\link{.cmsCell}}
#' @family helper functions
#' @return data.frame with two columns (index, distance) for filtered knn cells.
#'
#' @importFrom stats density
.filterLocMin <- function(knn_cell, k_min){
    distances <- density(knn_cell$distance)$x
    dist_density <- density(knn_cell$distance)$y
    #Find local minima
    loc_min <- distances[which(diff(sign(diff(dist_density)))== 2)+1]
    # Find first local minimum (> minimal threshold of cells (k_min))
    all_locmin <- unlist(lapply(loc_min, function(minim){
        ind <- length(which(knn_cell$distance <= minim))
    }))
    loc_th <- suppressWarnings(min(which(all_locmin > k_min)))
    loc_dist_th <- loc_min[loc_th]

    #filter cells before first local minimum (>k_min)
    if(length(loc_min) > 0 & !is.na(loc_dist_th)){
        knn_cell <- knn_cell[which(knn_cell$distance <= loc_dist_th),]
    }
    knn_cell
}

#' .filterKnn
#'
#' @param knn_cell Data frame with one column "distance" and one column named
#' by the group variable. Rows correspond to the knn cells and do not need
#' rownames.
#' @param batch_min Numeric. Minimum number of cells per batch to include.
#' @param group Character. Name of group/batch variable.
#' Needs to be one of \code{names(knn)}.
#' @param sce A \code{SingleCellExperiment} object with the combined data.
#'
#' @seealso \code{\link{.cmsCell}}
#' @family helper functions
#' @return data.frame with two columns (index, distance) for filtered knn cells.
#'
.filterKnn <- function(knn_cell, batch_min, group, sce){
    # Make sure at least batch_min cells are within knn
    b_ids <- levels(colData(sce)[, group])
    batch_k <- rep(batch_min, length(b_ids)) %>% set_names(b_ids)
    min_indices <- which(table(knn_cell[,group]) >= batch_min)
    max_indices <- which(table(knn_cell[,group]) < batch_min)
    if( length(max_indices) > 0 ){
        max_cell <- table(knn_cell[,group])[[max_indices]]
        batch_k[max_indices] <- max_cell
        warning(paste0("There are less than 'batch_min' cells of each batch in
                a reasonable sized neighbourhood. ", max_cell," number of cells
                       are used"))
    }
    #Find maximum of ith cells within knn
    batch_ind <- lapply(b_ids, function(batch){
        k_ind <- which(knn_cell[,group] %in% batch)[batch_k[batch]]
    }) %>% unlist()

    knn_cell <- knn_cell[seq_len(max(batch_ind, na.rm = TRUE)),]
}

#' .smoothCms
#'
#' Performs weighted smoothening of cms scores
#'
#' @param knn List with three elements. First "index" with indices of knn cells.
#' Second "distance" with distances to knn cells. Third a slot named by
#' \code{group} variable with group level of knn cells.
#' @param cms_raw Matrix with raw cms scores for all cells specified in
#' \code{cell_names} and \code{knn}. Colnames need to be "cms.
#' @param cell_names Character vector with cell names corresponding to the
#' rownames of the list elements in \code{knn} and \code{rownames(cms_raw)}.
#' @param k_min Numeric. Minimum number of knn to include.
#' Default is NA (see Details).
#' @param k Numeric. Number of k-nearest neighbours (knn) to use.
#'
#' @details Internal function to smooth cms scores. In a complete random setting
#'  cms scores are uniform distributed. To reduce the resulting random variance
#'  and enable visualization of local pattern cms scores can be smoothened
#'  assuming that within one region mixing is uniform. Generates smoothened cms
#'  scores using weigthed means of cms scores within the k-nearest
#'  neighbourhood. Reciprocal distances are used as weights.
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
.smoothCms <- function(knn, cms_raw, cell_names, k_min, k){

    # cms assignment of knn cells for each cell
    knn[["cms"]] <- cell_names %>%
        map(function(cell_id) cms_raw[knn[["index"]][cell_id, ], "cms"]) %>%
        bind_rows() %>% t()


    # how many cells to smooth over
    k_smooth <- ifelse(!is.na(k_min), k_min, k)

    # calculate weigthed mean
    cms_smooth <- cell_names %>%
        map(function(cell_id){
            #add 1 to for distances < 1
            weights <- c(1, (1/(knn[[2]][cell_id, seq_len(k_smooth)] + 1)))
            knn_cms <- c(cms_raw[cell_id, "cms"],
                         knn[["cms"]][cell_id, seq_len(k_smooth)])
            cms_new <- weighted.mean(knn_cms, weights, na.rm = TRUE)
        }) %>% bind_rows() %>% t() %>%
        set_colnames("cms_smooth")
    res_cms <- cbind(cms_smooth, cms_raw)
    na_cells <- which(is.na(res_cms[,"cms"]))
    res_cms[na_cells, "cms_smooth"] <- NA
    res_cms
}

## DefineSubspace
### See ldf_helper.R
