#' Toolbox to explore batch effects and data integration in scRNA data.
#'
#' \pkg{CellMixS} provides metrics and functions to evaluate batch effects,
#' data integration and batch effect correction in single cell trancriptome
#' data with single cell resolution. Results can be visualized and summarised
#' on different levels, e.g. on cell, celltype or dataset level.
#'
#' In particular, \pkg{CellMixS} includes two main metrics:
#' Cellspecific mixing scores to determine the probability of random mixing
#' in each cell's neighbourhood. It can be assesed via the \code{\link{cms}} function.
#' Local Density Factor Differences to evaluate the effect of data
#' integration methods on batch internal structures. It can be assesed via the
#' \code{\link{ldfDiff}} function.
#'
#' @import kSamples
#' @author Almut Lütge \email{almut.luetge@@uzh.ch}
#' @author Mark D Robinson \email{mark.robinson@@imls.uzh.ch}
#' @name CellMixS-package
#' @docType package
NULL
