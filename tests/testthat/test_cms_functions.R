library(SingleCellExperiment)
library(CellMixS)
#get simulated scRNA seq data with 3 unbalanced batches
sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
sce <- sim_list[[1]][, c(1:50,500:550)]

## Tests for cms fuction:
### Include internal functions as cmsCell, filterLocMin and cmsSmooth.

test_that("test that output of cms is correct",{
    sce_cms_smooth <- cms(sce, k = 20, group = "batch",
                      dim_red = "TSNE", assay_name = "counts", n_dim = 2)
    sce_cms_kmin <- cms(sce, k = 20, group = "batch", k_min = 15, res_name = "kmin")
    sce_cms_raw <- cms(sce, k = 20, group = "batch", smooth = FALSE)
    sce_cms_default <- cms(sce, k = 20, group = "batch")
    cms_only <- as.matrix(sce_cms_default$cms)
    sce_df <- as.data.frame(assay(sce))

    expect_is(sce_cms_smooth, "SingleCellExperiment")
    expect_equal(ncol(colData(sce_cms_smooth)) - ncol(colData(sce)), 2)
    expect_is(sce_cms_raw, "SingleCellExperiment")
    expect_equal(ncol(colData(sce_cms_raw)) - ncol(colData(sce)), 1)
    expect_is(sce_cms_kmin, "SingleCellExperiment")
    expect_is(sce_cms_kmin$cms_smooth.kmin, "numeric")
    expect_identical(cms_only, as.matrix(sce_cms_raw$cms))
    expect_error(cms(sce = sce, k=20, group = "batch",
                     dim_red = "TSNE", assay_name = "raw_counts"),
                 "Ambigious parameter: Please specify parameter for distance calculations.
         * If precalculated embeddings shall be used, keep 'assay_name' as default.
         * If a PCA based on 'assay_name' shall be used, keep 'dim_red' as default.", fixed = TRUE)
    expect_error(cms(sce = sce, k=20, group = "batch",
                     dim_red = "pca", assay_name = "raw_counts"),
                 "Parameter 'assay_name' not found: Please provide a valid value.",
                 fixed = TRUE)
    expect_warning(cms(sce = sce, k=20, group = "batch",
                       dim_red = "tsne", assay_name = "counts"),
                   "'dim_red' not found: PCA subspace is used to calculate distances.",
                   fixed = TRUE)
    expect_error(cms(sce = sce_df, k=20, group = "batch",
                     dim_red = "pca", assay_name = "raw_counts"),
                 "Input error: class('sce') must be 'SingleCellExperiment'.",
                 fixed = TRUE)
})


