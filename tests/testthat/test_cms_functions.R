library(SingleCellExperiment)
library(CellMixS)
#get simulated scRNA seq data with 3 unbalanced batches
sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
sce <- sim_list[[1]][, c(1:30,300:320)]

## Tests for cms fuction:
### Include internal functions as cmsCell, filterLocMin and cmsSmooth.

test_that("test that output of cms is correct",{
    sce_cms_smooth <- cms(sce, k = 20, group = "batch",
                      dim_red = "TSNE", assay_name = "counts", n_dim = 2)
    sce_cms_kmin <- cms(sce, k = 20, group = "batch", k_min = 15,
                        res_name = "kmin", n_dim = 2)
    sce_unbalanced <- cms(sce, k = 30, group = "batch", unbalanced = TRUE)
    sce_batch_min <- cms(sce, k = 30, group = "batch", batch_min = 5)
    sce_df <- as.data.frame(assay(sce))

    expect_is(sce_cms_smooth, "SingleCellExperiment")
    expect_equal(ncol(colData(sce_cms_smooth)) - ncol(colData(sce)), 2)
    expect_is(sce_cms_kmin, "SingleCellExperiment")
    expect_is(sce_cms_kmin$cms_smooth.kmin, "numeric")
    expect_length(which(is.na(sce_unbalanced$cms)), 4)
    expect_equal(length(which(is.na(sce_unbalanced$cms))),
                 length(which(is.na(sce_unbalanced$cms_smooth))))
    expect_error(cms(sce = sce, k=20, group = "batch",
                     dim_red = "TSNE", assay_name = "raw_counts"),
                 "Ambigious parameter: Specify subspace parameter.
         * For precalculated embeddings keep 'assay_name' as default.
         * For PCA based on 'assay_name' keep 'dim_red' as default.",
                 fixed = TRUE)
    expect_error(cms(sce = sce, k=20, group = "batch",
                     dim_red = "pca", assay_name = "raw_counts"),
                 "Parameter 'assay_name' not found: Provide a valid value.",
                 fixed = TRUE)
    expect_warning(cms(sce = sce, k=20, group = "batch",
                       dim_red = "tsne", assay_name = "counts"),
                   "'dim_red' not found:
            PCA subspace is used to calculate distances.", fixed = TRUE)
    expect_error(cms(sce = sce_df, k=20, group = "batch",
                     dim_red = "pca", assay_name = "raw_counts"),
                 "Error: 'sce' must be a 'SingleCellExperiment' object.",
                 fixed = TRUE)
    expect_error(cms(sce = sce, k=20, group = "batch2",
                     dim_red = "pca", assay_name = "raw_counts"),
                 "Error: 'group' variable must be in 'colData(sce)'",
                 fixed = TRUE)
    expect_warning(cms(sce = sce, k=100, group = "batch",
                     dim_red = "pca", assay_name = "counts"),
                 "'k' exceeds number of cells. Is set to max (all cells).",
                 fixed = TRUE)
    expect_error(cms(sce = sce, k=20, group = "batch", k_min = 10,
                     dim_red = "PCA", assay_name = "counts", batch_min = 10),
                 "Error: 'k_min' and 'batch_min'were set. Cms needs to be calculated
             based on neighbourhood density or batch occurence, not both.",
                 fixed = TRUE)
})


