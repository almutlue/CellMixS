library(SingleCellExperiment)
#get simulated scRNA seq data with 3 unbalanced batches
load(system.file("extdata/sim30.rda", package = "CellMixS"))
sce <- sim_30[[1]][, c(1:50,500:550)]

## Tests for cms fuction:
### Include internal functions as cmsCell, filterLocMin and cmsSmooth.

test_that("test that output of cms is correct",{
  cms_smooth <- cms(sce, k = 20, group = "batch",
                    dim_red = "TSNE", assay_name = "counts")
  cms_kmin <- cms(sce, k = 20, group = "batch", kmin = 15)
  cms_raw <- cms(sce, k = 20, group = "batch", smooth = FALSE)
  cms_default <- cms(sce, k = 20, group = "batch")
  cms_only <- as.matrix(cms_default[,"cms"])
  colnames(cms_only) <- "cms"

 expect_is(cms_smooth, "matrix")
  expect_equal(ncol(cms_smooth), 2)
  expect_is(cms_raw, "matrix")
  expect_is(cms_kmin, "matrix")
  expect_identical(cms_only, cms_raw)
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
})


