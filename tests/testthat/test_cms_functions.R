library(SingleCellExperiment)
data(sim_30)
sce <- sim_30[[1]][, c(1:100,500:600)]

## Tests for cms... fuctions:
### Include internal functions as cms.cell, filter.locmin and cms.smooth.

test_that("test that output of cms.pca is correct",{
  cms_smooth <- cms.pca(sce = sce, k=20, group = "batch")
  cms_raw <- cms.pca(sce = sce, k=20, group = "batch", smooth = FALSE)
  cms_only <- as.matrix(cms_smooth[,"cms"])
  colnames(cms_only) <- "cms"

  expect_is(cms_smooth, "matrix")
  expect_is(cms_raw, "matrix")
  expect_identical(cms_raw, cms_only)
  expect_error(cms.pca(sce = sce, k=20, group = "batch", cell_min = 5),
               "Error: 'cell_min' is < 10. Must be > 10 to estimate cms.")
})



test_that("test that output of cms.integration is correct",{
  cms_smooth <- cms.integration(sce, k = 20, group = "batch",
                    embedding = "TSNE", assay_name = "counts")
  cms_kmin <- cms.integration(sce, k = 20, group = "batch", k_min = 15)
  cms_raw <- cms.integration(sce, k = 20, group = "batch", smooth = FALSE)
  cms_default <- cms.integration(sce, k = 20, group = "batch")
  cms_pca <- cms.pca(sce = sce, k=20, group = "batch")

 expect_is(cms_smooth, "matrix")
  expect_equal(ncol(cms_smooth), 2)
  expect_is(cms_raw, "matrix")
  expect_is(cms_kmin, "matrix")
  expect_identical(cms_default, cms_pca)
  expect_error(cms.integration(sce = sce, k=20, group = "batch",
                              embedding = "TSNE", assay_name = "raw_counts"),
          "Ambigious parameter: Please specify parameter for distance calculations.
         * If precalculated embeddings shall be used, keep 'assay_name' as default.
         * If a PCA based on 'assay_name' shall be used, keep 'embedding' as default.", fixed = TRUE)
})

test_that("test output of cms summary",{
  cms_res <- cms.integration(sce, k = 20, group = "batch")
  cms_raw <- cms.integration(sce, k = 20, group = "batch", smooth = FALSE)
  summary_res <- cms.summary(cms_res)
  summary_raw <- cms.summary(cms_raw)
  summary_batch <- cms.summary(cms_res, sum_var = "batch", sce = sce)
  summary_raw_batch <- cms.summary(cms_res, sum_var = "batch", sce = sce)

  expect_is(summary_res, "data.frame")
  expect_is(summary_raw, "data.frame")
  expect_equal(nrow(summary_batch), length(table(colData(sce)$batch)))
  expect_equal(nrow(summary_batch), nrow(summary_raw_batch))
  expect_error(cms.summary(cms_res, sum_var = "batch"),
               "Missing variable: Please provide a 'sce' object.")
  expect_error(cms.summary(cms_raw, sum_var = "quatsch", sce = sce),
               "Missing variable: Could not find 'sum_var'.
           Please specify one of names(colData(sce)).", fixed = TRUE)
})
