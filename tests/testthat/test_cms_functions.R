library(SingleCellExperiment)
data(sim_30)
sce <- sim_30[[1]][, c(1:100,500:600)]

## Tests for cms fuction:
### Include internal functions as cms.cell, filter.locmin and cms.smooth.

test_that("test that output of cms is correct",{
  cms_smooth <- cms(sce, k = 20, group = "batch",
                    embedding = "TSNE", assay_name = "counts")
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
                              embedding = "TSNE", assay_name = "raw_counts"),
          "Ambigious parameter: Please specify parameter for distance calculations.
         * If precalculated embeddings shall be used, keep 'assay_name' as default.
         * If a PCA based on 'assay_name' shall be used, keep 'embedding' as default.", fixed = TRUE)
  expect_error(cms(sce = sce, k=20, group = "batch",
                   embedding = "pca", assay_name = "raw_counts"),
              "Parameter 'assay_name' not found: Please provide a valid value.",
               fixed = TRUE)
  expect_warning(cms(sce = sce, k=20, group = "batch",
                   embedding = "tsne", assay_name = "counts"),
                 "Embedding not found: PCA subspace is used to calculate distances.",
               fixed = TRUE)
})


### Test for cms_summary

test_that("test output of cms_summary",{
  cms_res <- cms(sce, k = 20, group = "batch")
  cms_raw <- cms(sce, k = 20, group = "batch", smooth = FALSE)
  summary_res <- cms_summary(cms_res)
  summary_raw <- cms_summary(cms_raw)
  summary_batch <- cms_summary(cms_res, sum_var = "batch", sce = sce)
  summary_raw_batch <- cms_summary(cms_res, sum_var = "batch", sce = sce)

  expect_is(summary_res, "data.frame")
  expect_is(summary_raw, "data.frame")
  expect_equal(nrow(summary_batch), length(table(colData(sce)$batch)))
  expect_equal(nrow(summary_batch), nrow(summary_raw_batch))
  expect_error(cms_summary(cms_res, sum_var = "batch"),
               "Missing variable: Please provide a 'sce' object.")
  expect_error(cms_summary(cms_raw, sum_var = "quatsch", sce = sce),
               "Missing variable: Could not find 'sum_var'.
           Please specify one of names(colData(sce)).", fixed = TRUE)
})
