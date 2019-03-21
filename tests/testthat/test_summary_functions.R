library(SingleCellExperiment)
library(CellMixS)
load(system.file("extdata/sim30.rda", package = "CellMixS"))
load(system.file("extdata/cms_sim30.rda", package = "CellMixS"))
sce <- sim_30[[1]][, c(1:50,500:550)]



### Test for cms_summary

test_that("test output of cms_summary",{
  cms_raw <- cms(sce, k = 20, group = "batch", smooth = FALSE)
  summary_res <- cmsSummary(cms_sim30)
  summary_raw <- cmsSummary(cms_raw)
  summary_batch <- cmsSummary(cms_sim30, sum_var = "batch", sce = sce)
  summary_raw_batch <- cmsSummary(cms_raw, sum_var = "batch", sce = sce)

  expect_is(summary_res, "data.frame")
  expect_is(summary_raw, "data.frame")
  expect_equal(nrow(summary_batch), length(table(droplevels(colData(sce)$batch))))
  expect_equal(nrow(summary_batch), nrow(summary_raw_batch))
  expect_error(cmsSummary(cms_res, sum_var = "batch"),
               "Missing variable: Please provide a 'sce' object.")
  expect_error(cmsSummary(cms_raw, sum_var = "quatsch", sce = sce),
               "Missing variable: Could not find 'sum_var'.
           Please specify one of names(colData(sce)).", fixed = TRUE)
})



### compareIntegration and compareCluster
test_that("test that compareIntegration and compareCluster work",{
  #compareIntegration
  cms_mnn <- cms(sce, k = 30, group = "batch", dim_red = "MNN")
  cms_list <- list("raw"= cms_sim30[,"cms"], "mnn" = cms_mnn[,"cms"])
  compare_methods <- compareIntegration(cms_list)
  compare_methods_df <- compareIntegration(cms_mnn)
  #compare groups
  colData(sce)$group <- sample(c("A", "B", "C"), ncol(sce), replace = TRUE)
  cms_group <- data.frame("cms"= cms_sim30[,"cms"],
                          "group" = sample(c("A", "B", "C"), nrow(cms_sim30), replace = TRUE))
  compare_group <- compareCluster(cms_sim30, "group", cms_var = "cms", sce = sce, violin = TRUE)
  compare_group2 <- compareCluster(data.frame("cms" = cms_sim30[,"cms"]), "batch", sce = sce)
  compare_group3 <- compareCluster(cms_group, "group")

  expect_is(compare_methods, "gg")
  expect_is(compare_methods_df, "gg")
  expect_is(compare_group, "gg")
  expect_is(compare_group2, "gg")
  expect_is(compare_group3, "gg")
})
