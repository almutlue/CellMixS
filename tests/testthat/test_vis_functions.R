library(SingleCellExperiment)
library(CellMixS)
load(system.file("extdata/sim30.rda", package = "CellMixS"))
load(system.file("extdata/cms_sim30.rda", package = "CellMixS"))
sce <- sim_30[[1]][, c(1:50,500:550)]

## Tests for visualization fuctions:

### visHist
test_that("test that visHist works",{
  hist <- visHist(cms_sim30)
  hist_cms <- visHist(data.frame("cms" = cms_sim30[, "cms"]))
  expect_is(hist, "ggplot")
  expect_is(hist_cms, "ggplot")

})

### visOverview

test_that("test that visOverview works",{
  overview <- visOverview(cms_sim30, sce, "batch")
  #change embeddings
  overview_mnn <- visOverview(cms_sim30, sce, "batch", dim_red = "MNN", log10_val = TRUE)
  #data frame as input
  overview_cms <- visOverview(data.frame("cms" = cms_sim30[, "cms"]), sce, "batch")
  #no smoothed results
  overview_cms2 <- visOverview(cms_sim30, sce, "batch", smooth = FALSE)
  #calculate new dim_red
  sce_noRedDim <- sce
  reducedDims(sce_noRedDim) <- NULL
  overview_tsne <- visOverview(cms_sim30, sce_noRedDim, "batch")

  expect_is(overview, "gg")
  expect_is(overview_mnn, "gg")
  expect_is(overview_cms, "gg")
  expect_is(overview_cms2, "gg")
  expect_error(visOverview(cms_sim30, sce, "batch", dim_red = "quatsch"),
               "Ambigous parameter 'dim_red', provide one of:
           * A dim_red method that is listed in reducedDimNames(sce).
           * Default('TSNE') will call runTSNE to calculate a subspace.",
               fixed = TRUE)
})


### visCMS and visGroup

test_that("test that visCMS and visGroup work",{
  #visCMS
  vis_cms <- visCms(cms_sim30, sce)
  vis_cms2 <- visCms(cms_sim30, sce, cms_var = "cms_smooth", dim_red = "MNN")
  #visGroup
  sce_num <- sce
  colData(sce_num)$batch <- as.numeric(colData(sce)$batch)
  vis_group <- visGroup(sce, "batch")
  vis_group2 <- visGroup(sce_num, "batch", dim_red = "MNN")

  expect_is(vis_cms, "gg")
  expect_is(vis_cms2, "gg")
  expect_is(vis_group, "gg")
  expect_is(vis_group2, "gg")
})

