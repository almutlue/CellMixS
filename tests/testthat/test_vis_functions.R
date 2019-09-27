library(SingleCellExperiment)
library(CellMixS)
sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
sce <- sim_list[["batch20"]][, c(1:30,300:320)]
sce_cms <- cms(sce, "batch", k = 20, n_dim = 2)

## Tests for visualization fuctions:

### visHist
test_that("test that visHist works",{
    hist_cms <- visHist(data.frame("cms" = sce_cms$cms))
    hist_sce <- visHist(sce_cms)
    expect_is(hist_cms, "ggplot")
    expect_is(hist_sce, "ggplot")
    expect_error(visHist(sce),
                 "Error: 'res_object' does not contain any metric results.
             Please continue by one of:
             * Run `cms` on your SingleCellExperiment object before plotting.
             * Specify colData(res_object) column to plot by `metric`.
             * Specify a matrix with results to plot as `res_object`.",
                 fixed = TRUE)

})

### visOverview

test_that("test that visOverview works",{
    overview <- visOverview(sce_cms, "batch")
    #change embeddings
    overview_mnn <- visOverview(sce_cms, "batch", dim_red = "MNN",
                                log10_val = TRUE)
    #add other vars (continous and discrete)
    overview_Vars <- visOverview(sce_cms, "batch",
                                 other_var = c("batch", "cms"))

    #calculate new dim_red
    sce_noRedDim <- sce_cms
    reducedDims(sce_noRedDim) <- NULL
    overview_tsne <- visOverview(sce_noRedDim, "batch")

    expect_is(overview, "gg")
    expect_is(overview_mnn, "gg")
    expect_is(overview_Vars, "gg")
    expect_is(overview_tsne, "gg")
    expect_error(visOverview(sce_cms, "batch", dim_red = "quatsch"),
                 "Ambigous parameter 'dim_red', provide one of:
                 * A dim_red method that is listed in reducedDimNames(sce_cms).
                 * Default('TSNE') will call runTSNE to calculate a subspace.",
                 fixed = TRUE)
    expect_error(visOverview(sce, "batch"),
                 "Error: 'sce_cms' does not contain any metric results.
             Please continue by one of:
             * Run `cms` on your SingleCellExperiment object before plotting.
             * Specify a colData(res_object) column in `metric`.",
               fixed = TRUE)
    expect_error(visOverview(assay(sce), "batch"),
                 "Error:'sce_cms' must be a 'SingleCellExperiment' object.")
    expect_error(visOverview(sce, "batch2"),
                 "Error: 'group' variable must be in 'colData(sce_cms)'",
                 fixed = TRUE)
})


### visMetric and visGroup

test_that("test that visMetric and visGroup work",{
  #visMetric
  vis_cms <- visMetric(sce_cms)
  vis_cms2 <- visMetric(sce_cms, metric_var = "cms_smooth", dim_red = "MNN")
  #visGroup
  sce_num <- sce
  colData(sce_num)$batch <- as.numeric(colData(sce)$batch)
  vis_group <- visGroup(sce, "batch")
  vis_group2 <- visGroup(sce_num, "batch", dim_red = "MNN")

  expect_is(vis_cms, "gg")
  expect_is(vis_cms2, "gg")
  expect_is(vis_group, "gg")
  expect_is(vis_group2, "gg")
  expect_error(visMetric(assay(sce_cms), "batch"),
               "Error:'sce_cms' must be a 'SingleCellExperiment' object.")
  expect_error(visMetric(sce_cms, "cms2"),
               "Error: 'metric_var' variable must be in 'colData(sce_cms)'",
               fixed = TRUE)
  expect_error(visGroup(assay(sce), "batch"),
               "Error:'sce' must be a 'SingleCellExperiment' object.")
  expect_error(visGroup(sce, "batch2"),
               "Error: 'group' variable must be in 'colData(sce)'",
               fixed = TRUE)
})

