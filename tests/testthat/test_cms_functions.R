test_that("test that output class of cms.pca is matrix",{
  data(sim_30)
  sce <- sim_30[[1]]
  sce <- sce[, c(1:100)]
  cms_smooth <- cms.pca(sce = sce, k=20, group = "batch")
  cms_raw <- cms.pca(sce = sce, k=20, group = "batch", smooth = FALSE)
  cms_only <- as.matrix(cms_smooth[,"cms"])
  colnames(cms_only) <- "cms"

  expect_equal(class(cms_smooth), "matrix")
  expect_equal(class(cms_raw), "matrix")
  expect_identical(cms_raw, cms_only)
  expect_error(cms.pca(sce = sce, k=20, group = "batch", cell_min = 5),
               "Error: 'cell_min' is < 10. Must be > 10 to estimate cms.")
})
