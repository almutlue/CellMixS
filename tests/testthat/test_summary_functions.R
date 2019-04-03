library(SingleCellExperiment)
library(CellMixS)
load(system.file("extdata/sim30.rda", package = "CellMixS"))
sce <- sim_30[[1]][, c(1:50,500:550)]
sce_cms <- cms(sce,"batch", k = 30, res_name = "unaligned")
sce_mnn <- cms(sce_cms,"batch", k = 30, dim_red = "MNN", res_name = "MNN")
cms_list <- list("unaligned"= sce_cms$cms.unaligned, "mnn" = sce_mnn$cms.MNN)
cms_df <- data.frame("unaligned"= sce_cms$cms.unaligned, "mnn" = sce_mnn$cms.MNN)



### Test for summary plots

### visIntegration and visCluster
test_that("test that visIntegration and visCluster work",{
    #visIntegration
    compare_Int_list <- visIntegration(cms_list)
    compare_Int_sce<- visIntegration(sce_mnn, metric_prefix = "cms.")
    compare_Int_df <- visIntegration(cms_df, violin = TRUE)

    #compare groups
    colData(sce_cms)$group <- sample(c("A", "B", "C"), ncol(sce), replace = TRUE)
    compare_group <- visCluster(sce_cms, "group", metric_var = "cms.unaligned")
    compare_group2 <- visCluster(sce_mnn, "batch", metric_var = "cms.MNN")


    expect_is(compare_Int_list, "gg")
    expect_is(compare_Int_sce, "gg")
    expect_is(compare_Int_df, "gg")
    expect_is(compare_group, "gg")
    expect_is(compare_group2, "gg")
    expect_error(visCluster(sce_mnn, "batch"))
})
