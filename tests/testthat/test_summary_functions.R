library(SingleCellExperiment)
library(CellMixS)
sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
sce <- sim_list[["batch20"]][, c(1:30,300:320)]
sce_cms <- cms(sce,"batch", k = 20, res_name = "unaligned", n_dim = 2)
sce_mnn <- cms(sce_cms,"batch", k = 20, dim_red = "MNN", res_name = "MNN",
               n_dim = 2)
cms_list <- list("unaligned"= sce_cms$cms.unaligned, "mnn" = sce_mnn$cms.MNN)
cms_df <- data.frame("unaligned"= sce_cms$cms.unaligned,
                     "mnn" = sce_mnn$cms.MNN)



### Test for summary plots

### visIntegration and visCluster
test_that("test that visIntegration and visCluster work",{
    #visIntegration
    compare_Int_list <- visIntegration(cms_list)
    compare_Int_sce<- visIntegration(sce_mnn, metric = "cms.",
                                     metric_name = "cms")
    compare_Int_df <- visIntegration(cms_df, violin = TRUE)

    #compare groups
    colData(sce_cms)$group <- sample(c("A", "B", "C"), ncol(sce),
                                     replace = TRUE)
    compare_group <- visCluster(sce_cms, "group", metric_var = "cms.unaligned")
    compare_group2 <- visCluster(sce_mnn, "batch", metric_var = "cms.MNN")


    expect_is(compare_Int_list, "gg")
    expect_is(compare_Int_sce, "gg")
    expect_is(compare_Int_df, "gg")
    expect_is(compare_group, "gg")
    expect_is(compare_group2, "gg")
    expect_error(visIntegration(sce))
    expect_error(visCluster(sce_mnn, "batch"),
                 "Error: 'metric_var' variable must be in 'colData(sce_cms)'",
                 fixed = TRUE)
    expect_error(visCluster(sce_mnn, "batch2"),
                 "Error: 'cluster_var' variable must be in 'colData(sce_cms)'",
                 fixed = TRUE)
    expect_error(visCluster(cms_df, "batch"),
                 "Error: 'sce_cms' must be a 'SingleCellExperiment' object.",
                 fixed = TRUE)
})
