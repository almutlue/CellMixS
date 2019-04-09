library(SingleCellExperiment)
library(CellMixS)

#get simulated scRNA seq data with 3 unbalanced batches
sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
sce <- sim_list[["batch20"]][, c(1:50,300:350)]
sce_batch1 <- sce[,colData(sce)$batch == "1"]
sce_batch2 <- sce[,colData(sce)$batch == "2"]
sce_pre_list <- list("1" = sce_batch1, "2" = sce_batch2)

## Tests for ldf fuctions:
### Include internal functions as ldfKnn and defineSubspace.

test_that("test that output of ldfDiff is correct",{
    sce_ldf_mnn <- ldfDiff(sce_pre_list, sce, k = 10, group = "batch",
                           dim_combined = "MNN", n_dim = 5, res_name = "MNN")
    sce_noDim <- sce
    reducedDims(sce_noDim) <-list()
    sce_ldf_sameDim <- ldfDiff(sce_pre_list, sce, k = 10, group = "batch",
                               n_dim = 5)
    sce_ldf_tsne <- ldfDiff(sce_pre_list, sce, k = 10, group = "batch",
                            dim_combined = "TSNE", n_dim = 2)

    expect_is(sce_ldf_mnn, "SingleCellExperiment")
    expect_equal(sum(sce_ldf_sameDim$diff_ldf), 0)
    expect_is(sce_ldf_mnn$diff_ldf.MNN, "numeric")
    expect_is(sce_ldf_tsne, "SingleCellExperiment")
    expect_error(ldfDiff(sce_pre_list, sce, k = 500, group = "batch",
                         n_dim = 5),
                 "Parameter 'k' is greater than dataset size:
            Please provide a valid value.")

    expect_error(ldfDiff(sce_pre_list, sce, k = 10, group = "batch",
                         dim_red = "tsne", assay_pre = "raw.counts"),
                 "Ambigious parameter: Specify subspace parameter.
         * For precalculated embeddings keep 'assay_name' as default.
         * For PCA based on 'assay_name' keep 'dim_red' as default.",
                 fixed = TRUE)

    expect_error(ldfDiff(sce_pre_list, sce, k = 10, group = "batch",
                         dim_combined = "tsne", assay_combined = "raw.counts"),
                 "Ambigious parameter: Specify subspace parameter.
         * For precalculated embeddings keep 'assay_name' as default.
         * For PCA based on 'assay_name' keep 'dim_red' as default.",
                 fixed = TRUE)

    expect_error(ldfDiff(sce_pre_list, sce_noDim, k = 10, group = "batch"),
                 "Parameter 'assay_name' not found: Provide a valid value.",
                 fixed = TRUE)
    expect_error(ldfDiff(sce_pre_list, sce, k = 10, group = "batch",
                         dim_red = "TSNE", n_dim = 5),
                 "Parameter 'n_dim' is greater than reduced dimensional space:
         Please provide a valid value.",
                 fixed = TRUE)
    expect_warning(ldfDiff(sce_pre_list, sce_noDim, k = 10, group = "batch",
                           assay_pre = "counts", assay_combined = "counts",
                           n_dim = 5),
                   "'dim_red' not found:
            PCA subspace is used to calculate distances.",
                   fixed = TRUE)

})



