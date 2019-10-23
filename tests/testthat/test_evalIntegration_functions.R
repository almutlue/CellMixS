library(SingleCellExperiment)
library(CellMixS)
#library(Seurat)
library(scater)
#get simulated scRNA seq data with 3 unbalanced batches
sim_list <- readRDS(system.file("extdata/sim50.rds", package = "CellMixS"))
sce <- sim_list[[1]][, c(1:30,300:320)]
sce_batch1 <- sce[,colData(sce)$batch == "1"]
sce_batch2 <- sce[,colData(sce)$batch == "2"]
pre <- list("1" = sce_batch1, "2" = sce_batch2)

### generate Mixing metric using seurat (to compare implementations)
#Run mixing metric
# sce <- logNormCounts(sce)
# rownames(sce) <- paste0("gene", seq_len(nrow(sce)))
# seurat <- as.Seurat(sce)
# mix_dist <- MixingMetric(seurat, grouping.var = "batch",
#                          reduction = "PCA",
#                          dims = seq_len(10),
#                          k = 5, max.k = 50)



## Tests for evalIntegration, isis and entropy fuction:

test_that("test that output of evalIntegration is correct",{
    #input
    sce_all <- evalIntegration(metrics = c("cms", "isi", "entropy", "ldfDiff"),
                               sce, sce_pre_list = pre, "batch", k = 10)
    sce_seur <- evalIntegration(metrics = c("localStructure", "mixingMetric"),
                                sce, "batch", k = 10, n_dim = 2, n_combined = 2,
                                assay_name = "counts")
    sce_all <- evalIntegration("isi", sce_all, "batch", k = 10, weight = FALSE,
                               res_name = "wisi")
    sce_isi <- isi(sce, "batch", k = 10)
    sce_entropy <- entropy(sce, "batch", k = 10)
    sce_mix <- mixMetric(sce, "batch", k = 50, k_pos = 5)
    sce_ls <- locStructure(sce, "batch", dim_combined = "TSNE",
                           k = 20, assay_name = "counts", n_combined = 2)

    #test output
    expect_is(sce_all, "SingleCellExperiment")
    expect_is(sce_seur, "SingleCellExperiment")
    expect_equal(ncol(colData(sce_all)) - ncol(colData(sce)), 6)
    expect_false(all(sce_all$isi == sce_all$wisi))
    expect_equal(sce_all$isi, sce_isi$isi)
    #expect_equal(mix_dist, sce_mix$mm)
    expect_equal(sce_all$entropy, sce_entropy$entropy)
    expect_false(length(levels(sce$batch)) == length(levels(sce_isi$batch)))
    expect_true(length(levels(sce_isi$batch)) == length(levels(sce_all$batch)))

    #test errors and warnings
    #entropy
    expect_error(entropy(sce$batch, "batch", k = 20),
                 "Error: 'sce' must be a 'SingleCellExperiment' object.",
                 fixed = TRUE)
    expect_error(entropy(sce, "quatsch", k = 20),
                 "Error: 'group' variable must be in 'colData(sce)'",
                 fixed = TRUE)
    expect_warning(entropy(sce, "batch", k = 100),
                   "'k' exceeds number of cells. Is set to max (all cells).",
                 fixed = TRUE)
    #isi
    expect_error(isi(sce$batch, "batch", k = 20),
                 "Error: 'sce' must be a 'SingleCellExperiment' object.",
                 fixed = TRUE)
    expect_error(isi(sce, "quatsch", k = 20),
                 "Error: 'group' variable must be in 'colData(sce)'",
                 fixed = TRUE)
    expect_warning(isi(sce, "batch", k = 100),
                   "'k' exceeds number of cells. Is set to max (all cells).",
                   fixed = TRUE)
    #evalIntegration
    expect_error(evalIntegration(metrics = c("isi", "entropy", "new_quatsch"),
                                 sce, "batch", k = 20),
                 "Error: 'metrics' is unknown. Please define one or more of 'cms', 'isi',
             'ldfDiff', 'mixingMetric', 'localStructure', 'entropy'",
                 fixed = TRUE)
    expect_error(evalIntegration("new_quatsch", sce, "batch", k = 20),
                 "Error: 'metrics' is unknown. Please define one or more of 'cms', 'isi',
             'ldfDiff', 'mixingMetric', 'localStructure', 'entropy'",
                 fixed = TRUE)
    expect_error(evalIntegration("mixingMetric", sce$batch, "batch", k = 20),
                 "Error: 'sce' must be a 'SingleCellExperiment' object.",
                 fixed = TRUE)
    expect_error(evalIntegration("mixingMetric", sce, "batch", k = 20,
                                 res_name = c("mm", 'mm_smooth')),
                 "Error: Define 'res_name' for all metrics to calculate'",
                 fixed = TRUE)
    expect_error(evalIntegration("entropy", sce, "quatsch", k = 20),
                 "Error: 'group' variable must be in 'colData(sce)'",
                 fixed = TRUE)
    expect_error(evalIntegration("cms", sce, "batch"),
                   "Please specify 'k', the number of nearest neigbours to check
                 for equal mixing, e.g. median of cells/celltype.",
                   fixed = TRUE)
    expect_error(evalIntegration("isi", sce, "batch"),
                 "Please specify 'k', the number of nearest neigbours to check
                 for equal mixing, e.g. median of cells/celltype.",
                 fixed = TRUE)
    expect_error(evalIntegration("ldfDiff", sce, "batch"),
                 "Please specify 'k', the number of nearest neigbours to check
                 for structual changes",
                 fixed = TRUE)
    expect_error(evalIntegration("entropy", sce, "batch"),
                 "Please specify 'k', the number of nearest neigbours to check
                 for equal mixing, e.g. median of cells/celltype.",
                 fixed = TRUE)
    expect_warning(evalIntegration("mixingMetric", sce, "batch", k = 100),
                   "'k' exceeds number of cells. Is set to max (all cells).",
                   fixed = TRUE)
    expect_error(evalIntegration("localStructure", sce, "batch",
                                 assay_name = "test", dim_combined = "MNN",
                                 k = 20),
                 "Error: 'assay_name' not found.", fixed = TRUE)
    expect_error(evalIntegration("localStructure", sce, "batch",
                                 dim_combined = "MNN"),
                 "'k' exceeds number of cells/batch.",
                 fixed = TRUE)
    expect_error(evalIntegration("localStructure", sce, "batch",
                                 assay_name = "counts", k = 20,
                                 dim_combined = "Quatsch"),
                 "Error: 'dim_combined' not found.",
                 fixed = TRUE)
    expect_error(locStructure( sce, "batch", dim_combined = "TSNE",
                                 assay_name = "counts", k = 20),
                 "Error: 'n_combined' exceeds the dimensions of 'dim_combined'.",
                 fixed = TRUE)
})


