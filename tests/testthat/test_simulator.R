library(SummarizedExperiment)
library(Matrix)

counts <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol = 300), sparse = TRUE)
tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol = 300), sparse = TRUE)
tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))

colnames(counts) <- paste0("cell_", rep(1:300))
colnames(tpm) <- paste0("cell_", rep(1:300))
rownames(counts) <- paste0("gene-", rep(1:1000))
rownames(tpm) <- paste0("gene-", rep(1:1000))

annotation <- data.frame(
  "ID" = paste0("cell_", rep(1:300)),
  "cell_type" = c(
    rep("T cells CD4", 50),
    rep("T cells CD8", 50),
    rep("Macrophages", 100),
    rep("NK cells", 10),
    rep("B cells", 70),
    rep("Monocytes", 20)
  ),
  "spikes" = runif(300)
)

dataset <- SimBu::dataset(
  annotation = annotation,
  count_matrix = counts,
  tpm_matrix = tpm,
  name = "test_dataset",
  spike_in_col = "spikes"
)

test_that("can generate different simulation scenarios & check multi threads", {
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "even", scaling_factor = "NONE", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "NONE", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "mirror_db", scaling_factor = "NONE", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "pure", pure_cell_type = "T cells CD4", scaling_factor = "NONE", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "weighted", weighted_cell_type = "T cells CD4", weighted_amount = .5, scaling_factor = "NONE", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")

  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "even", scaling_factor = "NONE", nsamples = 10, ncells = 100, BPPARAM = BiocParallel::MulticoreParam(workers = 2), run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
})

<<<<<<< HEAD
=======
test_that("test RNG", {
  # use even scenario with no variance between samples
  # with fixed seed, both simulations have exactly the same count values
  seed <- 123
  sim1 <- SimBu::simulate_bulk(data = dataset, scenario = "even", scaling_factor = "NONE", balance_even_mirror_scenario = 0, nsamples = 10, ncells = 100, run_parallel = FALSE, seed = seed)
  sim2 <- SimBu::simulate_bulk(data = dataset, scenario = "even", scaling_factor = "NONE", balance_even_mirror_scenario = 0, nsamples = 10, ncells = 100, run_parallel = FALSE, seed = seed)
  x1 <- Matrix::rowSums(assays(sim1$bulk)[['bulk_counts']])
  x2 <- Matrix::rowSums(assays(sim2$bulk)[['bulk_counts']])
  testthat::expect_equal(x1,x2)
  
  # test that samples inside one simulation still are different
  sample1 <- sim1$bulk[,1]
  sample2 <- sim1$bulk[,2]
  testthat::expect_false(all(assays(sample1)[['bulk_counts']] == assays(sample2)[['bulk_counts']]))

})
>>>>>>> f6de93e (use correct function in test)

test_that("test different scaling factor calculations + mRNA bias removal from counts", {
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "census", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "spike_in", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "read_number", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "expressed_genes", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "expressed_genes", nsamples = 10, ncells = 100, run_parallel = FALSE, scaling_factor_single_cell = FALSE)[["bulk"]], "SummarizedExperiment")

  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "epic", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_warning(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "abis", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "For some cell type")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "quantiseq", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_warning(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "custom", custom_scaling_vector = c("B cells" = 10), nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "For some cell type")

  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "NONE", remove_bias_in_counts = TRUE, remove_bias_in_counts_method = "read-number", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "NONE", remove_bias_in_counts = FALSE, nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
  testthat::expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = "random", scaling_factor = "NONE", remove_bias_in_counts_method = "gene-number", nsamples = 10, ncells = 100, run_parallel = FALSE)[["bulk"]], "SummarizedExperiment")
})




test_that("test census", {
  tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol = 300), sparse = TRUE)
  tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))

  cen <- expect_equal(length(SimBu::census(tpm, exp_capture_rate = 0.25, expr_threshold = 0.1, run_parallel = FALSE)), 300)
})

# test_that('test different simulation parameters', {
#   counts <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
#   tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
#   tpm <- Matrix::t(1e6*Matrix::t(tpm)/Matrix::colSums(tpm))
#
#   colnames(counts) <- paste0("cell_",rep(1:300))
#   colnames(tpm) <- paste0("cell_",rep(1:300))
#   rownames(counts) <- paste0("gene_",rep(1:1000))
#   rownames(tpm) <- paste0("gene_",rep(1:1000))
#
#   annotation <- data.frame("ID"=paste0("cell_",rep(1:300)),
#                            "cell_type"=c(rep("T cells CD4",50),
#                                          rep("T cells CD8",50),
#                                          rep("Macrophages",100),
#                                          rep("NK cells",10),
#                                          rep("B cells",70),
#                                          rep("Monocytes",20)),
#                            "spikes" = runif(300))
#
#   dataset <- SimBu::dataset(annotation = annotation,
#                             count_matrix = counts,
#                             tpm_matrix = tpm,
#                             name = "test_dataset",
#                             spike_in_col = 'spikes')
#
#   sim_depth <- SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'NONE',total_read_counts = 100000, norm_counts = FALSE ,nsamples = 10, ncells = 100, ncores = 1)
#   expect_true(as.integer(colSums(assays(sim_depth$bulk)[['bulk_counts']]))[1] == 1e5)
#
#   #sim_depth <- SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'NONE',total_read_counts = 100, norm_counts = FALSE ,nsamples = 10, ncells = 1000, ncores = 1)
#   #expect_true(any(as.integer(colSums(assays(sim_depth$bulk)[['bulk_counts']]))==100))
#
#
# })
