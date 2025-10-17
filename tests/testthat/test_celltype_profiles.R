library(testthat)
library(SummarizedExperiment)
library(Matrix)

set.seed(123)
counts <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
tpm <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))

colnames(counts) <- paste0("cell_", rep(1:300))
colnames(tpm) <- paste0("cell_", rep(1:300))
rownames(counts) <- paste0("gene_", rep(1:1000))
rownames(tpm) <- paste0("gene_", rep(1:1000))

annotation <- data.frame(
  "ID" = paste0("cell_", rep(1:300)),
  "cell_type" = c(
    rep("T cells CD4", 50),
    rep("T cells CD8", 50),
    rep("Macrophages", 100),
    rep("NK cells", 10),
    rep("B cells", 70),
    rep("Monocytes", 20)
  )
)

dataset <- SimBu::dataset(
  annotation = annotation,
  count_matrix = counts,
  tpm_matrix = tpm,
  name = "test_dataset"
)

# Test backward compatibility - simulate_bulk without cell-type profiles
test_that("simulate_bulk works without cell-type profiles (backward compatibility)", {
  
  # Test without cell-type profiles (default behavior)
  simulation <- SimBu::simulate_bulk(
    dataset,
    scenario = "even",
    scaling_factor = "NONE",
    nsamples = 3,
    ncells = 100,
    generate_celltype_profiles = FALSE
  )
  
  # Should have traditional output structure
  expect_true("bulk" %in% names(simulation))
  expect_true("cell_fractions" %in% names(simulation))
  expect_true("scaling_vector" %in% names(simulation))
  expect_false("celltype_profiles" %in% names(simulation))
  
  # Check bulk data structure
  expect_s4_class(simulation$bulk, "SummarizedExperiment")
  expect_equal(ncol(simulation$bulk), 3)
  expect_equal(nrow(simulation$bulk), 1000)
})

# Test new feature - simulate_bulk with cell-type profiles
test_that("simulate_bulk works with cell-type profiles enabled", {

  # Test with cell-type profiles enabled
  simulation <- SimBu::simulate_bulk(
    dataset,
    scenario = "even",
    scaling_factor = "NONE",
    nsamples = 3,
    ncells = 100,
    generate_celltype_profiles = TRUE
  )
  
  # Should have extended output structure
  expect_true("bulk" %in% names(simulation))
  expect_true("cell_fractions" %in% names(simulation))
  expect_true("scaling_vector" %in% names(simulation))
  expect_true("celltype_profiles" %in% names(simulation))
  
  # Check cell-type profiles structure
  expect_true("counts" %in% names(simulation$celltype_profiles))
  expect_true("tpm" %in% names(simulation$celltype_profiles))
  
  # Check that we have profiles for each cell type
  expected_celltypes <- c("T cells CD4", "T cells CD8", "Macrophages", "NK cells", "B cells", "Monocytes")
  expect_true(all(expected_celltypes %in% names(simulation$celltype_profiles$counts)))
  expect_true(all(expected_celltypes %in% names(simulation$celltype_profiles$tpm)))
  
  # Check dimensions of cell-type specific matrices
  for (celltype in expected_celltypes) {
    expect_equal(nrow(simulation$celltype_profiles$counts[[celltype]]), 1000)
    expect_equal(ncol(simulation$celltype_profiles$counts[[celltype]]), 3)
    expect_equal(nrow(simulation$celltype_profiles$tpm[[celltype]]), 1000)
    expect_equal(ncol(simulation$celltype_profiles$tpm[[celltype]]), 3)
  }
})

# Test mathematical consistency
test_that("Cell-type profiles sum to bulk expression", {
   simulation <- SimBu::simulate_bulk(
    dataset,
    scenario = "even",
    scaling_factor = "NONE",
    nsamples = 2,
    ncells = 100,
    generate_celltype_profiles = TRUE
  )
  
  # Test that cell-type profiles sum to bulk expression (counts)
  bulk_counts <- as.matrix(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]])
  celltype_sum_counts <- Reduce(`+`, lapply(simulation$celltype_profiles$counts, as.matrix))
  
  # Allow for small numerical differences due to floating point arithmetic
  expect_true(all(abs(bulk_counts - celltype_sum_counts) < 1e-10))
  
  # Test that cell-type profiles sum to bulk expression (TPM)
  bulk_tpm <- as.matrix(SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]])
  celltype_sum_tpm <- Reduce(`+`, lapply(simulation$celltype_profiles$tpm, as.matrix))
  
  # Allow for small numerical differences due to floating point arithmetic and normalization
  expect_true(all(abs(bulk_tpm - celltype_sum_tpm) < 1e-6))
})

# Test with custom scenario
test_that("Cell-type profiles work with custom scenarios", {
  
  # Custom scenario with specific fractions
  fractions <- data.frame(
    "T cells CD4" = c(0.3, 0.2),
    "Macrophages" = c(0.5, 0.6),
    "B cells" = c(0.2, 0.2),
    check.names = FALSE
  )
  
  simulation <- SimBu::simulate_bulk(
    dataset,
    scenario = "custom",
    custom_scenario_data = fractions,
    scaling_factor = "NONE",
    nsamples = 2,
    ncells = 100,
    generate_celltype_profiles = TRUE
  )
  
  # Check that only the specified cell types have non-zero profiles
  specified_celltypes <- c("T cells CD4", "Macrophages", "B cells")
  for (celltype in specified_celltypes) {
    expect_true(celltype %in% names(simulation$celltype_profiles$counts))
    # Check that profiles are not all zeros
    expect_true(any(as.matrix(simulation$celltype_profiles$counts[[celltype]]) > 0))
  }
})

# Test edge cases
test_that("Cell-type profiles handle edge cases correctly", {
  
  # Test pure scenario (only one cell type)
  simulation_pure <- SimBu::simulate_bulk(
    dataset,
    scenario = "pure",
    pure_cell_type = "B cells",
    scaling_factor = "NONE",
    nsamples = 2,
    ncells = 100,
    generate_celltype_profiles = TRUE
  )
  
  # Only B cells should have non-zero profiles
  expect_true("B cells" %in% names(simulation_pure$celltype_profiles$counts))
  expect_true(any(as.matrix(simulation_pure$celltype_profiles$counts[["B cells"]]) > 0))
  
  # Other cell types should have zero profiles if they exist
  other_celltypes <- setdiff(names(simulation_pure$celltype_profiles$counts), "B cells")
  for (celltype in other_celltypes) {
    expect_true(all(as.matrix(simulation_pure$celltype_profiles$counts[[celltype]]) == 0))
  }
})
