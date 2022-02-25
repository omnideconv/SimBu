library(SummarizedExperiment)
library(Matrix)

test_that('can generate different simulation scenarios', {
  counts <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  tpm <- Matrix::t(1e6*Matrix::t(tpm)/Matrix::colSums(tpm))

  colnames(counts) <- paste0("cell_",rep(1:300))
  colnames(tpm) <- paste0("cell_",rep(1:300))
  rownames(counts) <- paste0("gene_",rep(1:1000))
  rownames(tpm) <- paste0("gene_",rep(1:1000))

  annotation <- data.frame("ID"=paste0("cell_",rep(1:300)),
                           "cell_type"=c(rep("T cells CD4",50),
                                         rep("T cells CD8",50),
                                         rep("Macrophages",100),
                                         rep("NK cells",10),
                                         rep("B cells",70),
                                         rep("Monocytes",20)))

  dataset <- SimBu::dataset(annotation = annotation,
                            count_matrix = counts,
                            tpm_matrix = tpm,
                            name = "test_dataset")

  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'even', scaling_factor = 'NONE', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'NONE', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'mirror_db', scaling_factor = 'NONE', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'unique', unique_cell_type = 'T cells CD4', scaling_factor = 'NONE', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'controlled', controlled_cell_type = 'T cells CD4', controlled_amount = .5, scaling_factor = 'NONE', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
})


test_that('test different scaling factor calculations + mRNA bias removal from counts', {
  counts <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  tpm <- Matrix::t(1e6*Matrix::t(tpm)/Matrix::colSums(tpm))

  colnames(counts) <- paste0("cell_",rep(1:300))
  colnames(tpm) <- paste0("cell_",rep(1:300))
  rownames(counts) <- paste0("gene_",rep(1:1000))
  rownames(tpm) <- paste0("gene_",rep(1:1000))

  annotation <- data.frame("ID"=paste0("cell_",rep(1:300)),
                           "cell_type"=c(rep("T cells CD4",50),
                                         rep("T cells CD8",50),
                                         rep("Macrophages",100),
                                         rep("NK cells",10),
                                         rep("B cells",70),
                                         rep("Monocytes",20)),
                           "spikes" = runif(300))

  dataset <- SimBu::dataset(annotation = annotation,
                            count_matrix = counts,
                            tpm_matrix = tpm,
                            name = "test_dataset",
                            spike_in_col = 'spikes')

  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'census', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'spike_in', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'read_number', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'expressed_genes', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')

  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'NONE',remove_bias_in_counts = F, nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')
  expect_s4_class(SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'NONE',remove_bias_in_counts_method = 'gene-number', nsamples = 10, ncells = 100, ncores = 1)[['bulk']], 'SummarizedExperiment')

})


test_that('test different simulation parameters', {
  counts <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  tpm <- Matrix::t(1e6*Matrix::t(tpm)/Matrix::colSums(tpm))

  colnames(counts) <- paste0("cell_",rep(1:300))
  colnames(tpm) <- paste0("cell_",rep(1:300))
  rownames(counts) <- paste0("gene_",rep(1:1000))
  rownames(tpm) <- paste0("gene_",rep(1:1000))

  annotation <- data.frame("ID"=paste0("cell_",rep(1:300)),
                           "cell_type"=c(rep("T cells CD4",50),
                                         rep("T cells CD8",50),
                                         rep("Macrophages",100),
                                         rep("NK cells",10),
                                         rep("B cells",70),
                                         rep("Monocytes",20)),
                           "spikes" = runif(300))

  dataset <- SimBu::dataset(annotation = annotation,
                            count_matrix = counts,
                            tpm_matrix = tpm,
                            name = "test_dataset",
                            spike_in_col = 'spikes')

  sim_depth <- SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'NONE',total_read_counts = 100000, norm_counts = F ,nsamples = 10, ncells = 100, ncores = 1)
  expect_true(any(as.integer(colSums(assays(sim_depth$bulk)[['bulk_counts']]))==1e5))

  #sim_depth <- SimBu::simulate_bulk(data = dataset, scenario = 'random', scaling_factor = 'NONE',total_read_counts = 100, norm_counts = F ,nsamples = 10, ncells = 1000, ncores = 1)
  #expect_true(any(as.integer(colSums(assays(sim_depth$bulk)[['bulk_counts']]))==100))


})
