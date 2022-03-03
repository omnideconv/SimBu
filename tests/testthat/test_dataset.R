library(Matrix)
library(SummarizedExperiment)

test_that('can create dataset with counts or counts+tpm', {
  counts <-  Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  tpm <- Matrix::t(1e6*Matrix::t(tpm)/Matrix::colSums(tpm))

  colnames(counts) <- paste0("cell_",rep(1:300))
  colnames(tpm) <- paste0("cell_",rep(1:300))
  rownames(counts) <- paste0("gene_",rep(1:1000))
  rownames(tpm) <- paste0("gene_",rep(1:1000))

  annotation <- data.frame("ID" = paste0("cell_",rep(1:300)),
                           "cell_type" = c(rep("T cells CD4",300)))

  testthat::expect_s4_class(SimBu::dataset(annotation = annotation, count_matrix = counts, tpm_matrix = tpm, name = "test_dataset"), 'SummarizedExperiment')
  testthat::expect_s4_class(SimBu::dataset(annotation = annotation, count_matrix = counts, name = "test_dataset"), 'SummarizedExperiment')
  testthat::expect_error(SimBu::dataset(annotation = annotation, tpm_matrix = tpm, name = "test_dataset"))
})


test_that('carry over additional columns from annotation + have nReads_SimBu and nGenes_SimBu in anno', {
  counts <-  Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
  colnames(counts) <- paste0("cell_",rep(1:300))
  rownames(counts) <- paste0("gene_",rep(1:1000))

  annotation <- data.frame("ID"=paste0("cell_",rep(1:300)),
                           "cell_type" = c(rep("T cells CD4",300)),
                           'spikes' = runif(300),
                           'add_1' = runif(300),
                           'add_2' = runif(300))

  ds <- SimBu::dataset(annotation = annotation,
                       count_matrix = counts,
                       name = "test_dataset",
                       spike_in_col = 'spikes',
                       additional_cols = c('add_1','add_2'))

  anno_ds <- data.frame(SummarizedExperiment::colData(ds))

  testthat::expect_true(all(c('add_1','add_2') %in% colnames(anno_ds)))
  testthat::expect_true('spike_in' %in% colnames(anno_ds))
  testthat::expect_true('nReads_SimBu' %in% colnames(anno_ds))
  testthat::expect_true('nGenes_SimBu' %in% colnames(anno_ds))

})
