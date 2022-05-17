library(Matrix)
library(SummarizedExperiment)

counts <-  Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol=300), sparse = TRUE)
tpm <- Matrix::t(1e6*Matrix::t(tpm)/Matrix::colSums(tpm))

colnames(counts) <- paste0("cell_",rep(1:300))
colnames(tpm) <- paste0("cell_",rep(1:300))
rownames(counts) <- paste0("gene-",rep(1:1000))
rownames(tpm) <- paste0("gene-",rep(1:1000))

test_that('can create dataset with counts or counts+tpm', {

  annotation <- data.frame("ID" = paste0("cell_",rep(1:300)),
                           "cell_type" = c(rep("T cells CD4",300)))

  testthat::expect_s4_class(SimBu::dataset(annotation = annotation, count_matrix = counts, tpm_matrix = tpm, name = "test_dataset"), 'SummarizedExperiment')
  testthat::expect_s4_class(SimBu::dataset(annotation = annotation, count_matrix = counts, name = "test_dataset"), 'SummarizedExperiment')
  testthat::expect_error(SimBu::dataset(annotation = annotation, tpm_matrix = tpm, name = "test_dataset"))
})


test_that('carry over additional columns from annotation + have nReads_SimBu and nGenes_SimBu in anno', {

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

test_that('can create dataset from seurat object', {
  
  annotation <- data.frame("ID"=paste0("cell_",rep(1:300)),
                           "cell_type"=c(rep("T cells CD4",50),
                                         rep("T cells CD8",50),
                                         rep("Macrophages",100),
                                         rep("NK cells",10),
                                         rep("B cells",70),
                                         rep("Monocytes",20)),
                           row.names = paste0("cell_",rep(1:300)))
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, assay = 'counts', meta.data = annotation)
  tpm_assay <- Seurat::CreateAssayObject(counts = tpm)
  seurat_obj[['tpm']] <- tpm_assay
  
  testthat::expect_s4_class(SimBu::dataset_seurat(seurat_obj = seurat_obj, count_assay = "counts",cell_id_col = 'ID', cell_type_col = 'cell_type', tpm_assay = 'tpm', name = "seurat_dataset"), 'SummarizedExperiment')
})


test_that('can load h5ad file with cells in obs and var', {
  h5 <- system.file('extdata', 'anndata.h5ad', package='SimBu')
  h5_rev <- system.file('extdata', 'anndata_rev.h5ad', package='SimBu')

  testthat::expect_s4_class(SimBu::dataset_h5ad(h5ad_file_counts = h5,name = "h5ad_dataset",cell_id_col = 'id',cell_type_col = 'group', cells_in_obs = TRUE), 'SummarizedExperiment')
  testthat::expect_s4_class(SimBu::dataset_h5ad(h5ad_file_counts = h5_rev,name = "h5ad_dataset",cell_id_col = 'id',cell_type_col = 'group', cells_in_obs = FALSE), 'SummarizedExperiment')
})


# test_that('check sfaira connection', {
#   setup_list <- SimBu::setup_sfaira(tempdir())
#   testthat::expect_s4_class(SimBu::dataset_sfaira(sfaira_id = 'homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x',sfaira_setup = setup_list, name = "test_dataset"), 'SummarizedExperiment')
# })
