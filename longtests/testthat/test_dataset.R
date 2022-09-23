library(Matrix)
library(SummarizedExperiment)


setup_list <- SimBu::setup_sfaira(tempdir())

test_that('check sfaira for single dataset', {
  testthat::expect_s4_class(SimBu::dataset_sfaira(sfaira_id = 'homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x',sfaira_setup = setup_list, name = "test_dataset"), 'SummarizedExperiment')
})


test_that('check sfaira for multiple datasets', {
  testthat::expect_s4_class(SimBu::dataset_sfaira_multiple(sfaira_setup = setup_list,
                                                           organisms = "Homo sapiens",
                                                           tissues = "lung parenchyma",
                                                           assay = "10x 3' v2",
                                                           name = "human_lung"), 'SummarizedExperiment')
})
