#' Generate SummarizedExperiment using multiple parameters
#'
#' @param annotation (mandatory) dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param count_matrix (mandatory) sparse count matrix; raw count data is expected with genes in rows, cells in columns
#' @param tpm_matrix sparse count matrix; TPM like count data is expected with genes in rows, cells in columns
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param spike_in_col which column in annotation contains information on spike_in counts, which can be used to re-scale counts; mandatory for spike_in scaling factor in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return Return a \link[SummarizedExperiment]{SummarizedExperiment} object
#'
#' @keywords internal
generate_summarized_experiment <- function(annotation, count_matrix, tpm_matrix, name, spike_in_col, additional_cols, filter_genes, variance_cutoff, type_abundance_cutoff, scale_tpm) {
  if (is.null(count_matrix)) {
    stop("A count_matrix is required to generate a dataset.")
  }

  if (is.null(count_matrix) && is.null(tpm_matrix)) {
    stop("count_matrix and tpm_matrix are both empty. At least one is required.")
  }

  # check if both matrices (if provided) have the same dimensions, row-names & col-names
  if (!is.null(count_matrix) && !is.null(tpm_matrix)) {
    if (any(dim(count_matrix) != dim(tpm_matrix))) {
      stop("count_matrix and tpm_matrix have different dimensions. Please check again.")
    }
    if (!all(colnames(count_matrix) %in% colnames(tpm_matrix))) {
      stop("Column names of count_matrix and tpm_matrix are not in common. Please check again.")
    }
    if (!all(rownames(count_matrix) %in% rownames(tpm_matrix))) {
      stop("Row names of count_matrix and tpm_matrix are not in common. Please check again.")
    }
  }

  # convert matrices to sparse matrix
  tryCatch(
    {
      count_matrix <- methods::as(count_matrix, "dgCMatrix")
    },
    error = function(e) {
      em <- paste0("Cannot convert count matrix in sparse matrix (dgCMatrix):", e$message)
      stop(em)
    }
  )
  if (!is.null(tpm_matrix)) {
    tryCatch(
      {
        tpm_matrix <- methods::as(tpm_matrix, "dgCMatrix")
      },
      error = function(e) {
        em <- paste0("Cannot convert tpm matrix in sparse matrix (dgCMatrix):", e$message)
        stop(em)
      }
    )
  }

  # generate new IDs for the cells
  n_cells <- dim(count_matrix)[2]
  cells_old <- colnames(count_matrix)
  new_ids <- paste0(name, "_", seq_len(n_cells))

  #### build annotation table ####
  anno_df <- data.frame(
    cell_ID = new_ids,
    cell_ID.old = cells_old,
    cell_type = annotation[["cell_type"]]
  )


  # add additional column with name "spike_in" if this data is available
  if (!is.null(spike_in_col)) {
    if (!spike_in_col %in% colnames(annotation)) {
      stop("Could not find spike_in column in annotation file.")
    }
    anno_df <- cbind(anno_df, spike_in = annotation[[spike_in_col]])
  }
  # add all additional columns mentioned in "additional_cols"
  if (!is.null(additional_cols)) {
    if (!any(additional_cols %in% colnames(annotation))) {
      stop("Not all columns mentioned in additional_cols are present in the provided annotation file.")
    }
    if (any(colnames(anno_df) %in% additional_cols)) {
      stop("At least one of the columns mentioned in additional_cols is already present in the dataset annotation object. Please rename.")
    }
    anno_df <- cbind(anno_df, annotation[, additional_cols])
    # rename new column(s)
    # colnames(anno_df)[ncols(anno_df)-length(additional_cols), ncols(anno_df)] <- additional_cols
  }

  # filter cells by type abundance
  if (type_abundance_cutoff > 0) {
    low_abundant_types <- names(which(table(anno_df$cell_type) < type_abundance_cutoff))
    low_abundant_cells <- anno_df[anno_df$cell_type %in% low_abundant_types, ]$cell_ID
    anno_df <- anno_df[!anno_df$cell_ID %in% low_abundant_cells, ]
  } else {
    low_abundant_cells <- NULL
  }

  #### handle matrices ####

  assays <- list()

  # filter genes
  matrices <- filter_matrix(count_matrix, tpm_matrix, filter_genes, variance_cutoff)
  count_matrix <- matrices$m1
  tpm_matrix <- matrices$m2

  # handle count matrix
  count_matrix <- compare_matrix_with_annotation(count_matrix, annotation)
  colnames(count_matrix) <- new_ids
  # remove low abundant cells
  count_matrix <- count_matrix[, which(!colnames(count_matrix) %in% low_abundant_cells)]

  assays <- append(assays, c(counts = count_matrix))

  # handle TPM matrix (if present)
  if (!is.null(tpm_matrix)) {
    if (!check_if_tpm(tpm_matrix)) {
      message("Some cells in your TPM matrix are not scaled between 7e5 and 1e6, as it would be expected for TPM data. Moving on..")
    }
    tpm_matrix <- compare_matrix_with_annotation(tpm_matrix, annotation)
    colnames(tpm_matrix) <- new_ids
    # remove low abundant cells
    tpm_matrix <- tpm_matrix[, which(!colnames(tpm_matrix) %in% low_abundant_cells)]
    if (scale_tpm) {
      tpm_matrix <- cpm_normalize(tpm_matrix)
    }

    assays <- append(assays, c(tpm = tpm_matrix))
  }

  # add additional columns to annotation based on count matrix [nReads, nGenes]
  # add number of reads per cell
  anno_df <- cbind(anno_df, nReads_SimBu = Matrix::colSums(count_matrix))
  # add number of expressed genes per cell (number of genes - number of genes with 0 expression)
  anno_df <- cbind(anno_df, nGenes_SimBu = (nrow(count_matrix) - proxyC::colZeros(count_matrix)))

  # create the actual SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = anno_df
  )

  message("Created dataset.")

  return(se)
}




#' Build \link[SummarizedExperiment]{SummarizedExperiment} using local annotation and count matrix R objects
#'
#' @param annotation (mandatory) dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param count_matrix (mandatory) sparse count matrix; raw count data is expected with genes in rows, cells in columns
#' @param tpm_matrix sparse count matrix; TPM like count data is expected with genes in rows, cells in columns
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param spike_in_col which column in annotation contains information on spike_in counts, which can be used to re-scale counts; mandatory for spike_in scaling factor in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff (default = 0)
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return Return a \link[SummarizedExperiment]{SummarizedExperiment} object
#' @export
#'
#' @examples
#'
#' counts <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
#' tpm <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
#' tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))
#'
#' colnames(counts) <- paste0("cell_", rep(1:300))
#' colnames(tpm) <- paste0("cell_", rep(1:300))
#' rownames(counts) <- paste0("gene_", rep(1:1000))
#' rownames(tpm) <- paste0("gene_", rep(1:1000))
#'
#' annotation <- data.frame(
#'   "ID" = paste0("cell_", rep(1:300)),
#'   "cell_type" = c(rep("T cells CD4", 300))
#' )
#'
#' ds <- SimBu::dataset(annotation = annotation, count_matrix = counts, tpm_matrix = tpm, name = "test_dataset")
#'
dataset <- function(annotation, count_matrix = NULL, tpm_matrix = NULL, name = "SimBu_dataset",
                    spike_in_col = NULL, additional_cols = NULL, filter_genes = TRUE, variance_cutoff = 0, type_abundance_cutoff = 0, scale_tpm = TRUE) {
  generate_summarized_experiment(
    annotation = annotation,
    count_matrix = count_matrix,
    tpm_matrix = tpm_matrix,
    name = name,
    spike_in_col = spike_in_col,
    additional_cols = additional_cols,
    filter_genes = filter_genes,
    variance_cutoff = variance_cutoff,
    type_abundance_cutoff = type_abundance_cutoff,
    scale_tpm = scale_tpm
  )
}

#' Merge multiple \link[SummarizedExperiment]{SummarizedExperiment} datasets into one
#'
#' The objects need to have the same number of assays in order to work.
#'
#' @param dataset_list  (mandatory) list of \link[SummarizedExperiment]{SummarizedExperiment} objects
#' @param name name of the new dataset
#' @param spike_in_col which column in annotation contains information on spike_in counts, which can be used to re-scale counts; mandatory for spike_in scaling factor in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return \link[SummarizedExperiment]{SummarizedExperiment} object
#' @export
#'
#' @examples
#'
#' counts <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
#' tpm <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
#' tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))
#'
#' colnames(counts) <- paste0("cell_", rep(1:300))
#' colnames(tpm) <- paste0("cell_", rep(1:300))
#' rownames(counts) <- paste0("gene_", rep(1:1000))
#' rownames(tpm) <- paste0("gene_", rep(1:1000))
#'
#' annotation <- data.frame(
#'   "ID" = paste0("cell_", rep(1:300)),
#'   "cell_type" = c(rep("T cells CD4", 300))
#' )
#'
#' ds1 <- SimBu::dataset(annotation = annotation, count_matrix = counts, tpm_matrix = tpm, name = "test_dataset1")
#' ds2 <- SimBu::dataset(annotation = annotation, count_matrix = counts, tpm_matrix = tpm, name = "test_dataset2")
#' ds_merged <- SimBu::dataset_merge(list(ds1, ds2))
dataset_merge <- function(dataset_list, name = "SimBu_dataset", spike_in_col = NULL,
                          additional_cols = NULL, filter_genes = TRUE, variance_cutoff = 0, type_abundance_cutoff = 0, scale_tpm = TRUE) {
  if (length(dataset_list) <= 1) {
    stop("You need at least 2 datasets to merge them into one!")
  }

  n_assays <- unlist(lapply(dataset_list, function(x) {
    return(length(names(SummarizedExperiment::assays(x))))
  }))

  if (length(unique(n_assays)) != 1) {
    stop("The datasets you want to merge have different numbers of assays. Stopping.")
  }

  # merge SEs
  merged_se <- do.call(SummarizedExperiment::cbind, dataset_list)

  if ("counts" %in% names(SummarizedExperiment::assays(merged_se))) {
    counts <- Matrix::Matrix(SummarizedExperiment::assays(merged_se)[["counts"]], sparse = TRUE)
  } else {
    counts <- NULL
  }
  if ("tpm" %in% names(SummarizedExperiment::assays(merged_se))) {
    tpm <- Matrix::Matrix(SummarizedExperiment::assays(merged_se)[["tpm"]], sparse = TRUE)
  } else {
    tpm <- NULL
  }

  # combine all annotation dataframes to a single dataframe
  anno_df <- data.frame(SummarizedExperiment::colData(merged_se))
  anno_df$ID <- anno_df[["cell_ID"]]

  # check if there is a spike_in value for each sample; else cannot use it
  if (sum(is.na(anno_df[["spike_in"]])) > 0) {
    anno_df[["spike_in"]] <- NULL
  }

  generate_summarized_experiment(
    annotation = anno_df,
    count_matrix = counts,
    tpm_matrix = tpm,
    name = name,
    spike_in_col = spike_in_col,
    additional_cols = additional_cols,
    filter_genes = filter_genes,
    variance_cutoff = variance_cutoff,
    type_abundance_cutoff = type_abundance_cutoff,
    scale_tpm = scale_tpm
  )
}

#' Build \link[SummarizedExperiment]{SummarizedExperiment} using a h5ad file for the counts
#'
#' @param h5ad_file_counts (mandatory) h5ad file with raw count data
#' @param h5ad_file_tpm h5ad file with TPM count data
#' @param cell_id_col (mandatory) name of column in Seurat meta.data with unique cell ids; 0 for rownames
#' @param cell_type_col (mandatory) name of column in Seurat meta.data with cell type name
#' @param cells_in_obs boolean, if TRUE, cell identifiers are taken from `obs` layer in anndata object; if FALSE, they are taken from `var`
#' @param name name of the dataset; will be used for new unique IDs of cells#' @param spike_in_col which column in annotation contains information on spike_in counts, which can be used to re-scale counts; mandatory for spike_in scaling factor in simulation
#' @param spike_in_col which column in annotation contains information on spike_in counts, which can be used to re-scale counts; mandatory for spike_in scaling factor in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return Return a \link[SummarizedExperiment]{SummarizedExperiment} object
#' @export
#'
#' @examples
#' # h5 <- system.file("extdata", "anndata.h5ad", package = "SimBu")
#' # ds_h5ad <- SimBu::dataset_h5ad(
#' #  h5ad_file_counts = h5,
#' #  name = "h5ad_dataset",
#' #  cell_id_col = "id", # this will use the 'id' column of the metadata as cell identifiers
#' #  cell_type_col = "group", # this will use the 'group' column of the metadata as cell type info
#' #  cells_in_obs = TRUE
#' # ) # in case your cell information is stored in the var layer, switch to FALSE
dataset_h5ad <- function(h5ad_file_counts, h5ad_file_tpm = NULL, cell_id_col = "ID", cell_type_col = "cell_type",
                         cells_in_obs = TRUE, name = "SimBu_dataset", spike_in_col = NULL, additional_cols = NULL, filter_genes = TRUE, variance_cutoff = 0, type_abundance_cutoff = 0, scale_tpm = TRUE) {
  if (all(is.null(c(h5ad_file_counts, h5ad_file_tpm)))) {
    stop("You need to provide at least one h5ad file.")
  }

  if (!is.null(h5ad_file_counts) && !file.exists(h5ad_file_counts)) {
    stop("Incorrect path to counts file given; file does not exist.")
  } else if (!is.null(h5ad_file_counts)) {
    h5ad_file_counts <- normalizePath(h5ad_file_counts)

    file_type <- tools::file_ext(h5ad_file_counts)
    if (file_type == "h5ad") {
      h5ad_data <- h5ad_to_adata(h5ad_file_counts, cells_in_obs)
      count_matrix <- h5ad_data$mm
      anno_counts <- h5ad_data$anno

      # find cell id information in annotation dataframe
      if (!cell_id_col %in% colnames(anno_counts)) {
        # use rownames as ID if wished
        if (cell_id_col == 0) {
          anno_counts$ID <- rownames(anno_counts)
          cell_id_col <- "ID"
        } else {
          em <- paste0('Cannot find "', cell_id_col, '" column in cell annotation of h5ad_file_counts.')
          stop(em)
        }
      }

      # find cell type information in annotation dataframe
      if (!cell_type_col %in% colnames(anno_counts)) {
        em <- paste0('Cannot find "', cell_type_col, '" column in cell annotation of h5ad_file_counts.')
        stop(em)
      }
      colnames(anno_counts)[which(colnames(anno_counts) == cell_id_col)] <- "ID"
      colnames(anno_counts)[which(colnames(anno_counts) == cell_type_col)] <- "cell_type"
    } else {
      stop("No valid file type; only h5ad is permitted")
    }
  } else {
    count_matrix <- NULL
    anno_counts <- NULL
  }

  if (!is.null(h5ad_file_tpm) && !file.exists(h5ad_file_tpm)) {
    stop("Incorrect path to counts file given; file does not exist.")
  } else if (!is.null(h5ad_file_tpm)) {
    h5ad_file_tpm <- normalizePath(h5ad_file_tpm)

    file_type <- tools::file_ext(h5ad_file_tpm)
    if (file_type == "h5ad") {
      h5ad_data <- h5ad_to_adata(h5ad_file_tpm, cells_in_obs)
      tpm_matrix <- h5ad_data$mm
      anno_tpm <- h5ad_data$anno

      # find cell id information in annotation dataframe
      if (!cell_id_col %in% colnames(anno_tpm)) {
        # use rownames as ID if wished
        if (cell_id_col == 0) {
          anno_tpm$ID <- rownames(anno_tpm)
          cell_id_col <- "ID"
        } else {
          em <- paste0('Cannot find "', cell_id_col, '" column in cell annotation of h5ad_file_tpm.')
          stop(em)
        }
      }

      # find cell type information in annotation dataframe
      if (!cell_type_col %in% colnames(anno_tpm)) {
        em <- paste0('Cannot find "', cell_type_col, '" column in cell annotation of h5ad_file_tpm.')
        stop(em)
      }
      colnames(anno_tpm)[which(colnames(anno_tpm) == cell_id_col)] <- "ID"
      colnames(anno_tpm)[which(colnames(anno_tpm) == cell_type_col)] <- "cell_type"
    } else {
      stop("No valid file type; only h5ad is permitted")
    }
  } else {
    tpm_matrix <- NULL
    anno_tpm <- NULL
  }

  # check if both annotations are equal
  if (!is.null(anno_tpm)) {
    if (!all(anno_counts[["ID"]] %in% anno_tpm[["ID"]])) {
      warning("Not all cells are present in counts and tpm annotation. Intersection will be used.")
      intersect_cells <- Reduce(intersect, list(anno_counts[["ID"]], anno_tpm[["ID"]]))
      annotation <- anno_counts[anno_counts[["ID"]] %in% intersect_cells, ]
    } else {
      annotation <- anno_counts
    }
  } else {
    annotation <- anno_counts
  }

  generate_summarized_experiment(
    annotation = annotation,
    count_matrix = count_matrix,
    tpm_matrix = tpm_matrix,
    name = name,
    spike_in_col = spike_in_col,
    additional_cols = additional_cols,
    filter_genes = filter_genes,
    variance_cutoff = variance_cutoff,
    type_abundance_cutoff = type_abundance_cutoff,
    scale_tpm = scale_tpm
  )
}

#' Build \link[SummarizedExperiment]{SummarizedExperiment} using a \link[Seurat]{Seurat} object
#'
#' @param seurat_obj (mandatory) \link[Seurat]{Seurat} object with TPM counts
#' @param count_assay (mandatory) name of assay in Seurat object which contains count data in 'counts' slot
#' @param cell_id_col (mandatory) name of column in Seurat meta.data with unique cell ids
#' @param cell_type_col (mandatory) name of column in Seurat meta.data with cell type name
#' @param tpm_assay name of assay in Seurat object which contains TPM data in 'counts' slot
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param spike_in_col which column in annotation contains information on spike_in counts, which can be used to re-scale counts; mandatory for spike_in scaling factor in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return Return a \link[SummarizedExperiment]{SummarizedExperiment} object
#' @export
#'
#' @examples
#' counts <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
#' tpm <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
#' tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))
#'
#' colnames(counts) <- paste0("cell-", rep(1:300))
#' colnames(tpm) <- paste0("cell-", rep(1:300))
#' rownames(counts) <- paste0("gene-", rep(1:1000))
#' rownames(tpm) <- paste0("gene-", rep(1:1000))
#'
#' annotation <- data.frame(
#'   "ID" = paste0("cell-", rep(1:300)),
#'   "cell_type" = c(
#'     rep("T cells CD4", 50),
#'     rep("T cells CD8", 50),
#'     rep("Macrophages", 100),
#'     rep("NK cells", 10),
#'     rep("B cells", 70),
#'     rep("Monocytes", 20)
#'   ),
#'   row.names = paste0("cell-", rep(1:300))
#' )
#'
#' seurat_obj <- Seurat::CreateSeuratObject(counts = counts, assay = "counts", meta.data = annotation)
#' tpm_assay <- Seurat::CreateAssayObject(counts = tpm)
#' seurat_obj[["tpm"]] <- tpm_assay
#'
#' ds_seurat <- SimBu::dataset_seurat(
#'   seurat_obj = seurat_obj,
#'   count_assay = "counts",
#'   cell_id_col = "ID",
#'   cell_type_col = "cell_type",
#'   tpm_assay = "tpm",
#'   name = "seurat_dataset"
#' )
dataset_seurat <- function(seurat_obj, count_assay, cell_id_col, cell_type_col, tpm_assay = NULL, name = "SimBu_dataset",
                           spike_in_col = NULL, additional_cols = NULL, filter_genes = TRUE, variance_cutoff = 0, type_abundance_cutoff = 0, scale_tpm = TRUE) {
  if (is.null(count_assay)) {
    stop("You have to provide the name of the assay in the Seurat object which contains count data.")
  }

  if (!count_assay %in% names(seurat_obj@assays)) {
    stop("The provided count_assay name was not found in the Seurat object.")
  }

  if (!is.null(tpm_assay)) {
    if (!count_assay %in% names(seurat_obj@assays)) {
      stop("The provided count_assay name was not found in the Seurat object.")
    }
  }

  annotation <- seurat_obj@meta.data
  if (is.null(annotation)) {
    stop("No annotation found in the Seurat object.")
  }

  if (!cell_id_col %in% colnames(annotation)) {
    stop("Did not find cell_id_col in Seurat metadata.")
  }

  if (!cell_type_col %in% colnames(annotation)) {
    stop("Did not find cell_type_col in Seurat metadata.")
  }


  tryCatch(
    {
      count_matrix <- seurat_obj@assays[[count_assay]]@counts
    },
    error = function(e) {
      em <- paste("Could not access count matrix from Seurat object (counts): ", e)
      stop(em)
      return(NULL)
    }
  )


  if (!is.null(tpm_assay)) {
    tryCatch(
      {
        tpm_matrix <- seurat_obj@assays[[tpm_assay]]@counts
      },
      error = function(e) {
        em <- paste("Could not access count matrix from Seurat object (tpm): ", e)
        stop(em)
        return(NULL)
      }
    )
  } else {
    tpm_matrix <- NULL
  }

  generate_summarized_experiment(
    annotation = annotation,
    count_matrix = count_matrix,
    tpm_matrix = tpm_matrix,
    name = name,
    spike_in_col = spike_in_col,
    additional_cols = additional_cols,
    filter_genes = filter_genes,
    variance_cutoff = variance_cutoff,
    type_abundance_cutoff = type_abundance_cutoff,
    scale_tpm = scale_tpm
  )
}

#' Build \link[SummarizedExperiment]{SummarizedExperiment} using a single sfaira entry ID
#'
#' @param sfaira_id (mandatory) ID of a sfaira dataset
#' @param sfaira_setup (mandatory) the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param spike_in_col which column in annotation contains information on spike_in counts, which can be used to re-scale counts
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param force boolean, if TRUE, datasets without annotation will be downloaded, FALSE otherwise (default)
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return dataset object
#' @export
#'
#' @examples
#' \donttest{
#' setup_list <- SimBu::setup_sfaira(tempdir())
#' ds <- SimBu::dataset_sfaira(
#'   sfaira_id = "homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x",
#'   sfaira_setup = setup_list,
#'   name = "test_dataset"
#' )
#' }
dataset_sfaira <- function(sfaira_id, sfaira_setup, name = "SimBu_dataset",
                           spike_in_col = NULL, additional_cols = NULL, force = FALSE, filter_genes = TRUE, variance_cutoff = 0, type_abundance_cutoff = 0, scale_tpm = TRUE) {
  if (is.null(sfaira_setup)) {
    warning("You need to setup sfaira first; please use setup_sfaira() to do so.")
    return(NULL)
  }
  message("Starting to download dataset from Sfaria with id: ")
  message(sfaira_id)
  sfaira_data <- download_sfaira(setup_list = sfaira_setup, ids = sfaira_id, force = force)
  count_matrix <- Matrix::t(sfaira_data$X)
  if (!is.null(sfaira_data$var$gene_symbol)) {
    rownames(count_matrix) <- sfaira_data$var$gene_symbol
  }
  annotation <- check_annotation(sfaira_data$obs)
  if (is.null(colnames(count_matrix)) && dim(count_matrix)[2] == dim(annotation)[1]) {
    colnames(count_matrix) <- annotation[["ID"]]
  }
  if (is.null(colnames(count_matrix)) && dim(count_matrix)[2] != dim(annotation)[1]) {
    warning("The count matrix has no column names and has different dimensions than the annotation. Cannot load this as a dataset.")
    return(NULL)
  }

  generate_summarized_experiment(
    annotation = annotation,
    count_matrix = count_matrix,
    tpm_matrix = NULL,
    name = name,
    spike_in_col = spike_in_col,
    additional_cols = additional_cols,
    filter_genes = filter_genes,
    variance_cutoff = variance_cutoff,
    type_abundance_cutoff = type_abundance_cutoff,
    scale_tpm = scale_tpm
  )
}


#' Build \link[SummarizedExperiment]{SummarizedExperiment} using multiple sfaira entries
#'
#' You can apply different filters on the whole data-zoo of sfaria; the resulting single-cell datasets will
#' be combined into a single dataset which you can use for simulation
#' Note: only datasets in sfaira with annotation are considered!
#'
#' @param organisms (mandatory) list of organisms (only human and mouse available)
#' @param tissues (mandatory) list of tissues
#' @param assays (mandatory) list of assays
#' @param sfaira_setup (mandatory) the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param spike_in_col which column in annotation contains information on spike_in counts, which can be used to re-scale counts
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return dataset object
#' @export
#'
#' @examples
#' \donttest{
#' setup_list <- SimBu::setup_sfaira(tempdir())
#' ds_human_lung <- SimBu::dataset_sfaira_multiple(
#'   sfaira_setup = setup_list,
#'   organisms = "Homo sapiens",
#'   tissues = "lung parenchyma",
#'   assay = "10x 3' v2",
#'   name = "human_lung"
#' )
#' }
#'
dataset_sfaira_multiple <- function(organisms = NULL, tissues = NULL, assays = NULL, sfaira_setup, name = "SimBu_dataset",
                                    spike_in_col = NULL, additional_cols = NULL, filter_genes = TRUE, variance_cutoff = 0, type_abundance_cutoff = 0, scale_tpm = TRUE) {
  if (is.null(sfaira_setup)) {
    warning("You need to setup sfaira first; please use setup_sfaira() to do so.")
    return(NULL)
  }
  message("Starting to download datasets from Sfaria...")
  sfaira_data <- download_sfaira_multiple(sfaira_setup, organisms, tissues, assays)
  count_matrix <- Matrix::t(sfaira_data$X)
  if (!is.null(sfaira_data$var$gene_symbol)) {
    rownames(count_matrix) <- sfaira_data$var$gene_symbol
  }
  annotation <- check_annotation(sfaira_data$obs)
  if (is.null(colnames(count_matrix)) && dim(count_matrix)[2] == dim(annotation)[1]) {
    colnames(count_matrix) <- annotation[["ID"]]
  }
  if (is.null(colnames(count_matrix)) && dim(count_matrix)[2] != dim(annotation)[1]) {
    warning("The count matrix has no column names and has different dimensions than the annotation. Cannot load this as a dataset.")
    return(NULL)
  }
  generate_summarized_experiment(
    annotation = annotation,
    count_matrix = count_matrix,
    tpm_matrix = NULL,
    name = name,
    spike_in_col = spike_in_col,
    additional_cols = additional_cols,
    filter_genes = filter_genes,
    variance_cutoff = variance_cutoff,
    type_abundance_cutoff = type_abundance_cutoff,
    scale_tpm = scale_tpm
  )
}


#' Use basilisk environment to read h5ad file and access anndata object
#'
#' @param h5ad_path path to h5ad file
#' @param cells_in_obs boolean, if TRUE, cell identifiers are taken from `obs` layer in anndata object; if FALSE, they are taken from `var`
#'
#' @return matrix contained on h5ad file as dgCMatrix
#'
#' @keywords internal
h5ad_to_adata <- function(h5ad_path, cells_in_obs) {
  h5ad_path <- normalizePath(h5ad_path)

  # create conda environment with anndata
  proc <- basilisk::basiliskStart(SimBu_env)
  on.exit(basilisk::basiliskStop(proc))

  tryCatch(
    {
      # initialize environment
      h5ad_data <- basilisk::basiliskRun(proc, function(h5ad_path, cells_in_obs) {
        sp <- reticulate::import("scanpy")
        adata <- sp$read_h5ad(h5ad_path)
        if (!cells_in_obs) {
          adata <- adata[["T"]]
        }
        mm <- Matrix::t(methods::as(methods::as(adata[["X"]], "CsparseMatrix"), "dgCMatrix"))
        colnames(mm) <- rownames(data.frame(adata[["obs"]]))
        rownames(mm) <- rownames(data.frame(adata[["var"]]))
        return(list(
          mm = mm,
          anno = data.frame(adata[["obs"]], stringsAsFactors = FALSE)
        ))
      }, h5ad_path = h5ad_path, cells_in_obs = cells_in_obs)
      return(h5ad_data)
    },
    error = function(e) {
      message("Could not access h5ad file: ", h5ad_path)
      return(NULL)
    }
  )
}


#' check for correct column names in annotation file and replace them if necessary
#'
#' @param annotation dataframe; annotation dataframe
#' @param cell_column name of cell-type column; default is "cell_type"
#' @param id_column name of cell ID column; default is 1, which uses the rownames
#'
#' @return annotation dataframe with correct column names
#'
#' @keywords internal
check_annotation <- function(annotation, cell_column = "cell_type", id_column = 1) {
  # check the ID column
  if (id_column == 1) {
    message("Using rownames for cell-IDs.")
    annotation$ID <- rownames(annotation)
  } else if (!"ID" %in% colnames(annotation)) {
    message("No \'ID\' column in the annotation file. Will use the supplied id_column name.")
    if (!id_column %in% colnames(annotation)) {
      warning("Supplied id_column name does not exist in annotation. Possible column names are:")
      message(colnames(annotation))
      return(NULL)
    } else {
      colnames(annotation)[which(colnames(annotation) == id_column)] <- "ID"
    }
  }

  # check the cell_type column
  if (!"cell_type" %in% colnames(annotation)) {
    message("No \'cell_type\' column in the annotation file. Will use the supplied cell_column name.")
    if (!cell_column %in% colnames(annotation)) {
      warning("Supplied cell_column name does not exist in annotation. Possible column names are:")
      message(colnames(annotation))
      return(NULL)
    } else {
      colnames(annotation)[which(colnames(annotation) == cell_column)] <- "cell_type"
    }
  }

  return(annotation)
}


#' Checks, if a matrix is TPM-like (columns sum up to 1e6)
#'
#' @param tpm_matrix matrix to check
#' @param lower_limit the lowest sum value, a cell may have
#'
#' @return boolean
#'
#' @keywords internal
check_if_tpm <- function(tpm_matrix, lower_limit = 7e5) {
  checks <- lapply(as.integer(Matrix::colSums(tpm_matrix)), function(x) {
    return(x <= 1e6 && x > lower_limit)
  })
  return(all(unlist(checks)))
}


#' Check if annotation and matrix have same cells
#'
#' Otherwise intersection of both is used
#'
#' @param m matrix, column names are cells
#' @param annotation data.frame, rownames are genes, cell names are in ID column
#'
#' @return intersected matrix
#'
#' @keywords internal
compare_matrix_with_annotation <- function(m, annotation) {
  cells_m <- colnames(m)
  cells_a <- annotation[["ID"]]

  # check if annotation and matrix have same cells; use intersection otherwise (apply intersection on all additional_counts as well)
  if (length(cells_a) != length(cells_m)) {
    warning("Unequal number of cells in annotation and count matrix. Intersection of both will be used!")
    cells_it <- Reduce(intersect, list(cells_a, cells_m))
    annotation <- annotation[annotation[["ID"]] %in% cells_it, ]
    m <- m[, cells_it]
    warning("Remaining number of cells:")
    warning(length(cells_it))
  }
  if (!all(cells_a %in% cells_m)) {
    stop("The cell IDs in the annotation and count_matrix do not correspond.")
  }

  return(m)
}


#' filter one (or two) expression matrix by genes
#'
#' @param m1 Matrix 1
#' @param m2 Matrix 2 (optional)
#' @param filter_genes boolean
#' @param variance_cutoff numeric, genes below this variance value are removed
#'
#' @return filtered matrix
#'
#' @keywords internal
filter_matrix <- function(m1, m2 = NULL, filter_genes = TRUE, variance_cutoff = 0) {
  genes <- rownames(m1)

  if (filter_genes) {
    message("Filtering genes...")
    # filter by expression
    low_expressed_genes_1 <- rownames(m1[which(Matrix::rowSums(m1) == 0), ])
    if (!is.null(m2)) {
      low_expressed_genes_2 <- rownames(m2[which(Matrix::rowSums(m2) == 0), ])
    } else {
      low_expressed_genes_2 <- genes
    }
    low_expressed_genes <- unlist(Reduce(intersect, list(low_expressed_genes_1, low_expressed_genes_2)))

    # filter by variance
    m1_m <- methods::as(m1, "dgCMatrix")
    if (!is.null(m2)) {
      m2_m <- methods::as(m2, "dgCMatrix")
    }
    low_variance_genes_1 <- rownames(m1_m[which(sparseMatrixStats::rowVars(m1_m) < variance_cutoff), ])
    if (!is.null(m2)) {
      low_variance_genes_2 <- rownames(m2_m[which(sparseMatrixStats::rowVars(m2_m) < variance_cutoff), ])
    } else {
      low_variance_genes_2 <- genes
    }
    low_variance_genes <- unlist(Reduce(intersect, list(low_variance_genes_1, low_variance_genes_2)))
    genes_to_keep <- genes[which(!genes %in% unique(c(low_expressed_genes, low_variance_genes)))]

    # remove low expressed and low variance genes from count matrix
    m1 <- m1[which(genes %in% genes_to_keep), ]
    if (!is.null(m2)) {
      m2 <- m2[which(genes %in% genes_to_keep), ]
    }
  }

  return(list(
    m1 = m1,
    m2 = m2
  ))
}
