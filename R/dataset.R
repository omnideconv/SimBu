#' Generate SummarizedExperiment using multiple parameters
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param count_matrix sparse count matrix; raw count data is expected with genes in rows, cells in columns
#' @param tpm_matrix sparse count matrix; TPM like count data is expected with genes in rows, cells in columns
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts; mandatory for spike-in scaling factor in simulation
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return Return a \link[SummarizedExperiment]{SummarizedExperiment} object
#'
generate_summarized_experiment <- function(annotation, count_matrix, tpm_matrix, name, spike_in_col, read_number_col, additional_cols, filter_genes, variance_cutoff, type_abundance_cutoff, scale_tpm){

  if(is.null(count_matrix) && is.null(tpm_matrix)){
    stop("count_matrix and tpm_matrix are both empty. At least one is required.")
  }

  # check if both matrices (if provided) have the same dimensions, row-names & col-names
  if(!is.null(count_matrix) && !is.null(tpm_matrix)){
    if(any(dim(count_matrix) != dim(tpm_matrix))){
      stop("count_matrix and tpm_matrix have different dimensions. Please check again.")
    }
    if(!all(colnames(count_matrix) %in% colnames(tpm_matrix))){
      stop("Column names of count_matrix and tpm_matrix are not in common. Please check again.")
    }
    if(!all(rownames(count_matrix) %in% rownames(tpm_matrix))){
      stop("Row names of count_matrix and tpm_matrix are not in common. Please check again.")
    }
  }

  # generate new IDs for the cells
  n_cells <- if(!is.null(count_matrix)){dim(count_matrix)[2]}else{dim(tpm_matrix)[2]}
  cells_old <- if(!is.null(count_matrix)){colnames(count_matrix)}else{colnames(tpm_matrix)}
  new_ids <- paste0(name,"_", rep(1:n_cells))

  #### build annotation table ####
  anno_df <- data.frame(cell_ID = new_ids,
                        cell_ID.old = cells_old,
                        cell_type = annotation[["cell_type"]])

  # add additional column with name "spike_in" if this data is available
  if(!is.null(spike_in_col)){
    if(!spike_in_col %in% colnames(annotation)){
      stop("Could not find spike-in column in annotation file.")
    }
    anno_df <- cbind(anno_df, spike_in=annotation[[spike_in_col]])
  }
  # add additional column with name "read_number" if this data is available
  if(!is.null(read_number_col)){
    if(!read_number_col %in% colnames(annotation)){
      stop("Could not find read_number_col column in annotation file.")
    }
    anno_df <- cbind(anno_df, read_number=annotation[[read_number_col]])
  }
  # add all additional columns mentioned in "additional_cols"
  if(!is.null(additional_cols)){
    if(!any(additional_cols %in% colnames(annotation))){
      stop("Not all columns mentioned in additional_cols are present in the provided annotation file.")
    }
    if(any(colnames(anno_df) %in% additional_cols)){
      stop("At least one of the columns mentioned in additional_cols is already present in the dataset annotation object. Please rename.")
    }
    anno_df <- cbind(anno_df, annotation[, additional_cols])
    # rename new column(s)
    # colnames(anno_df)[ncols(anno_df)-length(additional_cols), ncols(anno_df)] <- additional_cols
  }

  # filter cells by type abundance
  if(type_abundance_cutoff > 0){
    low_abundant_types <- names(which(table(anno_df$cell_type) < type_abundance_cutoff))
    low_abundant_cells <- anno_df[anno_df$cell_type %in% low_abundant_types,]$cell_ID
    anno_df <- anno_df[!anno_df$cell_ID %in% low_abundant_cells,]
  }else{
    low_abundant_cells <- NULL
  }

  #### handle matrices ####

  assays <- list()

  # filter genes
  matrices <- filter_matrix(count_matrix, tpm_matrix, filter_genes, variance_cutoff)
  count_matrix <- matrices$m1
  tpm_matrix <- matrices$m2

  if(!is.null(count_matrix)){
    count_matrix <- compare_matrix_with_annotation(count_matrix, annotation)
    colnames(count_matrix) <- new_ids
    # remove low abundant cells
    count_matrix <- count_matrix[,which(!colnames(count_matrix) %in% low_abundant_cells)]
    anno_df$total_counts_custom <- Matrix::colSums(count_matrix)

    assays <- append(assays, c(counts = count_matrix))
  }

  if(!is.null(tpm_matrix)){
    if(!check_if_tpm(tpm_matrix)){warning("Warning: Some cells in your TPM matrix are not scaled between 7e5 and 1e6, as it would be expected for TPM data.")}
    tpm_matrix <- compare_matrix_with_annotation(tpm_matrix, annotation)
    colnames(tpm_matrix) <- new_ids
    # remove low abundant cells
    tpm_matrix <- tpm_matrix[,which(!colnames(tpm_matrix) %in% low_abundant_cells)]
    if(scale_tpm){
      tpm_matrix <- cpm_normalize(tpm_matrix)
    }

    assays <- append(assays, c(tpm = tpm_matrix))
  }

  # create the actual SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                   colData = anno_df)

  print("Created dataset.")

  return(se)
}




#' Default function to generate a \link[SummarizedExperiment]{SummarizedExperiment}
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param count_matrix sparse count matrix; raw count data is expected with genes in rows, cells in columns
#' @param tpm_matrix sparse count matrix; TPM like count data is expected with genes in rows, cells in columns
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts; mandatory for spike-in scaling factor in simulation
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return Return a \link[SummarizedExperiment]{SummarizedExperiment} object
#' @export
#'
dataset <- function(annotation, count_matrix = NULL, tpm_matrix = NULL, name = "SimBu_dataset",spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0, scale_tpm=T){

  generate_summarized_experiment(annotation=annotation,
                                 count_matrix=count_matrix,
                                 tpm_matrix=tpm_matrix,
                                 name=name,
                                 spike_in_col=spike_in_col,
                                 read_number_col=read_number_col,
                                 additional_cols=additional_cols,
                                 filter_genes=filter_genes,
                                 variance_cutoff=variance_cutoff,
                                 type_abundance_cutoff=type_abundance_cutoff,
                                 scale_tpm=scale_tpm
  )
}

#' Merge multiple \link[SummarizedExperiment]{SummarizedExperiment} datasets into one
#'
#' The objects need to have the same number of assays in order to work.
#'
#' @param dataset_list  list of \link[SummarizedExperiment]{SummarizedExperiment} objects
#' @param name name of the new dataset
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts; mandatory for spike-in scaling factor in simulation
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return \link[SummarizedExperiment]{SummarizedExperiment} object
#' @export
#'
dataset_merge <- function(dataset_list, name = "SimBu_dataset", spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL,filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0, scale_tpm=T){
  if(length(dataset_list) <= 1){
    stop("You need at least 2 datasets to merge them into one!")
  }

  n_assays <- unlist(lapply(dataset_list, function(x){
    return(length(names(SummarizedExperiment::assays(x))))
  }))

  if(length(unique(n_assays)) != 1){
    stop("The datasets you want to merge have different numbers of assays. Stopping.")
  }

  # merge SEs
  merged_se <- do.call(cbind, dataset_list)

  if("counts" %in% names(SummarizedExperiment::assays(merged_se))){
    counts <- Matrix::Matrix(SummarizedExperiment::assays(merged_se)[["counts"]], sparse = T)
  }else{counts <- NULL}
  if("tpm" %in% names(SummarizedExperiment::assays(merged_se))){
    tpm <- Matrix::Matrix(SummarizedExperiment::assays(merged_se)[["tpm"]], sparse = T)
  }else{tpm <- NULL}

  # combine all annotation dataframes to a single dataframe
  anno_df <- data.frame(SummarizedExperiment::colData(merged_se))

  # check if there is a spike-in value for each sample; else cannot use it
  #if(sum(is.na(anno_df[["spike_in"]])) > 0){
  #  anno_df[["spike_in"]] <- NULL
  #}

  generate_summarized_experiment(annotation=anno_df,
                                 count_matrix=counts,
                                 tpm_matrix=tpm,
                                 name=name,
                                 spike_in_col=spike_in_col,
                                 read_number_col=read_number_col,
                                 additional_cols=additional_cols,
                                 filter_genes=filter_genes,
                                 variance_cutoff=variance_cutoff,
                                 type_abundance_cutoff=type_abundance_cutoff,
                                 scale_tpm=scale_tpm
  )

}

#' Function to generate a \link[SummarizedExperiment]{SummarizedExperiment} using a h5ad file for the counts
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param h5ad_file_counts h5ad file with raw count data
#' @param h5ad_file_tpm h5ad file with TPM count data
#' @param name name of the dataset; will be used for new unique IDs of cells#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts; mandatory for spike-in scaling factor in simulation
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts; mandatory for spike-in scaling factor in simulation
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return Return a \link[SummarizedExperiment]{SummarizedExperiment} object
#' @export
#'
dataset_h5ad <- function(annotation, h5ad_file_counts = NULL, h5ad_file_tpm = NULL, name = "SimBu_dataset",spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0, scale_tpm=T){

  if(all(is.null(c(h5ad_file_counts, h5ad_file_tpm)))){
    stop("You need to provide at least one h5ad file.")
  }

  if(!is.null(h5ad_file_counts) && !file.exists(h5ad_file_counts)){
    stop("Incorrect path to counts file given; file does not exist.")
  }else if(!is.null(h5ad_file_counts)){
    h5ad_file_counts <- normalizePath(h5ad_file_counts)

    file_type <- tools::file_ext(h5ad_file_counts)
    if(file_type == "h5ad"){
      ad <- anndata::read_h5ad(h5ad_file_counts)
      ad <- ad$transpose()
      count_matrix <- Matrix::Matrix(as.matrix(ad$X), sparse = T)
      rownames(count_matrix) <- ad$obs_names
      colnames(count_matrix) <- ad$var_names
    }else{
      stop("No valid file type; only h5ad is permitted")
    }
  }

  if(!is.null(h5ad_file_tpm) && !file.exists(h5ad_file_tpm)){
    stop("Incorrect path to tpm file given; file does not exist.")
  }else if(!is.null(h5ad_file_tpm)){
    h5ad_file_tpm <- normalizePath(h5ad_file_tpm)

    file_type <- tools::file_ext(h5ad_file_tpm)
    if(file_type == "h5ad"){
      ad <- anndata::read_h5ad(h5ad_file_tpm)
      ad <- ad$transpose()
      tpm_matrix <- Matrix::Matrix(as.matrix(ad$X), sparse = T)
      rownames(tpm_matrix) <- ad$obs_names
      colnames(tpm_matrix) <- ad$var_names
    }else{
      stop("No valid file type; only h5ad is permitted")
    }
  }

  generate_summarized_experiment(annotation=annotation,
                                 count_matrix=count_matrix,
                                 tpm_matrix=tpm_matrix,
                                 name=name,
                                 spike_in_col=spike_in_col,
                                 read_number_col=read_number_col,
                                 additional_cols=additional_cols,
                                 filter_genes=filter_genes,
                                 variance_cutoff=variance_cutoff,
                                 type_abundance_cutoff=type_abundance_cutoff,
                                 scale_tpm=scale_tpm
  )
}

#' Function to generate a \link[SummarizedExperiment]{SummarizedExperiment} using a \link[Seurat]{Seurat} object
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param seurat_obj_counts \link[Seurat]{Seurat} object with raw counts
#' @param seurat_obj_tpm \link[Seurat]{Seurat} object with TPM counts
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts; mandatory for spike-in scaling factor in simulation
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return Return a \link[SummarizedExperiment]{SummarizedExperiment} object
#' @export
#'
dataset_seurat <- function(annotation, seurat_obj_counts=NULL, seurat_obj_tpm=NULL, name = "SimBu_dataset",spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0, scale_tpm=T){

  if(all(is.null(c(seurat_obj_counts, seurat_obj_tpm)))){
    stop("You need to provide at least one Seurat object.")
  }

  if(!is.null(seurat_obj_counts)){
    tryCatch({
      count_matrix <- seurat_obj_counts@assays$RNA@counts
    }, error=function(e){
      stop(paste("Could not access count matrix from Seurat object (counts): ", e))
      return(NULL)
    })
  }

  if(!is.null(seurat_obj_tpm)){
    tryCatch({
      tpm_matrix <- seurat_obj_tpm@assays$RNA@counts
    }, error=function(e){
      stop(paste("Could not access count matrix from Seurat object (tpm): ", e))
      return(NULL)
    })
  }

  generate_summarized_experiment(annotation=annotation,
                                 count_matrix=count_matrix,
                                 tpm_matrix=tpm_matrix,
                                 name=name,
                                 spike_in_col=spike_in_col,
                                 read_number_col=read_number_col,
                                 additional_cols=additional_cols,
                                 filter_genes=filter_genes,
                                 variance_cutoff=variance_cutoff,
                                 type_abundance_cutoff=type_abundance_cutoff,
                                 scale_tpm=scale_tpm
  )
}

#' Build a dataset using a single sfaira entry ID
#'
#' @param sfaira_id ID of a sfaira dataset
#' @param sfaira_setup the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
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
dataset_sfaira <- function(sfaira_id, sfaira_setup, name,
                           spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, force=F, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0, scale_tpm=T){

  if(is.null(sfaira_setup)){
    warning("You need to setup sfaira first; please use setup_sfaira() to do so.")
    return(NULL)
  }
  adata <- download_sfaira(setup_list = sfaira_setup, id = sfaira_id, force = force)
  count_matrix <- Matrix::t(Matrix::as.matrix(adata$X))
  if(!is.null(adata$var$gene_symbol)){
    rownames(count_matrix) <- adata$var$gene_symbol
  }
  annotation <- check_annotation(adata$obs)
  if(is.null(colnames(count_matrix)) && dim(count_matrix)[2] == dim(annotation)[1]){
    colnames(count_matrix) <- annotation[["ID"]]
  }
  if(is.null(colnames(count_matrix)) && dim(count_matrix)[2] != dim(annotation)[1]){
    warning("The count matrix has no column names and has different dimensions than the annotation. Cannot load this as a dataset.")
    return(NULL)
  }

  generate_summarized_experiment(annotation=annotation,
                                 count_matrix=count_matrix,
                                 tpm_matrix=NULL,
                                 name=name,
                                 spike_in_col=spike_in_col,
                                 read_number_col=read_number_col,
                                 additional_cols=additional_cols,
                                 filter_genes=filter_genes,
                                 variance_cutoff=variance_cutoff,
                                 type_abundance_cutoff=type_abundance_cutoff,
                                 scale_tpm=scale_tpm
  )

}


#' Build a dataset using multiple sfaira entries
#'
#' You can apply different filters on the whole data-zoo of sfaria; the resulting single-cell datasets will
#' be combined into a single dataset which you can use for simulation
#' Note: only datasets in sfaira with annotation are considered!
#'
#' @param organisms list of organisms (only human and mouse available)
#' @param tissues list of tissues
#' @param assays list of assays
#' @param sfaira_setup the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#' @param scale_tpm boolean, if TRUE (default) the cells in tpm_matrix will be scaled to sum up to 1e6
#'
#' @return dataset object
#' @export
#'
dataset_sfaira_multiple <- function(organisms=NULL, tissues=NULL, assays=NULL, sfaira_setup, name,
                                    spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0, scale_tpm=T){
  if(is.null(sfaira_setup)){
    warning("You need to setup sfaira first; please use setup_sfaira() to do so.")
    return(NULL)
  }
  adata <- download_sfaira_multiple(sfaira_setup, organisms, tissues, assays)
  count_matrix <- Matrix::t(Matrix::as.matrix(adata$X))
  if(!is.null(adata$var$gene_symbol)){
    rownames(count_matrix) <- adata$var$gene_symbol
  }
  annotation <- check_annotation(adata$obs)
  if(is.null(colnames(count_matrix)) && dim(count_matrix)[2] == dim(annotation)[1]){
    colnames(count_matrix) <- annotation[["ID"]]
  }
  if(is.null(colnames(count_matrix)) && dim(count_matrix)[2] != dim(annotation)[1]){
    warning("The count matrix has no column names and has different dimensions than the annotation. Cannot load this as a dataset.")
    return(NULL)
  }
  generate_summarized_experiment(annotation=annotation,
                                 count_matrix=count_matrix,
                                 tpm_matrix=NULL,
                                 name=name,
                                 spike_in_col=spike_in_col,
                                 read_number_col=read_number_col,
                                 additional_cols=additional_cols,
                                 filter_genes=filter_genes,
                                 variance_cutoff=variance_cutoff,
                                 type_abundance_cutoff=type_abundance_cutoff,
                                 scale_tpm=scale_tpm
  )
}


#' check for correct column names in annotation file and replace them if neccesary
#'
#' @param annotation dataframe; annotation dataframe
#' @param cell_column name of cell-type column; default is "cell_type"
#' @param id_column name of cell ID column; default is 1, which uses the rownames
#'
#' @return annotation dataframe with correct column names
#'
check_annotation <- function(annotation, cell_column="cell_type", id_column=1){

  # check the ID column
  if(id_column == 1){
    print("Using rownames for cell-IDs.")
    annotation$ID <- rownames(annotation)
  }else if(!"ID" %in% colnames(annotation)){
    print("No \'ID\' column in the annotation file. Will use the supplied id_column name.")
    if(!id_column %in% colnames(annotation)){
      warning("Supplied id_column name does not exist in annotation. Possible column names are:")
      print(colnames(annotation))
      return(NULL)
    }else{
      colnames(annotation)[which(colnames(annotation) == id_column)] <- "ID"
    }
  }

  # check the cell_type column
  if(!"cell_type" %in% colnames(annotation)){
    print("No \'cell_type\' column in the annotation file. Will use the supplied cell_column name.")
    if(!cell_column %in% colnames(annotation)){
      warning("Supplied cell_column name does not exist in annotation. Possible column names are:")
      print(colnames(annotation))
      return(NULL)
    }else{
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
check_if_tpm <- function(tpm_matrix, lower_limit=7e5){
  checks <- lapply(Matrix::colSums(tpm_matrix), function(x){
    return(x <= 1e6 && x > lower_limit)
  })
  return(all(unlist(checks)))
}


compare_matrix_with_annotation <- function(m, annotation){
  genes_m <- rownames(m)
  cells_m <- colnames(m)
  cells_a <- annotation[["ID"]]

  # check if annotation and matrix have same cells; use intersection otherwise (apply intersection on all additional_counts as well)
  if(length(cells_a) != length(cells_m)){
    warning("Unequal number of cells in annotation and count matrix. Intersection of both will be used!")
    cells_it <- Reduce(intersect, list(cells_a, cells_m))
    annotation <- annotation[annotation[["ID"]] %in% cells_it,]
    m <- m[, cells_it]
    warning(paste0("Keeping ", length(cells_it), " cells."))
  }
  if(!all(cells_a %in% cells_m)){
    stop("The cell IDs in the annotation and count_matrix do not correspond.")
  }

  return(m)
}

# ugliest function ever
# TODO please redo this at some point...
filter_matrix <- function(m1, m2, filter_genes, variance_cutoff){

  if(!is.null(m1)){genes <- rownames(m1)}else{genes <- rownames(m2)}

  if(filter_genes){
    print("Filtering genes...")
    # filter by expression
    if(!is.null(m1)){low_expressed_genes_1 <- rownames(m1[which(Matrix::rowSums(m1) == 0),])}else{low_expressed_genes_1=genes}
    if(!is.null(m2)){low_expressed_genes_2 <- rownames(m2[which(Matrix::rowSums(m2) == 0),])}else{low_expressed_genes_2=genes}
    low_expressed_genes <- unlist(Reduce(intersect, list(low_expressed_genes_1, low_expressed_genes_2)))

    #filter by variance
    if(!is.null(m1)){m1_m <- methods::as(m1, "dgCMatrix")}
    if(!is.null(m2)){m2_m <- methods::as(m2, "dgCMatrix")}
    if(!is.null(m1)){low_variance_genes_1 <- rownames(m1_m[which(sparseMatrixStats::rowVars(m1_m) < variance_cutoff),])}else{low_variance_genes_1=genes}
    if(!is.null(m2)){low_variance_genes_2 <- rownames(m2_m[which(sparseMatrixStats::rowVars(m2_m) < variance_cutoff),])}else{low_variance_genes_2=genes}
    low_variance_genes <- unlist(Reduce(intersect, list(low_variance_genes_1, low_variance_genes_2)))
    genes_to_keep <- genes[which(!genes %in% unique(c(low_expressed_genes, low_variance_genes)))]

    # remove low expressed and low variance genes from count matrix
    if(!is.null(m1)){m1 <- m1[which(genes %in% genes_to_keep),]}
    if(!is.null(m2)){m2 <- m2[which(genes %in% genes_to_keep),]}
  }

  return(list(m1=m1,
              m2=m2))
}
