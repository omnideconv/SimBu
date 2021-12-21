#' class for a dataset
#'
#' @description
#' This is the main class to store single-cell datasets and use them for pseudo-bulk simulations.
#'
#' @slot annotation data.frame. a single dataframe with cell-type annotation of each cell in dataset
#' @slot counts Matrix. a sparse matrix for the counts
#' @slot name character. name of the dataset
#'
#' @return a datasets object
#' @export
#'
#' @examples
setClass("dataset", slots=list(annotation="data.frame",
                               counts="Matrix",
                               name="character"))

# constructor for dataset
setMethod(
  f="initialize",
  signature="dataset",
  definition = function(.Object, annotation, count_matrix, name, count_type, spike_in_col, read_number_col, additional_cols, filter_genes, variance_cutoff, type_abundance_cutoff){
    genes <- rownames(count_matrix)
    annotation <- as.data.table(annotation)
    cells_a <- annotation[["ID"]]
    cells_m <- colnames(count_matrix)

    if(length(cells_a) != length(cells_m)){
      warning("Unequal number of cells in annotation and count matrix. Intersection of both will be used!")
      cells_it <- Reduce(intersect, list(cells_a, cells_m))
      annotation <- annotation[annotation[["ID"]] %in% cells_it,]
      count_matrix <- count_matrix[, cells_it]
      warning(paste0("Keeping ", length(cells_it), " cells."))
    }
    if(!all(cells_a %in% cells_m)){
      stop("The cell IDs in the annotation and count matrix do not correspond.")
    }

    # filter expression matrix if wanted (remove 0 genes and genes with no variance)
    if(filter_genes){
      print("Filtering genes...")
      # filter by expression
      count_matrix <- count_matrix[which(Matrix::rowSums(count_matrix) > 0),]
      #filter by variance
      count_matrix <- as(count_matrix, "dgCMatrix")
      count_matrix <- count_matrix[which(sparseMatrixStats::rowVars(count_matrix) > variance_cutoff),]
    }

    # generate new IDs for the cells and replace them in the count table
    new_ids <- paste0(name,"_", rep(1:length(cells_m)))
    colnames(count_matrix) <- new_ids

    # build annotation table
    anno_df <- data.frame(cell_ID = new_ids,
                          cell_ID.old = cells_m,
                          cell_type = annotation[["cell_type"]],
                          count_type = count_type)

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
      count_matrix <- count_matrix[,which(!colnames(count_matrix) %in% low_abundant_cells)]
    }

    # add additional column with total read counts/TPMs per sample
    anno_df$total_counts_custom <- Matrix::colSums(count_matrix)
    print("Created dataset.")

    .Object@annotation <- anno_df
    .Object@counts <- count_matrix
    .Object@name <- name
    return(.Object)
  }
)

#' basic constructor for a dataset
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param count_matrix sparse matrix; genes in rows, cells in columns
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param count_type what type of counts are in the dataset; default is 'TPM'
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts; mandatory for spike-in scaling factor in simulation
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param additional_cols list of column names in annotation, that should be stored as well in dataset object
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset <- function(annotation, count_matrix, name, count_type="raw", spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0){
  methods::new(Class="dataset",
               annotation=annotation,
               count_matrix=count_matrix,
               name=name,
               count_type=count_type,
               spike_in_col=spike_in_col,
               read_number_col=read_number_col,
               additional_cols=additional_cols,
               filter_genes=filter_genes,
               variance_cutoff=variance_cutoff,
               type_abundance_cutoff=type_abundance_cutoff
  )
}

#' constructor to merge multiple datasets into one
#'
#' @param dataset_list  list of dataset objects
#' @param name name of the new dataset
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset_multiple <- function(dataset_list, name, count_type="raw", spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL,filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0){
  if(length(dataset_list) <= 1){
    stop("You need at least 2 datasets to merge them into one!")
  }

  # only use genes which appear in all databases
  genes_list <- lapply(dataset_list, function(x){rownames(x@counts)})
  genes_it <- Reduce(intersect, genes_list)
  if(length(genes_it)==0){
    stop("There are no common genes between your datasets; cannot build a database from this.")
  }

  # subset count matrices accordingly
  count_matrices <- lapply(dataset_list, function(x){
    x@counts[genes_it, ]
  })

  # combine count matrices to a single one
  count_matrix <- do.call(cbind, count_matrices)

  # combine all annotation dataframes to a single dataframe
  anno_df <- rbindlist(lapply(dataset_list, function(x){x@annotation}), fill = T)

  # check if there is a spike-in value for each sample; else cannot use it
  if(sum(is.na(anno_df[["spike_in"]])) > 0){
    anno_df[["spike_in"]] <- NULL
  }

  methods::new(Class="dataset",
               annotation=annotation,
               count_matrix=count_matrix,
               name=name,
               count_type=count_type,
               spike_in_col=spike_in_col,
               read_number_col=read_number_col,
               additional_cols=additional_cols,
               filter_genes=filter_genes,
               variance_cutoff=variance_cutoff,
               type_abundance_cutoff=type_abundance_cutoff
  )

}

#' constructor for a dataset using a h5ad file for the counts
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param h5ad_file h5ad file with count data
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param count_type what type of counts are in the dataset; default is 'TPM'
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset_h5ad <- function(annotation, h5ad_file, name, count_type="raw", spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL,filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0){

  #TODO check for valid file
  h5ad_file <- normalizePath(h5ad_file)

  file_type <- tools::file_ext(h5ad_file)
  if(file_type == "h5ad"){
    ad <- anndata::read_h5ad(h5ad_file)
    ad <- ad$transpose()
    X_mat <- ad$X
    rownames(X_mat) <- ad$obs_names
    colnames(X_mat) <- ad$var_names
  }else{
    stop("No valid file type; only h5ad is permitted")
  }

  methods::new(Class="dataset",
               annotation=annotation,
               count_matrix=X_mat,
               name=name,
               count_type=count_type,
               spike_in_col=spike_in_col,
               read_number_col=read_number_col,
               additional_cols=additional_cols,
               filter_genes=filter_genes,
               variance_cutoff=variance_cutoff,
               type_abundance_cutoff=type_abundance_cutoff
  )
}

#' constructor for dataset using a seurat object for the counts
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param seurat_obj Seurat object; needs to have counts at position assays$RNA$counts
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param count_type what type of counts are in the dataset; default is 'TPM'
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset_seurat <- function(annotation, seurat_obj, name, count_type ="raw",spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0){

  tryCatch({
    count_matrix <- seurat_obj@assays$RNA@counts
  }, error=function(e){
    stop(paste("Could not access counts or annotation from Seurat object: ", e))
    return(NULL)
  })

  methods::new(Class="dataset",
               annotation=annotation,
               count_matrix=count_matrix,
               name=name,
               count_type=count_type,
               spike_in_col=spike_in_col,
               read_number_col=read_number_col,
               additional_cols=additional_cols,
               filter_genes=filter_genes,
               variance_cutoff=variance_cutoff,
               type_abundance_cutoff=type_abundance_cutoff
  )
}

#' Build a dataset using a single sfaira entry ID
#'
#' @param sfaira_id ID of a sfaira dataset
#' @param sfaira_setup the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param count_type what type of counts are in the dataset; default is 'TPM'
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param read_number_col which column in annotation contains information on total read numbers in a cell; mandatory for "spike-in" scaling factor and "read-number" in simulation
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset_sfaira <- function(sfaira_id, sfaira_setup, name, count_type ="raw",
                           spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, force=F, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0){

  if(is.null(sfaira_setup)){
    warning("You need to setup sfaira first; please use setup_sfaira() to do so.")
    return(NULL)
  }
  adata <- download_sfaira(setup_list = sfaira_setup, id = sfaira_id, force = force)
  count_matrix <- Matrix::t(adata$X)
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

  methods::new(Class="dataset",
               annotation=annotation,
               count_matrix=count_matrix,
               name=name,
               count_type=count_type,
               spike_in_col=spike_in_col,
               read_number_col=read_number_col,
               additional_cols=additional_cols,
               filter_genes=filter_genes,
               variance_cutoff=variance_cutoff,
               type_abundance_cutoff=type_abundance_cutoff
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
#' @param count_type what type of counts are in the dataset; default is 'TPM'
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param filter_genes boolean, if TRUE, removes all genes with 0 expression over all samples & genes with variance below \code{variance_cutoff}
#' @param variance_cutoff numeric, is only applied if \code{filter_genes} is TRUE: removes all genes with variance below the chosen cutoff
#' @param type_abundance_cutoff numeric, remove all cells, whose cell-type appears less then the given value. This removes low abundant cell-types
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset_sfaira_multiple <- function(organisms=NULL, tissues=NULL, assays=NULL, sfaira_setup, name,
                                    count_type ="raw",spike_in_col=NULL, read_number_col=NULL, additional_cols=NULL, filter_genes=T, variance_cutoff=0, type_abundance_cutoff=0){
  if(is.null(sfaira_setup)){
    warning("You need to setup sfaira first; please use setup_sfaira() to do so.")
    return(NULL)
  }
  adata <- download_sfaira_multiple(sfaira_setup, organisms, tissues, assays)
  count_matrix <- Matrix::t(adata$X)
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
  methods::new(Class="dataset",
               annotation=annotation,
               count_matrix=count_matrix,
               name=name,
               count_type=count_type,
               spike_in_col=spike_in_col,
               read_number_col=read_number_col,
               additional_cols=additional_cols,
               filter_genes=filter_genes,
               variance_cutoff=variance_cutoff,
               type_abundance_cutoff=type_abundance_cutoff
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
#' @examples
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
