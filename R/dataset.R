#' class for a dataset
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
  definition = function(.Object, annotation, count_matrix, name, count_type, spike_in_col, whitelist){
    genes <- rownames(count_matrix)
    cells_m <- colnames(count_matrix)
    cells_a <- annotation[["ID"]]

    if(length(cells_a) != length(cells_m)){
      warning("Unequal number of cells in annotation and count matrix. Intersection of both will be used!")
      cells_it <- Reduce(intersect, list(cells_a, cells_m))
      annotation <- annotation[annotation[["ID"]] %in% cells_it]
      count_matrix <- count_matrix[, cells_it]
      warning(paste0("Keeping ", length(cells_it), " cells."))
    }
    if(!all(cells_a %in% cells_m)){
      stop("The cell IDs in the annotation and count matrix do not correspond.")
    }
    # remove all cells which are not in the whitelist of cell-types from annotation & count matrix
    if(!is.null(whitelist)){
      annotation <- annotation[annotation[["cell_type"]] %in% whitelist]
      if(length(annotation) == 0){
        stop("No cells are left after using this whitelist; please check that the correct names are used.")
      }
      remaining_cells <- annotation[["ID"]]
      count_matrix <- count_matrix[,remaining_cells]
      cells_m <- colnames(count_matrix)
      genes <- rownames(count_matrix)
    }

    # generate new IDs for the cells and replace them in the count table
    new_ids <- paste0(name,"_", rep(1:length(cells_m)))
    colnames(count_matrix) <- new_ids

    # build annotation table
    anno_df <- data.frame(cell_ID = new_ids,
                          cell_ID.old = cells_m,
                          cell_type = annotation[["cell_type"]],
                          count_type = count_type)

    if(!is.null(spike_in_col)){
      anno_df <- cbind(anno_df, spike_in=annotation[[spike_in_col]])
    }

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
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param whitelist list of cell-types, which you want to use later for simulation; NULL to simulate all cell-types
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset <- function(annotation, count_matrix, name, count_type="TPM", spike_in_col=NULL, whitelist=NULL){
  new(Class="dataset",
      annotation=annotation,
      count_matrix=count_matrix,
      name=name,
      count_type=count_type,
      spike_in_col=spike_in_col,
      whitelist=whitelist
     )
}

#' constructor for a dataset using a h5ad file for the counts
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param h5ad_file h5ad file with count data
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param count_type what type of counts are in the dataset; default is 'TPM'
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param whitelist list of cell-types, which you want to use later for simulation; NULL to simulate all cell-types
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset_h5ad <- function(annotation, h5ad_file, name, count_type="TPM", spike_in_col=NULL, whitelist=NULL){

  #TODO check for valid file
  h5ad_file <- normalizePath(h5ad_file)

  file_type <- file_ext(h5ad_file)
  if(file_type == "h5ad"){
    ad <- anndata::read_h5ad(h5ad_file)
    ad <- ad$transpose()
    X_mat <- ad$X
    rownames(X_mat) <- ad$obs_names
    colnames(X_mat) <- ad$var_names
  }else{
    stop("No valid file type; only h5ad is permitted")
  }

  new(Class="dataset",
      annotation=annotation,
      count_matrix=X_mat,
      name=name,
      count_type=count_type,
      spike_in_col=spike_in_col,
      whitelist=whitelist
  )
}

#' constructor for dataset using a seurat object for the counts
#'
#' @param annotation dataframe; needs columns 'ID' and 'cell_type'; 'ID' needs to be equal with cell_names in count_matrix
#' @param seurat_obj \code{\link{[Seurat]} object; needs to have counts at position assays$RNA$counts
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param count_type what type of counts are in the dataset; default is 'TPM'
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param whitelist list of cell-types, which you want to use later for simulation; NULL to simulate all cell-types
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset_seurat <- function(annotation, seurat_obj, name, count_type ="TPM",spike_in_col=NULL, whitelist=NULL){

  tryCatch({
    count_matrix <- seurat_obj@assays$RNA@counts
  }, error=function(e){
    stop(paste("Could not access counts or annotation from Seurat object: ", e))
    return(NULL)
  })

  new(Class="dataset",
      annotation=annotation,
      count_matrix=count_matrix,
      name=name,
      count_type=count_type,
      spike_in_col=spike_in_col,
      whitelist=whitelist
  )
}

#' constructor for a dataset using a sfaira ID
#'
#' @param sfaira_id ID of a sfaira dataset
#' @param sfaira_setup the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param name name of the dataset; will be used for new unique IDs of cells
#' @param count_type what type of counts are in the dataset; default is 'TPM'
#' @param spike_in_col which column in annotation contains information on spike-in counts, which can be used to re-scale counts
#' @param whitelist list of cell-types, which you want to use later for simulation; NULL to simulate all cell-types
#' @param annotation_column name of column in annotation with cell-types; should not be needed
#' @param id_column name of column in annotation with cell IDs; should not be needed
#'
#' @return dataset object
#' @export
#'
#' @examples
dataset_sfaira <- function(sfaira_id, sfaira_setup, name, count_type ="TPM",
                           spike_in_col=NULL, whitelist=NULL,
                           annotation_column, id_column){

  if(is.null(sfaira_setup)){
    warning(paste0("You need to setup sfaira first; please use setup_sfaira() to do so."))
    return(NULL)
  }
  adata <- download_sfaira(sfaira_setup, sfaira_id)
  count_matrix <- Matrix::t(adata$X)
  annotation <- check_annotation(adata$obs, cell_column = annotation_column, id_column=id_column)

  new(Class="dataset",
      annotation=annotation,
      count_matrix=count_matrix,
      name=name,
      count_type=count_type,
      spike_in_col=spike_in_col,
      whitelist=whitelist
  )

}

#' check for correct column names in annotation file and replace them if neccesary
#'
#' @param annotation dataframe; annotation dataframe
#' @param cell_column name of cell-type column; default is "cell_type"
#' @param id_column name of cell ID column; default is 1, which uses the rownames
#'
#' @return annotation dataframe with correct column names
#' @export
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
