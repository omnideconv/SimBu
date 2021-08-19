###### dataset ######

# class for the dataset
# annotation_db = a single dataframe with cell-type annotation of each cell in DB
# counts_db = a sparse matrix for the counts
# name = name of the dataset
setClass("dataset", slots=list(annotation="data.frame", 
                               counts="Matrix", 
                               name="character"))

# constructor for dataset
setMethod(
  f="initialize",
  signature="dataset",
  definition = function(.Object, annotation, count_matrix, name, count_type, spike_in_col){
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
    
    # generate new IDs for the cells and replace them in the count table
    new_ids <- paste0(name,"_", rep(1:length(cells_m)))
    colnames(count_matrix) <- new_ids
    
    # build annotation table
    anno_df <- data.frame(cell_ID = new_ids, 
                          cell_ID.old = cells_m, 
                          cell_type = annotation[["cell_type"]],
                          spike_in = ifelse(!is.null(spike_in_col), annotation[[spike_in_col]], NULL),
                          count_type = count_type)
    
    .Object@annotation <- anno_df
    .Object@counts <- count_matrix
    .Object@name <- name
    return(.Object)
  }
)

# constructor for a dataset
# mat = sparse matrix; genes in rows, cells in columns 
# annotation = dataframe; needs columns ID and cell_type
#              ID needs to be equal to cell-name in mat
# name = name of dataset, will be used for unique ID of cells
# count_type = which types of counts are you providing? (TPM, gene_counts, RPKM, read_counts)
dataset <- function(annotation, count_matrix, name, count_type="TPM", spike_in_col=NULL){
  new(Class="dataset", 
      annotation=annotation, 
      count_matrix=count_matrix, 
      name=name, 
      count_type=count_type, 
      spike_in_col=spike_in_col
     )
}

# constructor for a dataset
# h5ad_file = h5ad file with count data 
# annotation = dataframe; needs columns ID and cell_type
#              ID needs to be equal to cell-name in mat
# name = name of dataset, will be used for unique ID of cells
# count_type = which types of counts are you providing? (TPM, gene_counts, RPKM, read_counts)
dataset_h5ad <- function(annotation, h5ad_file, name, count_type="TPM", spike_in_col=NULL){
  
  #TODO check for valid file
  
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
      spike_in_col=spike_in_col
  )
}
