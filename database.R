###### database ######

# class for the database (multiple datasets)
# annotation_db = a single dataframe with cell-type annotation of each cell in DB
# counts_db = a sparse matrix for the counts
# size = how many datasets are stored in the database object
setClass("database", slots=list(annotation="data.frame", 
                                counts="Matrix", 
                                size="numeric"))

# constructor for database
setMethod(
  f="initialize",
  signature = "database",
  definition = function(.Object, dataset_list){
    
    if(length(dataset_list) <= 1){
      stop("You need at least 2 datasets to create a database!")
    }
    
    # only use genes which appear in all databases
    genes_list <- lapply(dataset_list, function(x){rownames(x@counts)})
    genes_it <- Reduce(intersect, genes_list)
    
    # subset count matrices accordingly
    count_matrices <- lapply(dataset_list, function(x){
      x@counts[genes_it, ]
    })
    
    # combine count matrices to a single one
    count_matrix <- do.call(cbind, count_matrices)
    
    # combine all annotation dataframes to a single dataframe
    anno_df <- rbindlist(lapply(dataset_list, function(x){x@annotation}))
    
    #TODO make cell-type names uniform over all datasets
    
    .Object@annotation <- anno_df
    .Object@counts <- count_matrix
    .Object@size <- length(dataset_list)
    
    return(.Object)
  }
  
)


database <- function(dataset_list){
  new(Class="database", dataset_list = dataset_list)
}


# add a single dataset to an existing database
add_dataset <- function(database, dataset){
  #TODO
}