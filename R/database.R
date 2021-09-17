#' class for a database (multiple datasets)
#'
#' @slot counts Matrix. a sparse matrix for the counts
#' @slot annotation data.frame. a single dataframe with cell-type annotation of each cell in database
#' @slot size numeric. how many datasets are stored in the database object
#'
#' @return a database object
#' @export
#'
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
    anno_df <- rbindlist(lapply(dataset_list, function(x){x@annotation}), fill = T)

    # check if there is a spike-in value for each sample; else cannot use it
    if(sum(is.na(anno_df[["spike_in"]])) > 0){
      anno_df[["spike_in"]] <- NULL
    }

    #TODO make cell-type names uniform over all datasets

    .Object@annotation <- anno_df
    .Object@counts <- count_matrix
    .Object@size <- length(dataset_list)

    return(.Object)
  }

)


#' create a database class from multiple datasets
#'
#' @param dataset_list list of datasets
#'
#' @return database object
#' @export
#'
#' @examples
database <- function(dataset_list){
  methods::new(Class="database", dataset_list = dataset_list)
}


# add a single dataset to an existing database
add_dataset <- function(database, dataset){
  #TODO
}

merge_cell_types <- function(database){

}
