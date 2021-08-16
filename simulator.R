###### dataset ######

# class for the dataset
# annotation_db = a single dataframe with cell-type annotation of each cell in DB
# counts_db = a sparse matrix for the counts
# name = name of the dataset
setClass("dataset", slots=list(annotation_ds="data.frame", 
                               counts_ds="Matrix", 
                               name="character"))

# constructor for dataset
setMethod(
  f="initialize",
  signature="dataset",
  definition = function(.Object, annotation, count_matrix, name, count_type){
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
                          count_type = count_type)
    
    .Object@annotation_ds <- anno_df
    .Object@counts_ds <- count_matrix
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
dataset <- function(annotation, count_matrix, name, count_type="TPM"){
  new(Class="dataset", annotation=annotation, count_matrix=count_matrix, name=name, count_type=count_type)
}


###### database ######

# class for the database (multiple datasets)
# annotation_db = a single dataframe with cell-type annotation of each cell in DB
# counts_db = a sparse matrix for the counts
# size = how many datasets are stored in the database object
setClass("database", slots=list(annotation_db="data.frame", 
                                counts_db="Matrix", 
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
    genes_list <- lapply(dataset_list, function(x){rownames(x@counts_ds)})
    genes_it <- Reduce(intersect, genes_list)
    
    # subset count matrices accordingly
    count_matrices <- lapply(dataset_list, function(x){
      x@counts_ds[genes_it, ]
    })
    
    # combine count matrices to a single one
    count_matrix <- do.call(cbind, count_matrices)
    
    # combine all annotation dataframes to a single dataframe
    anno_df <- rbindlist(lapply(dataset_list, function(x){x@annotation_ds}))
    
    #TODO make cell-type names uniform over all datasets
    
    .Object@annotation_db <- anno_df
    .Object@counts_db <- count_matrix
    .Object@size <- length(dataset_list)
    
    return(.Object)
  }
  
)


database <- function(dataset_list){
  new(Class="database", dataset_list = dataset_list)
}


# add a single dataset to an existing database
add_database <- function(database, dataset){
  #TODO
}


###### simulation ######

# function to sample cells according to given cell-type fractions
# returns a sparse count matrix and an annotation dataframe
# database = database object
# simulation_vector = named vector with wanted cell-types and their fractions
# total_cells = number of cells you want to sample
simulate <- function(database, simulation_vector, total_cells){
  
  if(!all(names(simulation_vector) %in% unique(database@annotation_db[["cell_type"]]))){
    stop("Some cell-types in the provided simulation vector are not in the annotation.")
  }
  
  if(sum(simulation_vector) != 1){
    stop("The cell-type fractions need to sum up to 1.")
  }
  
  # loop over all wanted cell-types and sample to have the final amount
  sampled_cells <- lapply(seq_along(simulation_vector), function(x){
    # get all cells with the current type
    cells_of_type_x <- database@annotation_db[database@annotation_db[["cell_type"]] == names(simulation_vector[x])]
    
    # how many cells of this type do we need?
    n_cells <- total_cells * simulation_vector[x]
    
    out <- sample(cells_of_type_x[["cell_ID"]], n_cells, replace = T)
    return(out)
  })
  
  # annotation 
  names(sampled_cells) <- names(simulation_vector)
  simulated_annotation <- stack(sampled_cells)
  
  # get the corresponding columns from the count matrix in the database
  simulated_count_matrix <- database@counts_db[, unlist(sampled_cells)]
  
  return(list(simulated_annotation = simulated_annotation,
              simulated_count_matrix = simulated_count_matrix))
  
}

#### tests #####

# tpm <- fread("/home/Data/Maynard/datasets_annotated/maynard_2020_annotated_fine/X_tpm.csv")
# mat <- Matrix(as.matrix(tpm), sparse = T)
# rm(tpm)
# genes <- fread("/home/Data/Maynard/datasets_annotated/maynard_2020_annotated_fine/var.csv")
# cells <- fread("/home/Data/Maynard/datasets_annotated/maynard_2020_annotated_fine/obs.csv")
# cellnames <- cells$Run
# genenames <- genes$symbol
# dimnames(mat)<-list(cellnames, genenames)
# mat <- t(mat)
# 
# meta_cells <- fread("/home/Data/Maynard/annotation_obs.csv")[,c("Run","cell_type")]
# colnames(meta_cells) <- c("ID", "cell_type")

ds1 <- dataset(meta_cells, mat, name = "Maynard1", count_type = "TPM")
ds2 <- dataset(meta_cells, mat, name = "Maynard2", count_type = "TPM")
ds3 <- dataset(meta_cells, mat, name = "Maynard3", count_type = "TPM")

db <- database(list(ds1,ds2,ds3))

simulation_vector <- c("Macrophage"=0.3, "T cell CD8" = 0.2, "T cell CD4" = 0.2, "NK cell" = 0.3)

sim<-simulate(db, simulation_vector,total_cells = 5000)
