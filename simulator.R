library(parallel)
library(data.table)
library(Matrix)
library(Seurat)
library(anndata)
library(tools)
library(Biobase)
source("scripts/census.R")
source("dataset.R")
source("database.R")
#library(tidyverse)


###### simulation ######

# function to sample cells according to given cell-type fractions
# returns a vector with mean count/TPM values over all selected samples per gene
# data = database or dataset object
# simulation_vector = named vector with wanted cell-types and their fractions
# total_cells = number of cells you want to sample
simulate_sample <- function(data, scaling_factor, simulation_vector, total_cells, ncores){
  
  if(!all(names(simulation_vector) %in% unique(data@annotation[["cell_type"]]))){
    stop("Some cell-types in the provided simulation vector are not in the annotation.")
  }
  
  if(sum(simulation_vector) != 1){
    stop("The cell-type fractions need to sum up to 1.")
  }
  
  # loop over all wanted cell-types and sample to have the final amount
  sampled_cells <- lapply(seq_along(simulation_vector), function(x){
    # get all cells with the current type
    cells_of_type_x <- data@annotation[data@annotation[["cell_type"]] == names(simulation_vector[x]),]
    
    # how many cells of this type do we need?
    n_cells <- ceiling(total_cells * simulation_vector[x])
    
    out <- sample(cells_of_type_x[["cell_ID"]], n_cells, replace = T)
    return(out)
  })
  
  # annotation 
  names(sampled_cells) <- names(simulation_vector)
  simulated_annotation <- stack(sampled_cells)
  
  # get the corresponding columns from the count matrix in the data
  m <- data@counts[, unlist(sampled_cells)]
  
  # apply selected scaling factor on each cell in matrix
  # this calculates a scaling vector (one value per cell) which will be applied to the matrix
  if(scaling_factor == "census"){
    census_vector <- census(m, ncores = ncores, method="monocle", expr_threshold = 0.1)
    m <- t(t(m)/census_vector)
  }else if(scaling_factor == "spike_in"){
    # if you want to transform your counts by spike-in data, an additional column in the annotation table is needed
    # with name "spike-in"; the matrix counts will then be transformed accordingly
    if(!"spike_in" %in% colnames(data@annotation)){
      stop("No column with spike-in information in annotation data. Check your dataset/database again!")
    }else{
      # get subset of spike-in counts from the sampled cells
      spike_in_vector <- data@annotation[["spike_in"]]
      names(spike_in_vector) <- data@annotation$cell_ID
      spike_in_vector <- spike_in_vector[simulated_annotation$values]
      m <- t(t(m) * spike_in_vector)
    }
  }else if(scaling_factor == "custom"){
    #TODO
  }

  # calculate the mean count/TPM value per gene to get a single pseudo-bulk sample
  simulated_count_vector <- rowMeans(m)
  
  return(simulated_count_vector = simulated_count_vector)
  
}


# function to simulate a whole pseudo-bulk dataset
# data = dataset or database object
# scenario = pre-defined cell-type fractions
# nsamples = number of pseudo-bulk samples to be generated
simulate_bulk <- function(data, 
                          scenario=c("noisy","random","spike-in","spill-over"), 
                          scaling_factor=c("NONE","census","spike-in","custom"),
                          spike_in_cell_type = NULL,
                          spike_in_amount = NULL,
                          spillover_cell_type = NULL,
                          nsamples=100, 
                          ncells=1000,
                          ncores = 1){
  
  ##### different cell-type scenarios #####
  
  # each existing cell-type will be appearing in equal amounts
  if(scenario == "noisy"){
    all_types <- unique(data@annotation[["cell_type"]])
    simulation_vector <- rep(1/length(all_types), length(all_types))
    names(simulation_vector) <- all_types
    # duplicate this nsamples amount of times
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){return(simulation_vector)})
    # give each sample a name
    sample_names <- paste0("noisy_sample",rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # generate random cell-type fractions (depending on appearence in database)
  if(scenario == "random"){
    # generate 'nsamples' random samples
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      # sample 'ncells' times from the different cell-types to get a random profile of all cell-types
      simulation_vector <- table(sample(unique(data@annotation[["cell_type"]]), size = ncells, replace = T))/ncells
      names <- names(simulation_vector)
      simulation_vector <- as.vector(simulation_vector)
      names(simulation_vector) <- names
      return(simulation_vector)
    })
    sample_names <- paste0("random_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # one cell-type will be highly over represented, the others are random
  if(scenario == "spike-in"){
    
    if(is.null(spike_in_cell_type) || is.null(spike_in_amount)){
      stop("The spike-in scenario requires you to select one cell-type which will be over represented")
    }
    if(spike_in_amount > 0.99 || spike_in_amount < 0){
      stop("The spike-in cell-type fraction needs to be between 0 and 0.99.")
    }
    if(!spike_in_cell_type %in% unique(data@annotation[["cell_type"]])){
      stop("The spike-in cell-type could not be found in your dataset/database.")
    }
    
    n_spike_in_cells <- ncells * spike_in_amount   # this is how many cells will be of type spike-in
    random_cells <- ncells - n_spike_in_cells      # this is how many cells will be of type random
    # generate random samples 
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      # sample 'random_cells' times from the different cell-types to get a random profile of all cell-types except spike-in type
      possible_cell_types <- setdiff(unique(data@annotation[["cell_type"]]), spike_in_cell_type)
      simulation_vector <- (table(sample(possible_cell_types, size = random_cells, replace = T))/random_cells)*(1-spike_in_amount)
      names <- names(simulation_vector)
      simulation_vector <- as.vector(simulation_vector)
      names(simulation_vector) <- names
      simulation_vector <-c(simulation_vector, spike_in_amount)  
      names(simulation_vector)[length(simulation_vector)] <- as.character(spike_in_cell_type)
      return(simulation_vector) 
    })
    sample_names <- paste0("spike_in_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # spill-over: only simulate a single cell-type
  if(scenario == "spill-over"){
    if(is.null(spillover_cell_type)){
      stop("The spill-over scenario requires you to select a cell-type which will be simulated")
    }
    
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      simulation_vector <- c(1)
      names(simulation_vector) <- as.character(spillover_cell_type)
      return(simulation_vector)
    })
    sample_names <- paste0("spillover_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  
  ##### generate the samples #####
  
  # sample cells and generate pseudo-bulk profiles
  idx <- 1
  cell_fractions <- data.frame(types=unique(data@annotation[["cell_type"]]))
  bulk <- do.call(cbind, lapply(simulation_vector_list, function(x){
    sample<-simulate_sample(data=data, 
                            scaling_factor = scaling_factor,
                            simulation_vector = x,
                            total_cells = ncells,
                            ncores=ncores)
    # # add a new column to cell fractions dataframe for this sample
    cell_fractions <<- suppressWarnings(merge(cell_fractions, as.data.frame(x), all=T, by.x = "types", by.y=0))
    
    progress <- 100*(round(idx/length(simulation_vector_list), digits=3))
    print(paste(progress,"%"))
    idx <<- idx+1
    
    return(sample)
  }))
  
  colnames(bulk) <- sample_names
  
  rownames(cell_fractions) <- cell_fractions[["types"]]
  cell_fractions[["types"]] <- NULL
  colnames(cell_fractions) <- sample_names
  cell_fractions <- data.frame(t(cell_fractions))

  # build bioconductor expression set
  expr_set <- ExpressionSet(assayData = bulk,
                            phenoData = new("AnnotatedDataFrame", data=cell_fractions))
  
  return(list(pseudo_bulk = bulk,
              cell_fractions = cell_fractions,
              expression_set = expr_set))
  
}


