require(parallel)
require(data.table)
require(Matrix)
require(Seurat)
require(anndata)
require(tools)
require(Biobase)
require(reticulate)
require(tidyr)
require(tools)
require(methods)
require(dplyr)
#source("R/scripts/census.R")
#source("R/dataset.R")
#source("R/database.R")
# needs python library sfaira


###### simulation ######

#' function to sample cells according to given cell-type fractions
#' Note: if total_read_counts is used, the cell-fractions are applied to the number of counts, not the number of cells!
#'
#' @param data \code{\link{database}} or \code{\link{dataset}} object
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in}, \code{custom}
#' @param simulation_vector named vector with wanted cell-types and their fractions
#' @param total_cells numeric; number of total cells for this simulation
#' @param total_read_counts numeric; sets the total read count value for each sample
#' @param ncores numeric; number of cores used to create simulation
#'
#' @return returns a vector with mean count/TPM values over all selected samples per gene
#' @export
#'
#' @examples
simulate_sample <- function(data, scaling_factor, simulation_vector, total_cells, total_read_counts, ncores){

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
    if(is.null(total_read_counts)){
      out <- dplyr::sample_frac(cells_of_type_x, simulation_vector[x], replace = T)
      out <- out[["cell_ID"]]
    }else{
      # fill sample with cells of current cell-type until total_read_counts value is reached
      counts_per_cell_type <- ceiling(total_read_counts * simulation_vector[x]) # total count value that this cell-type can reach max
      read_counts_current <- 0
      out <- list()
      while(read_counts_current < counts_per_cell_type){
        sampled_cell <- dplyr::sample_n(cells_of_type_x, 1)
        read_counts_current <- read_counts_current + sampled_cell[["total_counts_custom"]]
        out <- append(out, values = sampled_cell[["cell_ID"]])
      }
    }
    return(out)
  })

  # annotation
  names(sampled_cells) <- names(simulation_vector)
  simulated_annotation <- utils::stack(sampled_cells)

  # get the corresponding columns from the count matrix in the data
  m <- data@counts[, unlist(sampled_cells)]

  # apply selected scaling factor on each cell in matrix
  # this calculates a scaling vector (one value per cell) which will be applied to the matrix
  if(scaling_factor == "census"){
    # TODO: does it make sense to calc census for all bulk samples combined? or keep it like this
    census_vector <- census(m, ncores = ncores, method="monocle", expr_threshold = 1)
    m <- m*(census_vector/10e6)
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
  simulated_count_vector <- rowMeans(as.matrix(m))

  return(simulated_count_vector = simulated_count_vector)

}


# function to simulate a whole pseudo-bulk dataset
# data = dataset or database object
# scenario = pre-defined cell-type fractions
# nsamples = number of pseudo-bulk samples to be generated
#' simulate whole pseudo-bulk RNAseq dataset
#'
#' @param data \code{\link{database}} or \code{\link{dataset}} object
#' @param scenario select on of the pre-defined cell-type fraction scenarios; possible are: \code{uniform},\code{random},\code{spike-in},\code{spill-over}
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in}, \code{custom}
#' @param spike_in_cell_type name of cell-type used for \code{spike-in} scenario
#' @param spike_in_amount fraction of cell-type used for \code{spike-in} scenario; must be between \code{0} and \code{0.99}
#' @param spillover_cell_type name of cell-type used for \code{spill-over} scenario
#' @param nsamples numeric; number of samples in pseudo-bulk RNAseq dataset
#' @param ncells numeric; number of cells in each dataset
#' @param total_read_counts numeric; sets the total read count value for each sample
#' @param ncores numeric; number of cores to use
#'
#' @return named list; \code{pseudo_bulk} is a sparse matrix with the simulated counts;
#' \code{cell-fractions} is a dataframe with the simulated cell-fractions per sample;
#' \code{expression_set} is a Bioconductor Expression Set \url{http://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf}
#' @export
#'
#' @examples
simulate_bulk <- function(data,
                          scenario=c("uniform","random","spike-in","spill-over"),
                          scaling_factor=c("NONE","census","spike-in","custom"),
                          spike_in_cell_type = NULL,
                          spike_in_amount = NULL,
                          spillover_cell_type = NULL,
                          nsamples=100,
                          ncells=1000,
                          total_read_counts = NULL,
                          whitelist = NULL,
                          ncores = 1){


  ##### different cell-type scenarios #####

  # each existing cell-type will be appearing in equal amounts
  if(scenario == "uniform"){
    all_types <- unique(data@annotation[["cell_type"]])
    simulation_vector <- rep(1/length(all_types), length(all_types))
    names(simulation_vector) <- all_types
    # duplicate this nsamples amount of times
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){return(simulation_vector)})
    # give each sample a name
    sample_names <- paste0("uniform_sample",rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # generate random cell-type fractions (depending on appearance in database)
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
  bulk <- do.call(cbind, parallel::mclapply(simulation_vector_list, function(x){
    sample<-simulate_sample(data=data,
                            scaling_factor = scaling_factor,
                            simulation_vector = x,
                            total_cells = ncells,
                            total_read_counts = total_read_counts,
                            ncores=ncores)

    progress <- 100*(round(idx/length(simulation_vector_list), digits=3))
    print(paste(progress,"%"))
    idx <<- idx+1

    return(sample)
  }, mc.cores = ncores))


  colnames(bulk) <- sample_names

  cell_fractions <- data.frame(t(data.frame(simulation_vector_list)))

  # normalize count matrix to have TPM-like values
  bulk <- tpm_normalize(bulk)

  # remove non-unique features from simulated dataset
  if(length(unique(rownames(bulk))) != dim(bulk)[1]){
    un <- unique(rownames(bulk))
    bulk <- bulk[un,]
  }

  # build bioconductor expression set
  expr_set <- Biobase::ExpressionSet(assayData = bulk,
                                     phenoData = new("AnnotatedDataFrame", data=cell_fractions))

  return(list(pseudo_bulk = bulk,
              cell_fractions = cell_fractions,
              expression_set = expr_set))

}


# normalize samples to one million -> TPM
tpm_normalize <- function(matrix){
  m <- Matrix::t(1e6*Matrix::t(matrix)/Matrix::colSums(matrix))
  m <- tidyr::replace_na(m, 0)

  return(m)
}

#' Save the expression matrix of a simulated pseudo-bulk dataset to a file
#'
#' @param simulation the result of simulate_bulk()
#' @param filename the filename where to save the expression matrix to
#'
#'
#' @examples
save_simulation <- function(simulation, filename){
  write.table(exprs(simulation$expression_set), filename, quote = F, sep="\t")
}
