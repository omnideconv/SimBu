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
require(sparseMatrixStats)
require(ggplot2)
require(scales)
require(matrixStats)


###### simulation ######

#' simulate single pseudo-bulk sample
#'
#' function to sample cells according to given cell-type fractions. This creates a single pseudo-bulk sample by calculating the
#' mean expression value per gene over all sampled cells.
#' Note: if total_read_counts is used, the cell-fractions are applied to the number of counts, not the number of cells!
#' @param data \code{\link{database}} or \code{\link{dataset}} object
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in}, \code{custom}
#' @param scaling_vector vector with scaling values for each cell; calculated by the \code{calc_scaling_vector} function
#' @param simulation_vector named vector with wanted cell-types and their fractions
#' @param scaling_factor_aggregation aggregation function how to apply the scaling factor, possible are: \code{multiply}(default) and \code{division}
#' @param sample_aggreation aggregation method on how to generate a sample; possible are: \code{mean}(default), \code{sum}, \code{median}
#' @param total_cells numeric; number of total cells for this simulation
#' @param total_read_counts numeric; sets the total read count value for each sample
#' @param ncores numeric; number of cores used to create simulation
#'
#' @return returns a vector with expression values for all genes in the provided dataset
#' @export
simulate_sample <- function(data,
                            scaling_factor,
                            scaling_vector,
                            simulation_vector,
                            scaling_factor_aggregation,
                            sample_aggreation,
                            total_cells,
                            total_read_counts, ncores){

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
      out <- dplyr::slice_sample(cells_of_type_x, n=total_cells*simulation_vector[x], replace=T)
      out <- out[["cell_ID"]]
    }else{
      # fill sample with cells of current cell-type until total_read_counts value is reached
      counts_per_cell_type <- ceiling(total_read_counts * simulation_vector[x]) # total count value that this cell-type can reach max
      read_counts_current <- 0
      out <- list()
      while(read_counts_current < counts_per_cell_type){
        sampled_cell <- dplyr::slice_sample(cells_of_type_x, n=1)
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

  # apply scaling vector on the sampled cells in the count matrix
  scaling_vector <- scaling_vector[unlist(sampled_cells)]
  switch (scaling_factor_aggregation,
    "multiply" = m <- t(t(m) * scaling_vector),
    "division" = m <- t(t(m) / scaling_vector)
  )

  # calculate the mean expression value per gene to get a single pseudo-bulk sample
  switch (sample_aggreation,
    "mean" = simulated_count_vector <- rowMeans(as.matrix(m)),
    "sum" = simulated_count_vector <- rowSums(as.matrix(m)),
    "median" = simulated_count_vector <- matrixStats::rowMedians(as.matrix(m))
  )


  return(simulated_count_vector = simulated_count_vector)

}


#' simulate whole pseudo-bulk RNAseq dataset
#'
#' This function allows you to create a full pseudo-bulk RNAseq dataset. You need to provide a \code{\link{dataset}} or \code{\link{database}} from which the cells
#' will be sampled for the simulation. Also a \code{scenario} has to be selected, where you can choose how the cells will be sampled and a
#' \code{scaling_factor} on how the read counts will be transformed proir to the simulation.
#'
#' @param data \code{\link{dataset}} object
#' @param scenario select on of the pre-defined cell-type fraction scenarios; possible are: \code{uniform},\code{random},\code{mirror_db},\code{unique},\code{spill-over}; you can also use the \code{custom} scenario, where you need to set the \code{custom_scenario_data} parameter.
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in} or \code{NONE} for no scaling factor
#' @param spike_in_cell_type name of cell-type used for \code{spike-in} scenario
#' @param spike_in_amount fraction of cell-type used for \code{spike-in} scenario; must be between \code{0} and \code{0.99}
#' @param spillover_cell_type name of cell-type used for \code{spill-over} scenario
#' @param custom_scenario_data dataframe; needs to be of size \code{nsamples} x number_of_cell_types, where each sample is a row and each entry is the cell-type fraction. Rows need to sum up to 1.
#' @param custom_scaling_vector named vector with custom scaling values for cell-types. Cell-types that do not occur in this vector but are present in the dataset will be set to 1
#' @param scaling_factor_aggregation aggregation function how to apply the scaling factor, possible are: \code{multiply}(default) and \code{division}
#' @param sample_aggreation aggregation method on how to generate a sample; possible are: \code{mean}(default), \code{sum}, \code{median}
#' @param nsamples numeric; number of samples in pseudo-bulk RNAseq dataset
#' @param ncells numeric; number of cells in each dataset
#' @param total_read_counts numeric; sets the total read count value for each sample
#' @param whitelist list; give a list of cell-types you want to keep for the simulation; if NULL, all are used
#' @param blacklist list; give a list of cell-types you want to remove for the simulation; if NULL, all are used; is applied after whitelist
#' @param ncores numeric; number of cores to use
#'
#' @return named list; \code{pseudo_bulk} is a sparse matrix with the simulated counts;
#' \code{cell-fractions} is a dataframe with the simulated cell-fractions per sample;
#' \code{scaling_vector} scaling value for each cell in dataset
#' \code{expression_set} is a Bioconductor Expression Set \url{http://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf}
#' @export
#'
#' @examples
#' \dontrun{
#' # this creates a basic dataset with uniform cell-type distribution and no additional transformation of the data with 10 samples and 2000 cells each
#' simulate_bulk(dataset, scenario="uniform", scaling_factor="NONE", nsamples=10, ncells=2000)
#'
#' # use the spike-in scenario to have 50% B cells per sample
#' simulate_bulk(dataset, scenario="spike_in", scaling_factor="NONE", nsamples=10, ncells=2000, spike_in_cell_type="Bcell", spike_in_amount=0.5)
#'
#' # use the unique scenario to only have B cells
#' simulate_bulk(dataset, scenario="unique", scaling_factor="NONE", nsamples=10, ncells=2000, unique_cell_type="Bcell")
#'
#' # simulate a dataset with custom cell-type fraction for each of the 3 samples
#' fractions <- data.frame("Bcell"=c(0.2,0.4,0.2),"Tcell"=c(0.4,0.2,0.1),"Macrophage"=c(0.4,0.4,0.7))
#' simulate_bulk(dataset, scenario="custom", scaling_factor="NONE", nsamples=3, ncells=2000, custom_scenario_data=fractions)
#'
#' # use a blacklist to exclude certain cell-types for the simulation
#' simulate_bulk(dataset, scenario="custom", scaling_factor="NONE", nsamples=3, ncells=2000, custom_scenario_data=fractions)
#' }
simulate_bulk <- function(data,
                          scenario=c("uniform","random","mirror_db","spike_in","unique", "custom"),
                          scaling_factor=c("NONE","census","spike-in", "custom"),
                          spike_in_cell_type = NULL,
                          spike_in_amount = NULL,
                          unique_cell_type = NULL,
                          custom_scenario_data = NULL,
                          custom_scaling_vector = NULL,
                          scaling_factor_aggregation = "multiply",
                          sample_aggreation = "mean",
                          nsamples=100,
                          ncells=1000,
                          total_read_counts = NULL,
                          whitelist = NULL,
                          blacklist = NULL,
                          ncores = 1){

  # keep only cell-types which are in whitelist in annotation & count matrix
  if(!is.null(whitelist)){
    if(!all(whitelist %in% data@annotation[["cell_type"]])){
      stop("Did not find all cell-types of whitelist in annotation.")
    }
    data@annotation <- data@annotation[data@annotation[["cell_type"]] %in% whitelist,]
    if(length(data@annotation) == 0){
      stop("No cells are left after using this whitelist; please check that the correct names are used.")
    }
    remaining_cells <- data@annotation[["cell_ID"]]
    data@counts <- as(data@counts[,remaining_cells], "sparseMatrix")
  }

  # remove cell-types which are in blacklist from annotation & count matrix
  if(!is.null(blacklist)){
    if(!all(blacklist %in% data@annotation[["cell_type"]])){
      stop("Did not find all cell-types of blacklist in annotation.")
    }
    data@annotation <- data@annotation[!data@annotation[["cell_type"]] %in% blacklist,]
    if(length(data@annotation) == 0){
      stop("No cells are left after using this blacklist; please check that the correct names are used.")
    }
    remaining_cells <- data@annotation[["cell_ID"]]
    data@counts <- as(data@counts[,remaining_cells], "sparseMatrix")
  }

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
      n_cell_types <- length(unique(data@annotation[["cell_type"]]))
      # generate n_cell_type amount of random fractions from the uniform distribution, which will sum up to 1
      m <- matrix(runif(n_cell_types, 0, 1), ncol=n_cell_types)
      m <- sweep(m, 1, rowSums(m), FUN="/")
      simulation_vector <- as.vector(m[1,])
      names(simulation_vector) <- unique(data@annotation[["cell_type"]])
      return(simulation_vector)
    })
    sample_names <- paste0("random_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # generate cell-type fractions, which mirror the fraction of each cell-type in the used dataset
  if(scenario == "mirror_db"){
    # generate 'nsamples' random samples
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      # each cell-type will be represented as many times as it occurs in the used dataset
      simulation_vector <- table(data@annotation$cell_type)/nrow(data@annotation)
      names <- names(simulation_vector)
      simulation_vector <- as.vector(simulation_vector)
      names(simulation_vector) <- names
      return(simulation_vector)
    })
    sample_names <- paste0("mirror_db_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # one cell-type will be highly over represented, the others are random
  if(scenario == "spike_in"){

    if(is.null(spike_in_cell_type) || is.null(spike_in_amount)){
      stop("The spike-in scenario requires you to select one cell-type which will be over represented")
    }
    if(spike_in_amount > 0.99 || spike_in_amount < 0){
      stop("The spike-in cell-type fraction needs to be between 0 and 0.99.")
    }
    if(!spike_in_cell_type %in% unique(data@annotation[["cell_type"]])){
      stop("The spike-in cell-type could not be found in your dataset/database.")
    }

    random_cell_types <- setdiff(unique(data@annotation[["cell_type"]]), spike_in_cell_type)
    all_cell_types <- c(random_cell_types, spike_in_cell_type)

    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      n_cell_types <- length(unique(data@annotation[["cell_type"]]))-1
      # generate n_cell_type amount of random fractions from the uniform distribution, which will sum up to 1
      m <- matrix(runif(n_cell_types, 0, 1), ncol=n_cell_types)
      simulation_vector <- as.vector(m[1,])
      simulation_vector <- spike_in_amount * simulation_vector/sum(simulation_vector)
      simulation_vector <- append(simulation_vector, spike_in_amount)
      names(simulation_vector) <- all_cell_types
      return(simulation_vector)
    })
    sample_names <- paste0("spike_in_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # unique: only simulate a single cell-type
  if(scenario == "unique"){
    if(is.null(unique_cell_type)){
      stop("The unique scenario requires you to select a cell-type which will be simulated")
    }
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      simulation_vector <- c(1)
      names(simulation_vector) <- as.character(unique_cell_type)
      return(simulation_vector)
    })
    sample_names <- paste0("unique_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # custom fractions
  if(scenario == "custom"){
    if(is.null(custom_scenario_data)){
      stop("You need to provide a dataframe with you custom cell-type fractions in this scenario with the custom_scenario_data parameter.")
    }
    # check dimensions and cell-types in dataframe; samples are rows
    if(nrow(custom_scenario_data) != nsamples){
      stop("The scenario data has a differnt amount of samples than you want to simulate.")
    }
    if(!all(colnames(custom_scenario_data) %in% unique(data@annotation[["cell_type"]]))){
      stop("Could not find all cell-types from scenario data in annotation.")
    }
    simulation_vector_list <- as.list(as.data.frame(t(custom_scenario_data)))
    simulation_vector_list <- lapply(simulation_vector_list, function(x){names(x)<-colnames(custom_scenario_data);x<-x})
    sample_names <- paste0("custom_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }

  ##### pre-calculate scaling factors ####

  scaling_vector <- calc_scaling_vector(data = data,
                                        scaling_factor = scaling_factor,
                                        custom_scaling_vector = custom_scaling_vector,
                                        ncores = ncores)

  ##### generate the samples #####

  # sample cells and generate pseudo-bulk profiles
  idx <- 1
  bulk <- do.call(cbind, parallel::mclapply(simulation_vector_list, function(x){
    sample<-simulate_sample(data=data,
                            scaling_factor = scaling_factor,
                            scaling_vector = scaling_vector,
                            simulation_vector = x,
                            scaling_factor_aggregation = scaling_factor_aggregation,
                            sample_aggreation = sample_aggreation,
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
              scaling_vector = scaling_vector,
              expression_set = expr_set))

}


#' Calculate scaling factor for a dataset
#'
#' @param data dataset object
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in} or \code{NONE} for no scaling factor
#' @param custom_scaling_vector named vector with custom scaling values for cell-types. Cell-types that do not occur in this vector but are present in the dataset will be set to 1
#' @param ncores number of cores
#'
#' @return a named vector with a scaling value for each cell in the dataset
#' @export
#'
calc_scaling_vector <- function(data, scaling_factor, custom_scaling_vector, ncores){

  m <- data@counts

  if(scaling_factor == "census"){
    scaling_vector <- census(m, ncores = ncores, method="monocle", expr_threshold = 0.1)
    scaling_vector <- scaling_vector/10e6
  }else if(scaling_factor == "spike-in"){
    # if you want to transform your counts by spike-in data, an additional column in the annotation table is needed
    # with name "spike-in"; the matrix counts will then be transformed accordingly
    if(!"spike_in" %in% colnames(data@annotation)){
      stop("No column with spike-in information in annotation data. Check your dataset again!")
    }
    if(!"read_number" %in% colnames(data@annotation)){
      stop("No column with total read number information in annotation data. Check your dataset again!")
    }
    # get subset of spike-in counts from the sampled cells
    tmp <- data@annotation[,c("cell_ID","spike_in","read_number")]
    #tmp <- merge(tmp, simulated_annotation, by.x="cell_ID",by.y="values", all.y=T) # need to merge to also get lines for duplicate cell entries
    scaling_vector <- (tmp$read_number - tmp$spike_in)/tmp$read_number
    names(scaling_vector) <- tmp$cell_ID

  }else if(scaling_factor == "reads"){
    if(!"read_number" %in% colnames(data@annotation)){
      stop("No column with total read number information in annotation data. Check your dataset again!")
    }
    tmp <- data@annotation[,c("cell_ID","read_number")]
    #tmp <- merge(tmp, simulated_annotation, by.x="cell_ID",by.y="values", all.y=T)
    scaling_vector <- tmp$read_number

  }else if (scaling_factor == "custom"){
    # needs vector with values for existing cell-types
    # cell-types that do not occur in this vector will have scaling-factor of 1
    if(is.null(custom_scaling_vector)){stop("For the custom scaling factor you need to provide a custom_scaling_vector!")}

    missing_cell_types <- as.vector(unique(data@annotation$cell_type)[which(!unique(data@annotation$cell_type) %in% names(custom_scaling_vector))])
    complete_vector <- rep(1, length(missing_cell_types))
    names(complete_vector) <- missing_cell_types
    complete_vector <- data.frame(value=append(complete_vector, custom_scaling_vector), check.names=F)
    scaling_vector <- merge(complete_vector, data@annotation, by.x=0,by.y="cell_type", all.y=T)[["value"]]
    names(scaling_vector) <- data@annotation$cell_ID

  }else if(scaling_factor == "NONE"){
    scaling_vector <- rep(1, nrow(data@annotation))
    names(scaling_vector) <- data@annotation$cell_ID
  }else{
    warning("No valid scaling factor method provided. Scaling all cells by 1.")
    scaling_vector <- rep(1, nrow(data@annotation))
    names(scaling_vector) <- data@annotation$cell_ID
  }

  return(scaling_vector)
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


#' Plot the cell-type fractions in your simulated dataset
#'
#' @param simulation a simulation object generated by \code{simulate_bulk}
#'
#' @return a gpplot2 barplot
#' @export
#'
plot_simulation <- function(simulation){
  fractions <- simulation$cell_fractions
  fractions$sample <- factor(rownames(fractions), levels = rownames(fractions))
  frac_long <- gather(fractions, cell_type, fraction, 1:length(fractions)-1)

  ggplot2::ggplot(frac_long, aes(x=fraction, y=sample, fill=cell_type))+
    geom_col()+
    ggtitle(paste0("Cell-type fractions for \n",nrow(fractions)," pseudo-bulk RNAseq samples"))+
    scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Paired"))(length(unique(frac_long$cell_type))))+
    theme_bw()
}


#' Combine multiple simulations into one result
#'
#' we recommend to only merge simulations from the same dataset object, otherwise the count matrices might not correspond on the gene level
#'
#' @param simulation_list a list of simulations
#'
#' @return named list; \code{pseudo_bulk} is a sparse matrix with the simulated counts;
#' \code{cell-fractions} is a dataframe with the simulated cell-fractions per sample;
#' \code{expression_set} is a Bioconductor Expression Set \url{http://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf}
#' @export
#'
merge_simulations <- function(simulation_list){
  # merge cell fractions dataframe
  sample_names <- unlist(lapply(simulation_list, function(x){rownames(x[["cell_fractions"]])}))
  sample_names <- paste0(sample_names,"_",rep(1:length(sample_names)))
  cell_fractions <- data.frame(rbindlist(lapply(simulation_list, function(x){x[["cell_fractions"]]}), fill=T, use.names=T), check.names = F)
  rownames(cell_fractions) <- sample_names
  cell_fractions[is.na(cell_fractions)] <- 0

  # merge count matrix
  intersect_genes <- Reduce(intersect, lapply(simulation_list, function(x){rownames(x[["pseudo_bulk"]])}))
  count_matrices <- lapply(simulation_list, function(x){
    m <- x[["pseudo_bulk"]]
    m[intersect_genes,]
  })
  bulk <- do.call("cbind",count_matrices)
  colnames(bulk) <- sample_names

  # build bioconductor expression set
  expr_set <- Biobase::ExpressionSet(assayData = bulk,
                                     phenoData = new("AnnotatedDataFrame", data=cell_fractions))

  return(list(pseudo_bulk = bulk,
              cell_fractions = cell_fractions,
              expression_set = expr_set))
}
