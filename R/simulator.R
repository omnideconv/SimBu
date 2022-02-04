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
#' @param data \link[SummarizedExperiment]{SummarizedExperiment} object
#' @param scaling_vector vector with scaling values for each cell; calculated by the \code{calc_scaling_vector} function
#' @param simulation_vector named vector with wanted cell-types and their fractions
#' @param total_cells numeric; number of total cells for this simulation
#' @param total_read_counts numeric; sets the total read count value for each sample
#' @param ncores numeric; number of cores used to create simulation
#'
#' @return returns two vectors (one based on counts, one based on tpm; depends on which matrices are present in data) with expression values for all genes in the provided dataset
simulate_sample <- function(data,
                            scaling_vector,
                            simulation_vector,
                            total_cells,
                            total_read_counts,
                            ncores){

  if(!all(names(simulation_vector) %in% unique(SummarizedExperiment::colData(data)[["cell_type"]]))){
    stop("Some cell-types in the provided simulation vector are not in the annotation.")
  }

  if(!all.equal(sum(simulation_vector), 1)){
    print(paste0("vector: ", unlist(simulation_vector),", sum: ", sum(simulation_vector)))
    stop("The sum of the cell-type fractions can not be larger than 1.")
  }

  # loop over all wanted cell-types and sample to have the final amount
  error_reads_total <- 0
  sampled_cells <- lapply(seq_along(simulation_vector), function(x){
    # get all cells with the current type
    cells_of_type_x <- data.frame(SummarizedExperiment::colData(data)[SummarizedExperiment::colData(data)[["cell_type"]] == names(simulation_vector[x]),])

    # how many cells of this type do we need?
    cells <- dplyr::slice_sample(cells_of_type_x, n=total_cells*simulation_vector[x], replace=T)
    cells <- cells[["cell_ID"]]

    return(cells)
  })

  # annotation
  names(sampled_cells) <- names(simulation_vector)
  simulated_annotation <- utils::stack(sampled_cells)

  # get the corresponding columns from the data
  data_sampled <- data[, unlist(sampled_cells)]

  # apply scaling vector on the sampled cells in the matrices
  scaling_vector <- scaling_vector[unlist(sampled_cells)]
  matrices <- lapply(names(SummarizedExperiment::assays(data_sampled)), function(x){})

  simulated_count_vector <- NULL
  simulated_tpm_vector <- NULL

  # apply scaling factor
  # + sum up counts/tpms of all cells per gene to get single bulk sample
  # + apply possible scaling of sample
  if("counts" %in% names(SummarizedExperiment::assays(data))){
    m <- Matrix::t(Matrix::t(SummarizedExperiment::assays(data_sampled)[["counts"]]) * scaling_vector)
    simulated_count_vector <- Matrix::rowSums(m)

    # if total_read_counts is used, rarefy/scale count data to get the desired sequencing depth
    if(!is.null(total_read_counts)){
      sum_counts <- sum(simulated_count_vector)

      if(sum_counts > total_read_counts){
        v <- phyloseq::otu_table(simulated_count_vector, taxa_are_rows = T)
        simulated_count_vector <- data.frame(phyloseq::rarefy_even_depth(physeq=otu, sample.size = 1e5, trimOTUs = F))[,1]

      }else if(sum_counts < total_read_counts){
        simulated_count_vector <- simulated_count_vector / sum(simulated_count_vector) * total_read_counts
      }
    }
  }
  if("tpm" %in% names(SummarizedExperiment::assays(data))){
    m <- Matrix::t(Matrix::t(SummarizedExperiment::assays(data_sampled)[["tpm"]]) * scaling_vector)
    simulated_tpm_vector <- Matrix::rowSums(m)
    # scale tpm based sample to 1e6
    simulated_tpm_vector <- (simulated_tpm_vector / sum(simulated_tpm_vector)) * 1e6
  }

  return(list(simulated_count_vector = simulated_count_vector,
              simulated_tpm_vector = simulated_tpm_vector))

}


#' simulate whole pseudo-bulk RNAseq dataset
#'
#' This function allows you to create a full pseudo-bulk RNAseq dataset. You need to provide a \link[SummarizedExperiment]{SummarizedExperiment} from which the cells
#' will be sampled for the simulation. Also a \code{scenario} has to be selected, where you can choose how the cells will be sampled and a
#' \code{scaling_factor} on how the read counts will be transformed proir to the simulation.
#'
#' @param data \link[SummarizedExperiment]{SummarizedExperiment} object
#' @param scenario select on of the pre-defined cell-type fraction scenarios; possible are: \code{uniform},\code{random},\code{mirror_db},\code{unique},\code{spill-over}; you can also use the \code{custom} scenario, where you need to set the \code{custom_scenario_data} parameter.
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in}, \code{read-number}, \code{custom} or \code{NONE} for no scaling factor
#' @param spike_in_cell_type name of cell-type used for \code{spike-in} scenario
#' @param spike_in_amount fraction of cell-type used for \code{spike-in} scenario; must be between \code{0} and \code{0.99}
#' @param unique_cell_type name of cell-type for \code{unique} scenario
#' @param spillover_cell_type name of cell-type used for \code{spill-over} scenario
#' @param custom_scenario_data dataframe; needs to be of size \code{nsamples} x number_of_cell_types, where each sample is a row and each entry is the cell-type fraction. Rows need to sum up to 1.
#' @param custom_scaling_vector named vector with custom scaling values for cell-types. Cell-types that do not occur in this vector but are present in the dataset will be set to 1; mandatory for \code{custom} scaling factor
#' @param balance_uniform_mirror_scenario balancing value for the \code{uniform} and \code{mirror_db} scenarios: increasing it will result in more diverse simulated fractions. To get the same fractions in each sample, set to 0. Default is 0.01.
#' @param nsamples numeric; number of samples in pseudo-bulk RNAseq dataset
#' @param ncells numeric; number of cells in each dataset
#' @param total_read_counts numeric; sets the total read count value for each sample
#' @param whitelist list; give a list of cell-types you want to keep for the simulation; if NULL, all are used
#' @param blacklist list; give a list of cell-types you want to remove for the simulation; if NULL, all are used; is applied after whitelist
#' @param ncores numeric; number of cores to use
#'
#' @return named list; \code{bulk} a \link[SummarizedExperiment]{SummarizedExperiment} object, where the assays store the simulated bulk RNAseq datasets. Can hold either one or two assays, depending on how many matrices were present in the dataset
#' \code{cell-fractions} is a dataframe with the simulated cell-fractions per sample;
#' \code{scaling_vector} scaling value for each cell in dataset
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
                          scaling_factor=c("NONE","census","spike-in", "custom", "read-number", "annotation_column"),
                          spike_in_cell_type = NULL,
                          spike_in_amount = NULL,
                          unique_cell_type = NULL,
                          custom_scenario_data = NULL,
                          custom_scaling_vector = NULL,
                          balance_uniform_mirror_scenario=0.01,
                          nsamples=100,
                          ncells=1000,
                          total_read_counts = NULL,
                          whitelist = NULL,
                          blacklist = NULL,
                          ncores = 1){

  # keep only cell-types which are in whitelist in annotation & count matrix
  if(!is.null(whitelist)){
    if(!all(whitelist %in% SummarizedExperiment::colData(data)[["cell_type"]])){
      stop("Did not find all cell-types of whitelist in annotation.")
    }
    data <- data[,SummarizedExperiment::colData(data)[SummarizedExperiment::colData(data)[["cell_type"]] %in% whitelist,][["cell_ID"]]]
    if(dim(data)[2] == 0){
      stop("No cells are left after using this whitelist; please check that the correct names are used.")
    }
  }

  # remove cell-types which are in blacklist from annotation & count matrix
  if(!is.null(blacklist)){
    if(!all(blacklist %in% SummarizedExperiment::colData(data)[["cell_type"]])){
      stop("Did not find all cell-types of blacklist in annotation.")
    }
    data <- data[,SummarizedExperiment::colData(data)[!SummarizedExperiment::colData(data)[["cell_type"]] %in% blacklist,][["cell_ID"]]]
    if(dim(data)[2] == 0){
      stop("No cells are left after using this blacklist; please check that the correct names are used.")
    }
  }

  ##### different cell-type scenarios #####

  # each existing cell-type will be appearing in equal amounts
  if(scenario == "uniform"){
    all_types <- unique(SummarizedExperiment::colData(data)[["cell_type"]])
    n_cell_types <- length(all_types)
    uniform_value <- 1/length(all_types)
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      m <- round(matrix(abs(rnorm(length(all_types), mean=1/length(all_types), sd=balance_uniform_mirror_scenario)), ncol=n_cell_types), 3)
      m <- sweep(m, 1, rowSums(m), FUN="/")
      simulation_vector <- as.vector(m[1,])
      names(simulation_vector) <- all_types
      return(simulation_vector)
    })
    # give each sample a name
    sample_names <- paste0("uniform_sample",rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # generate random cell-type fractions (depending on appearance in database)
  if(scenario == "random"){
    # generate 'nsamples' random samples
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      n_cell_types <- length(unique(SummarizedExperiment::colData(data)[["cell_type"]]))
      # generate n_cell_type amount of random fractions from the uniform distribution, which will sum up to 1
      m <- round(matrix(runif(n_cell_types, 0, 1), ncol=n_cell_types),3)
      m <- sweep(m, 1, rowSums(m), FUN="/")
      simulation_vector <- as.vector(m[1,])
      names(simulation_vector) <- unique(SummarizedExperiment::colData(data)[["cell_type"]])
      return(simulation_vector)
    })
    sample_names <- paste0("random_sample", rep(1:nsamples))
    names(simulation_vector_list) <- sample_names
  }
  # generate cell-type fractions, which mirror the fraction of each cell-type in the used dataset
  if(scenario == "mirror_db"){
    n_cell_types <- length(unique(SummarizedExperiment::colData(data)[["cell_type"]]))
    # generate 'nsamples' random samples
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      # each cell-type will be represented as many times as it occurs in the used dataset
      mirror_values <- table(SummarizedExperiment::colData(data)[["cell_type"]])/ncol(data)
      m <- unlist(lapply(mirror_values, function(y){
        return(abs(round(rnorm(1, mean=y, sd=balance_uniform_mirror_scenario),3)))
      }))
      m <- matrix(m, ncol=n_cell_types)
      m <- sweep(m, 1, rowSums(m), FUN="/")
      simulation_vector <- as.vector(m[1,])
      names(simulation_vector) <- unique(SummarizedExperiment::colData(data)[["cell_type"]])
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
    if(!spike_in_cell_type %in% unique(SummarizedExperiment::colData(data)[["cell_type"]])){
      stop("The spike-in cell-type could not be found in your dataset/database.")
    }

    random_cell_types <- setdiff(unique(SummarizedExperiment::colData(data)[["cell_type"]]), spike_in_cell_type)
    all_cell_types <- c(random_cell_types, spike_in_cell_type)

    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      n_cell_types <- length(unique(SummarizedExperiment::colData(data)[["cell_type"]]))-1
      # generate n_cell_type amount of random fractions from the uniform distribution, which will sum up to 1
      m <- matrix(runif(n_cell_types, 0, 1), ncol=n_cell_types)
      simulation_vector <- as.vector(m[1,])
      simulation_vector <- (1 - spike_in_amount) * simulation_vector/sum(simulation_vector)
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
  all_samples <- parallel::mclapply(simulation_vector_list, function(x){
    samples <- simulate_sample(data=data,
                             scaling_vector = scaling_vector,
                             simulation_vector = x,
                             total_cells = ncells,
                             total_read_counts = total_read_counts,
                             ncores=ncores)

    return(samples)
  }, mc.cores = ncores)


  if(length(names(SummarizedExperiment::assays(data))) == 2){
    bulk_counts <- Matrix(sapply(all_samples, "[[", 1), sparse=T) # 1st index in each sample is based on counts
    bulk_tpm <- Matrix(sapply(all_samples, "[[", 2), sparse=T) # 2nd index in each sample is based on TPM
    assays <- list(bulk_counts = bulk_counts, bulk_tpm = bulk_tpm)
  }else if("tpm" %in% names(SummarizedExperiment::assays(data))){
    bulk_tpm <- Matrix(sapply(all_samples, "[[", 1), sparse=T)
    assays <- list(bulk_tpm = bulk_tpm)
  }else if("counts" %in% names(SummarizedExperiment::assays(data))){
    bulk_counts <- Matrix(sapply(all_samples, "[[", 1), sparse=T)
    assays <- list(bulk_counts = bulk_counts)
  }

  # new SummarizedExperiment with bulk datasets based on counts or TPM or both
  se_bulk <- SummarizedExperiment::SummarizedExperiment(assays = assays)
  colnames(se_bulk) <- sample_names

  # remove non-unique features/genes from simulated dataset
  if(length(unique(rownames(se_bulk))) != dim(se_bulk)[1]){
    un <- unique(rownames(se_bulk))
    se_bulk <- se_bulk[un,]
  }

  # cell_fractions for all simulated samples
  cell_fractions <- data.frame(t(data.frame(simulation_vector_list)), check.names = F)

  print("Finished simulation.")

  return(list(bulk = se_bulk,
              cell_fractions = cell_fractions,
              scaling_vector = scaling_vector))

}


#' Calculate scaling factor for a dataset
#'
#' Each scaling factor has a default matrix it will try to use (counts or TPM). If the required matrix is not available, the other one is used and a warning is given.
#'
#' @param data dataset object
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in}, \code{read-number}, \code{custom} or \code{NONE} for no scaling factor
#' @param custom_scaling_vector named vector with custom scaling values for cell-types. Cell-types that do not occur in this vector but are present in the dataset will be set to 1
#' @param ncores number of cores
#'
#' @return a named vector with a scaling value for each cell in the dataset
#'
calc_scaling_vector <- function(data, scaling_factor, custom_scaling_vector, ncores){

  if(scaling_factor == "census"){
    if("tpm" %in% names(SummarizedExperiment::assays(data))){
      m <- SummarizedExperiment::assays(data)[["tpm"]]
    }else{
      warning("Scaling factor 'Census' requires TPM data, which is not available in the given dataset. Counts will be used instead.")
      m <- SummarizedExperiment::assays(data)[["counts"]]
    }
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
    if("tpm" %in% names(SummarizedExperiment::assays(data))){
      m <- SummarizedExperiment::assays(data)[["counts"]]
    }else{
      warning("Scaling factor 'spike-in' requires counts data, which is not available in the given dataset. TPMs will be used instead.")
      m <- SummarizedExperiment::assays(data)[["tpm"]]
    }
    # get subset of spike-in counts from the sampled cells
    tmp <- SummarizedExperiment::colData(data)[,c("cell_ID","spike_in","read_number")]
    #tmp <- merge(tmp, simulated_annotation, by.x="cell_ID",by.y="values", all.y=T) # need to merge to also get lines for duplicate cell entries
    scaling_vector <- (tmp$read_number - tmp$spike_in)/tmp$read_number
    names(scaling_vector) <- tmp$cell_ID

  }else if(scaling_factor == "read-number"){
    if(!"read_number" %in% colnames(data@annotation)){
      stop("The annotation in your dataset does not contain read-number information; you cannot apply the read-number scaling factor.")
    }
    tmp <- SummarizedExperiment::colData(data)[,c("cell_ID","read_number")]
    #tmp <- merge(tmp, simulated_annotation, by.x="cell_ID",by.y="values", all.y=T)
    scaling_vector <- tmp$read_number
    names(scaling_vector) <- tmp$cell_ID

  }else if (scaling_factor == "custom"){
    # needs vector with values for existing cell-types
    # cell-types that do not occur in this vector will have scaling-factor of 1
    if(is.null(custom_scaling_vector)){stop("For the custom scaling factor you need to provide a custom_scaling_vector!")}

    missing_cell_types <- as.vector(unique(SummarizedExperiment::colData(data)[["cell_type"]])[which(!unique(SummarizedExperiment::colData(data)[["cell_type"]]) %in% names(custom_scaling_vector))])
    complete_vector <- rep(1, length(missing_cell_types))
    names(complete_vector) <- missing_cell_types
    complete_vector <- data.frame(value=append(complete_vector, custom_scaling_vector), check.names=F)
    df <- merge(complete_vector, data.frame(SummarizedExperiment::colData(data)), by.x=0,by.y="cell_type", all.y=T)[,c("value","cell_ID")]
    scaling_vector <- df$value
    names(scaling_vector) <- df$cell_ID

  }else if(scaling_factor == "NONE"){
    scaling_vector <- rep(1, ncol(data))
    names(scaling_vector) <- SummarizedExperiment::colData(data)[["cell_ID"]]
  }else if(!is.null(scaling_factor)){
    message(paste0("Scaling by ", scaling_factor,"-column in annotation table; if no scaling is wished instead, use 'NONE'."))
    if(!scaling_factor %in% colnames(SummarizedExperiment::colData(data))){stop(paste0("A column with the name ", scaling_factor," cannot be found in the annotation."))}
    if(!is.numeric(SummarizedExperiment::colData(data)[,c(scaling_factor)])){stop(paste0("The column with the name ", scaling_factor," is not numeric and cannot be used for scaling."))}

    tmp <- SummarizedExperiment::colData(data)[,c("cell_ID",scaling_factor)]
    colnames(tmp) <- c("a","b")
    scaling_vector <- tmp$b
    names(scaling_vector) <- tmp$a

  }else{
    warning("No valid scaling factor method provided. Scaling all cells by 1.")
    scaling_vector <- rep(1, ncol(data))
    names(scaling_vector) <- SummarizedExperiment::colData(data)[["cell_ID"]]
  }

  return(scaling_vector)
}


# normalize samples to one million -> TPM
cpm_normalize <- function(matrix){
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
#' @return named list; \code{bulk} a \link[SummarizedExperiment]{SummarizedExperiment} object, where the assays store the simulated bulk RNAseq datasets. Can hold either one or two assays, depending on how many matrices were present in the dataset
#' \code{cell-fractions} is a dataframe with the simulated cell-fractions per sample;
#' \code{scaling_vector} scaling value for each cell in dataset
#' @export
#'
merge_simulations <- function(simulation_list){

  # merge SummarizedExperiments
  if(length(simulation_list) <= 1){
    stop("You need at least 2 simulations to merge them into one!")
  }

  n_assays <- unlist(lapply(simulation_list, function(x){
    return(length(names(SummarizedExperiment::assays(x$bulk))))
  }))

  if(length(unique(n_assays)) != 1){
    stop("The simulations you want to merge have different numbers of assays. Stopping.")
  }

  merged_se <- do.call(cbind, sapply(simulation_list, "[[", "bulk"))


  # merge cell fractions dataframe
  sample_names <- unlist(lapply(simulation_list, function(x){rownames(x[["cell_fractions"]])}))
  sample_names <- paste0(sample_names,"_",rep(1:length(sample_names)))
  cell_fractions <- data.frame(rbindlist(lapply(simulation_list, function(x){x[["cell_fractions"]]}), fill=T, use.names=T), check.names = F)
  rownames(cell_fractions) <- sample_names
  cell_fractions[is.na(cell_fractions)] <- 0


  # list of scaling vectors
  scaling_vector <- sapply(simulation_list, "[[", "scaling_vector")
  colnames(scaling_vector) <- paste0("simulation_",rep(1:length(simulation_list)))


  return(list(bulk = merged_se,
              cell_fractions = cell_fractions,
              scaling_vector = scaling_vector))
}
