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
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in}, \code{read-number},\code{custom}
#' @param scaling_vector vector with scaling values for each cell; calculated by the \code{calc_scaling_vector} function
#' @param simulation_vector named vector with wanted cell-types and their fractions
#' @param sample_aggregation aggregation method on how to generate a sample; possible are: \code{mean}(default), \code{sum}, \code{median}
#' @param total_cells numeric; number of total cells for this simulation
#' @param total_read_counts numeric; sets the total read count value for each sample
#' @param ncores numeric; number of cores used to create simulation
#'
#' @return returns a vector with expression values for all genes in the provided dataset
simulate_sample <- function(data,
                            scaling_factor,
                            scaling_vector,
                            simulation_vector,
                            sample_aggregation,
                            total_cells,
                            total_read_counts, ncores){

  if(!all(names(simulation_vector) %in% unique(data@annotation[["cell_type"]]))){
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
    cells_of_type_x <- data@annotation[data@annotation[["cell_type"]] == names(simulation_vector[x]),]

    # how many cells of this type do we need?
    if(is.null(total_read_counts)){
      cells <- dplyr::slice_sample(cells_of_type_x, n=total_cells*simulation_vector[x], replace=T)
      cells <- cells[["cell_ID"]]
    }else{
      total_counts_x <- total_read_counts * simulation_vector[x] # total count value that this cell-type can reach max
      out <- sample_cells_by_read_depth(cells = cells_of_type_x, depth = total_counts_x)
      cells <- out$cells
      missing_reads <- out$missing_reads
      error_reads_total <<- error_reads_total + missing_reads
    }
    return(cells)
  })

  if(!is.null(total_read_counts)){
    message(paste0("Total sequencing depth of ",total_read_counts, " was missed by ", round(error_reads_total/total_read_counts, digits = 5)*100,"%."))
  }
  # annotation
  names(sampled_cells) <- names(simulation_vector)
  simulated_annotation <- utils::stack(sampled_cells)

  # get the corresponding columns from the count matrix in the data
  m <- data@counts[, unlist(sampled_cells)]

  # apply scaling vector on the sampled cells in the count matrix
  scaling_vector <- scaling_vector[unlist(sampled_cells)]
  m <- Matrix::t(Matrix::t(m) * scaling_vector)

  # calculate the mean expression value per gene to get a single pseudo-bulk sample
  switch (sample_aggregation,
    "mean" = simulated_count_vector <- rowMeans(m),
    "sum" = simulated_count_vector <- Matrix::rowSums(m),
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
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in}, \code{read-number}, \code{custom} or \code{NONE} for no scaling factor
#' @param spike_in_cell_type name of cell-type used for \code{spike-in} scenario
#' @param spike_in_amount fraction of cell-type used for \code{spike-in} scenario; must be between \code{0} and \code{0.99}
#' @param unique_cell_type name of cell-type for \code{unique} scenario
#' @param spillover_cell_type name of cell-type used for \code{spill-over} scenario
#' @param custom_scenario_data dataframe; needs to be of size \code{nsamples} x number_of_cell_types, where each sample is a row and each entry is the cell-type fraction. Rows need to sum up to 1.
#' @param custom_scaling_vector named vector with custom scaling values for cell-types. Cell-types that do not occur in this vector but are present in the dataset will be set to 1; mandatory for \code{custom} scaling factor
#' @param balance_unique_mirror_scenario balancing value for the \code{uniform} and \code{mirror_db} scenarios: increasing it will result in more diverse simulated fractions. To get the same fractions in each sample, set to 0. Default is 0.01.
#' @param sample_aggregation aggregation method on how to generate a sample; possible are: \code{sum}(default), \code{mean}, \code{median}
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
                          scaling_factor=c("NONE","census","spike-in", "custom", "read-number", "annotation_column"),
                          spike_in_cell_type = NULL,
                          spike_in_amount = NULL,
                          unique_cell_type = NULL,
                          custom_scenario_data = NULL,
                          custom_scaling_vector = NULL,
                          balance_unique_mirror_scenario=0.01,
                          sample_aggregation = "sum",
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

  # read-number information is necessary if sampling by sequencing depth is used
  if((!is.null(total_read_counts))&(!"read_number" %in% colnames(data@annotation))){
    stop("The provided dataset has no read number information stored. Please add it during generation of the dataset if you want to use total_read_counts as a parameter!")
  }

  ##### different cell-type scenarios #####

  # each existing cell-type will be appearing in equal amounts
  if(scenario == "uniform"){
    all_types <- unique(data@annotation[["cell_type"]])
    n_cell_types <- length(all_types)
    uniform_value <- 1/length(all_types)
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      m <- round(matrix(abs(rnorm(length(all_types), mean=1/length(all_types), sd=balance_unique_mirror_scenario)), ncol=n_cell_types), 3)
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
      n_cell_types <- length(unique(data@annotation[["cell_type"]]))
      # generate n_cell_type amount of random fractions from the uniform distribution, which will sum up to 1
      m <- round(matrix(runif(n_cell_types, 0, 1), ncol=n_cell_types),3)
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
    n_cell_types <- length(unique(data@annotation[["cell_type"]]))
    # generate 'nsamples' random samples
    simulation_vector_list <- lapply(rep(1:nsamples), function(x){
      # each cell-type will be represented as many times as it occurs in the used dataset
      mirror_values <- table(data@annotation$cell_type)/nrow(data@annotation)
      m <- unlist(lapply(mirror_values, function(y){
        return(abs(round(rnorm(1, mean=y, sd=balance_unique_mirror_scenario),3)))
      }))
      m <- matrix(m, ncol=n_cell_types)
      m <- sweep(m, 1, rowSums(m), FUN="/")
      simulation_vector <- as.vector(m[1,])
      names(simulation_vector) <- unique(data@annotation[["cell_type"]])
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
  idx <- 1
  bulk <- do.call(cbind, parallel::mclapply(simulation_vector_list, function(x){
    sample<-simulate_sample(data=data,
                            scaling_factor = scaling_factor,
                            scaling_vector = scaling_vector,
                            simulation_vector = x,
                            sample_aggregation = sample_aggregation,
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

  # remove non-unique features/genes from simulated dataset
  if(length(unique(rownames(bulk))) != dim(bulk)[1]){
    un <- unique(rownames(bulk))
    bulk <- bulk[un,]
  }

  # if sequencing depth is provided: rarefy samples to equal sequencing depth
  if(!is.null(total_read_counts)){
    if(sample_aggregation != "sum"){
      warning("You are not using 'sum' as sample aggregation method. Rarefying the generated samples to the provided sequencing depth could lead to unwanted results.")
    }
    #TODO change this 'if'.. should not be based on un-checked user input
    # if raw counts are used and the single cell matrix was changed by a scaling factor, we need rarefaction
    if("raw" %in% data@annotation$count_type && scaling_factor != "NONE"){
      message("The provided expression data are raw counts and a scaling factor was applied. The samples will be rarefied to match the wanted sequencing depth.")
      bulk <- phyloseq::otu_table(bulk, taxa_are_rows = T)
      bulk <- as.data.frame(phyloseq::rarefy_even_depth(physeq = bulk, sample.size = total_read_counts, rngseed = 1310, verbose = F, trimOTUs = F))
    }
  }

  # normalize count matrix to CPMs
  tryCatch({
    bulk_cpm <- cpm_normalize(bulk)
  }, error=function(e){
    warning(e$message)
    stop()
  })


  # build bioconductor expression set
  expr_set <- Biobase::ExpressionSet(assayData = bulk_cpm,
                                     phenoData = new("AnnotatedDataFrame", data=cell_fractions))

  cat("Finished simulation.\n")

  return(list(pseudo_bulk = bulk_cpm,
              pseudo_bulk_raw = bulk,
              cell_fractions = cell_fractions,
              scaling_vector = scaling_vector,
              expression_set = expr_set))

}


#' Calculate scaling factor for a dataset
#'
#' @param data dataset object
#' @param scaling_factor name of scaling factor; possible are: \code{census}, \code{spike-in}, \code{read-number}, \code{custom} or \code{NONE} for no scaling factor
#' @param custom_scaling_vector named vector with custom scaling values for cell-types. Cell-types that do not occur in this vector but are present in the dataset will be set to 1
#' @param ncores number of cores
#'
#' @return a named vector with a scaling value for each cell in the dataset
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

  }else if(scaling_factor == "read-number"){
    if(!"read_number" %in% colnames(data@annotation)){
      stop("The annotation in your dataset does not contain read-number information; you cannot apply the read-number scaling factor.")
    }
    tmp <- data@annotation[,c("cell_ID","read_number")]
    #tmp <- merge(tmp, simulated_annotation, by.x="cell_ID",by.y="values", all.y=T)
    scaling_vector <- tmp$read_number
    names(scaling_vector) <- tmp$cell_ID

  }else if (scaling_factor == "custom"){
    # needs vector with values for existing cell-types
    # cell-types that do not occur in this vector will have scaling-factor of 1
    if(is.null(custom_scaling_vector)){stop("For the custom scaling factor you need to provide a custom_scaling_vector!")}

    missing_cell_types <- as.vector(unique(data@annotation$cell_type)[which(!unique(data@annotation$cell_type) %in% names(custom_scaling_vector))])
    complete_vector <- rep(1, length(missing_cell_types))
    names(complete_vector) <- missing_cell_types
    complete_vector <- data.frame(value=append(complete_vector, custom_scaling_vector), check.names=F)
    df <- merge(complete_vector, data@annotation, by.x=0,by.y="cell_type", all.y=T)[,c("value","cell_ID")]
    scaling_vector <- df$value
    names(scaling_vector) <- df$cell_ID

  }else if(scaling_factor == "NONE"){
    scaling_vector <- rep(1, nrow(data@annotation))
    names(scaling_vector) <- data@annotation$cell_ID
  }else if(!is.null(scaling_factor)){
    print(paste0("Scaling by ", scaling_factor,"-column in annotation table; if no scaling is wished instead, use 'NONE'."))
    if(!scaling_factor %in% colnames(data@annotation)){stop(paste0("A column with the name ", scaling_factor," cannot be found in the annotation."))}
    if(!is.numeric(data@annotation[,c(scaling_factor)])){stop(paste0("The column with the name ", scaling_factor," is not numeric and cannot be used for scaling."))}

    tmp <- data@annotation[,c("cell_ID",scaling_factor)]
    colnames(tmp) <- c("a","b")
    scaling_vector <- tmp$b
    names(scaling_vector) <- tmp$a

  }else{
    warning("No valid scaling factor method provided. Scaling all cells by 1.")
    scaling_vector <- rep(1, nrow(data@annotation))
    names(scaling_vector) <- data@annotation$cell_ID
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



#' Sample cells based of sequencing depth
#'
#' This function allows to sample cells from a dataframe including the number of reads per cell until a pre-defined
#' sequencing depth is met or no cells are left with low enough read numbers to reach this depth.
#'
#' @param cells dataframe; needs to contain a "read_number" column
#' @param total_read_counts numeric; the total number of reads the sampled cells will sum up to (approximatly)
#'
#' @return list of cell IDs
#'
sample_cells_by_read_depth <- function(cells, depth){
  read_counts_current <- 0
  out <- list()
  missing_reads <- 0
  # sample cells and fill up the read_counts
  while(TRUE){
    sampled_cell <- dplyr::slice_sample(cells, n=1)
    cell_reads <- sampled_cell[["read_number"]]

    if((read_counts_current + cell_reads) > depth){
      sampled_cell <- NULL
      # the next sampled cell has too many reads -> check if other cells exist, which "fit in"
      difference_to_fill <-  depth - read_counts_current
      remaining_cells <- cells[cells[["read_number"]] <= difference_to_fill, ]  # these are all cells of type x, which have less reads than are we still have "available"
      # check if any cells are left, that can fill the remaining reads
      if(dim(remaining_cells)[1] == 0){
        #message(paste0("No more cells to sample from with read depth of ", difference_to_fill, " or lower."))
        missing_reads <- missing_reads + difference_to_fill
        break
      }
      # look for cell with number of reads lower then difference
      while (TRUE) {
        cell_checked <- dplyr::slice_sample(remaining_cells, n=1)
        cell_reads <- cell_checked[["read_number"]]
        if(cell_reads <= difference_to_fill){
          sampled_cell <- cell_checked
          cell_reads <- sampled_cell[["read_number"]]
          break
        }
      }
    }
    if(!is.null(sampled_cell)){
      read_counts_current <- read_counts_current + cell_reads
      out <- append(out, values = sampled_cell[["cell_ID"]])
    }
  }

  return(list(cells=out, missing_reads=missing_reads))
}
