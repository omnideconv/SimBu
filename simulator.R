library(parallel)
library(data.table)
library(Matrix)
library(Seurat)
library(anndata)
library(tools)
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
  print("Finished sample ..")
  
  return(simulated_count_vector = simulated_count_vector)
  
}


# function to simulate a whole pseudo-bulk dataset
# data = dataset or database object
# scenario = pre-defined cell-type fractions
# nsamples = number of pseudo-bulk samples to be generated
simulate_bulk <- function(data, 
                          scenario=c("noisy"), 
                          scaling_factor=c("NONE","census","spike-in","custom"), 
                          nsamples=100, 
                          ncells=1000,
                          ncores = 1){
  
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
  
  # sample cells and generate pseudo-bulk profiles
  idx <- 1
  bulk <- do.call(cbind, mclapply(simulation_vector_list, function(x){
    sample<-simulate_sample(data=data, 
                            scaling_factor = scaling_factor,
                            simulation_vector = x,
                            total_cells = n_cells,
                            ncores=ncores)
    
    progress <- 100*(round(idx/length(simulation_vector_list), digits=3))
    print(paste(progress,"%"))
    idx <<- idx+1
    
    return(sample)
  }, mc.cores=ncores))
  
  colnames(bulk) <- sample_names

  
  return(list(pseudo_bulk = bulk,
              cell_fractions = simulation_vector_list))
  
}

#### tests #####

counts_maynard <- Matrix(as.matrix(fread("~/ma/data/Maynard/X_tpm.csv")), sparse = T)
genes <- fread("~/ma/data/Maynard/var.csv")
cells <- fread("~/ma/data/Maynard/obs.csv")
cellnames <- cells$Run
genenames <- genes$symbol
dimnames(counts_maynard)<-list(cellnames, genenames)
counts_maynard <- t(counts_maynard)

annotation_maynard <- fread("~/ma/data/Maynard/annotation_obs.csv")[,c("Run","cell_type")]
colnames(annotation_maynard) <- c("ID", "cell_type")

ds_m <- dataset(annotation_maynard, counts_maynard, name = "Maynard", count_type = "TPM")

annotation_travaglini <- fread("~/ma/data/Travaglini/obs_extended.csv")
counts_travaglini <- "~/ma/data/Travaglini/Travaglini_Krasnow_2020_Lung_SS2.h5ad"

ercc_cols <- grep("ERCC-",colnames(annotation_travaglini))
meta_red <- data.frame(annotation_travaglini)[ercc_cols]
meta_red$ercc_sum <- apply(meta_red,1,sum)
meta_red$cell_type <-annotation_travaglini$free_annotation
meta_red$total_reads <- annotation_travaglini$nReads
meta_red$genes <- annotation_travaglini$nGene
meta_red$ID <- annotation_travaglini$index
meta_red <- as.data.table(meta_red)

ds_t <- dataset_h5ad(meta_red, counts_travaglini, name = "Travaglini", spike_in_col = "ercc_sum")


db <- database(list(ds_m, ds_t))

sim<-simulate_bulk(db, scenario = "noisy", scaling_factor="census",ncores=4, ncells=5000) 

# sim$pseudo_bulk %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(-c(gene), names_to = "samples", values_to = "counts") %>%
#   ggplot(aes(x=samples,y=gene,fill=counts))+
#   geom_raster()+
#   scale_fill_viridis_c()
