library(monocle)
library(Seurat)

data_folder <- "/home/Data/Hao (CITEseq-PBMC)/hto_5p/"
mat <- Read10X(data_folder)

census_function <- function(expr_matrix, exp_capture_rate=0.25, expr_threshold=0.1){
  
  cells <- dim(expr_matrix)[2]
  idx <- 1
  # iterate over all cells
  total <- unlist(apply(expr_matrix, 2, function(x){
    # Find the most commonly occuring (log-transformed) TPM value in each cell above a threshold
    t_estimate <- 10^mean(dmode(log10(x[x > expr_threshold])))
    
    # only consider genes with TPM > 0.1; below this, no mRNA is believed to be present
    x <- x[x > 0.1]
    # calculate cumulative distribution function of gene expression values in cell
    P <- ecdf(x)
    # identify peak of distribution by looking at most common TPM value
    frac_x <- P(t_estimate) 
    
    # find all genes with single mRNA
    num_single_copy_genes <- sum(x <= t_estimate)

    #progress
    if(idx %% 1000 == 0){
      progress <- round(idx / cells, digits=3)
      print(paste(progress*100,"%..."))
    }
    idx <<- idx+1
    
    #final value (this is M_i)
    num_single_copy_genes / frac_x / exp_capture_rate
  }))
  
  return(total)
}


#use gaussian kernel to calculate the mode of transcript counts
#' @importFrom stats density
dmode <- function(x, breaks="Sturges") {
  if (length(x) < 2) return (0);
  den <- density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)] )
}



fractions <- census_function(mat)

