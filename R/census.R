# needs sparse matrix with cells in cols, genes in rows
#' applies the census count transformation on a count matrix; needs a sparse matrix with cells in columns and genes in rows
#'
#' @param matrix sparse count matrix; cells in columns, genes in rows
#' @param exp_capture_rate expected capture rate; default=0.25
#' @param expr_threshold expression threshold; default=0
#' @param ncores number of cores
#' @param method implementation of census; possible are: \code{monocle},\code{paper},\code{t_estimate}
#'
#' @return a vector for each cell-type, with a scaling factor which can be used to transform the counts of the matrix
#' @export
#'
#' @examples
census <- function(matrix, exp_capture_rate=0.25, expr_threshold=0, ncores=1, method=c("monocle","paper","t_estimate")){
  ncuts <- dim(matrix)[2]/1000

  #split matrix into cuts, each cut is a fractions of genes which are analyzed over all cells
  cuts<-split(seq_len(ncol(matrix)), cut(seq_len(ncol(matrix)), pretty(seq_len(ncol(matrix)), ncuts)))
  names(cuts) <- NULL

  idx <- 1
  out <- unlist(mclapply(cuts, function(x){
    x <- unlist(x, use.names = F)
    chunk <- matrix[,x]
    if(method == "monocle"){
      cen <- census_monocle(chunk, exp_capture_rate=exp_capture_rate, expr_threshold=expr_threshold)
    }else if (method == "paper"){
      cen <- census_paper(chunk, exp_capture_rate=exp_capture_rate, expr_threshold=expr_threshold)
    }else if (method == "t_estimate"){
      cen <- calc_xi(chunk, expr_threshold=expr_threshold)
    }

    progress <- 100*(round(idx/ncuts, digits=3))
    print(paste(progress,"%"))
    idx <<- idx+1

    return(cen)
  },mc.cores = ncores))

  return(out)
}

calc_xi <- function(expr_matrix, expr_threshold){
  cells <- dim(expr_matrix)[2]
  idx <- 1
  total <- unlist(apply(expr_matrix, 2, function(x){
    # Find the most commonly occuring (log-transformed) TPM value in each cell above a threshold
    10^mean(dmode(log10(x[x > expr_threshold])))
  }))

}


census_monocle <- function(expr_matrix, exp_capture_rate, expr_threshold){

  cells <- dim(expr_matrix)[1]
  i <- 1
  # iterate over all cells
  total <- unlist(apply(expr_matrix, 2, function(x){
    tryCatch({
      # Find the most commonly occurring (log-transformed) TPM value in each cell above a threshold
      t_estimate <- 10^mean(dmode(log10(x[x > expr_threshold])))
      # only consider genes with TPM > 0.1; below this, no mRNA is believed to be present
      #x <- x[x > 0.1]
      x <- x[x > expr_threshold]
      # calculate cumulative distribution function of gene expression values in cell
      P <- ecdf(x)
      # identify peak of distribution by looking at most common TPM value
      frac_x <- P(t_estimate)
      # find all genes with single mRNA
      num_single_copy_genes <- sum(x <= t_estimate)
      # counter
      i <<- i+1
      #final value (this is M_i)
      num_single_copy_genes / frac_x / exp_capture_rate
    }, error=function(e){
      print(x)
      print(e)
      break
    })
  }))

  return(total)
}

census_paper <- function(expr_matrix, exp_capture_rate, expr_threshold){
  cells <- dim(expr_matrix)[2]
  idx <- 1
  # iterate over all cells
  total <- unlist(apply(expr_matrix, 2, function(x){
      tryCatch({
        # Find the most commonly occuring (log-transformed) TPM value in each cell above a threshold
        x_star <- dmode(log10(x[x>expr_threshold]))
        # only consider genes with TPM > 0.1; below this, no mRNA is believed to be present
        x <- x[x>expr_threshold]
        # calculate cumulative distribution function of gene expression values in cell
        P <- ecdf(x)

        F_x_star <- P(x_star)
        F_x_epsilon <- P(expr_threshold)
        # find all genes with single mRNA
        num_single_copy_genes <- sum(x <= x_star)
        # counter
        idx <<- idx +1
        #final value (this is M_i)
        (1/exp_capture_rate) * (num_single_copy_genes / (F_x_star - F_x_epsilon))
      }, error=function(e){
        print(idx)
        print(e)
        break
      })
  }))
}


#use gaussian kernel to calculate the mode of transcript counts
#' @importFrom stats density
dmode <- function(x, breaks="Sturges") {
  if (length(x) < 2) return (0);
  den <- density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)] )
}



#census(t(m), expr_threshold = 0.1,ncores=1, method="monocle")
#maynard_census_1 <- census(t(mat), expr_threshold = 0.1,ncores=1, method="monocle")
