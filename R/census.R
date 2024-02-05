#' Applies the Census count transformation on a count matrix
#'
#' needs a sparse matrix with cells in columns and genes in rows. You can find the detailed explaination here:
#' \url{http://cole-trapnell-lab.github.io/monocle-release/docs/#census}
#'
#' @param matrix sparse count matrix; cells in columns, genes in rows
#' @param exp_capture_rate expected capture rate; default=0.25
#' @param expr_threshold expression threshold; default=0
#' @param BPPARAM BiocParallel::bpparam() by default; if specific number of threads x want to be used, insert: BiocParallel::MulticoreParam(workers = x)
#' @param run_parallel boolean, decide if multi-threaded calculation will be run. FALSE by default
#'
#' @return a vector for each cell-type, with a scaling factor which can be used to transform the counts of the matrix
#' @export
#'
#' @examples
#'
#' tpm <- Matrix::Matrix(matrix(rpois(3e5, 5), ncol = 300), sparse = TRUE)
#' tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))
#' cen <- SimBu::census(tpm)
#'
census <- function(matrix, exp_capture_rate = 0.25, expr_threshold = 0, BPPARAM = BiocParallel::bpparam(), run_parallel = FALSE) {
  # switch multi-threading on/off
  if (!run_parallel) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = 1)
  }

  # order_cells <- colnames(matrix)
  ncuts <- dim(matrix)[2] / 1000

  # split matrix into cuts, each cut is a fractions of genes which are analyzed over all cells
  cuts <- split(seq_len(ncol(matrix)), cut(seq_len(ncol(matrix)), pretty(seq_len(ncol(matrix)), ncuts)))
  names(cuts) <- NULL

  out <- unlist(BiocParallel::bplapply(cuts, function(x) {
    x <- unlist(x, use.names = FALSE)
    chunk <- Matrix::Matrix(matrix[, x])

    cen <- census_monocle(chunk, exp_capture_rate = exp_capture_rate, expr_threshold = expr_threshold)
    return(cen)
  }, BPPARAM = BPPARAM))


  return(out)
}

calc_xi <- function(expr_matrix, expr_threshold) {
  cells <- dim(expr_matrix)[2]
  idx <- 1
  total <- unlist(apply(expr_matrix, 2, function(x) {
    # Find the most commonly occuring (log-transformed) TPM value in each cell above a threshold
    10^mean(dmode(log10(x[x > expr_threshold])))
  }))
}


#' Census calculation as implemented in monocle
#'
#' Implementation taken from Monocle2: https://github.com/cole-trapnell-lab/monocle-release/blob/master/R/normalization.R#L140
#'
#' @param expr_matrix TPM matrix
#' @param exp_capture_rate expected capture rate; default=0.25
#' @param expr_threshold expression threshold; default=0
#'
#' @return vector with estimated mRNA values per cell in expr_matrix
#'
#' @keywords internal
census_monocle <- function(expr_matrix, exp_capture_rate, expr_threshold) {
  # iterate over all cells
  total <- unlist(apply(expr_matrix, 2, function(x) {
    tryCatch(
      {
        # Find the most commonly occurring (log-transformed) TPM value in each cell above a threshold
        t_estimate <- 10^mean(dmode(log10(x[x > expr_threshold])))
        # only consider genes with TPM > 0.1; below this, no mRNA is believed to be present
        # x <- x[x > 0.1]
        x <- x[x > expr_threshold]
        # calculate cumulative distribution function of gene expression values in cell
        P <- stats::ecdf(x)
        # identify peak of distribution by looking at most common TPM value
        frac_x <- P(t_estimate)
        # find all genes with single mRNA
        num_single_copy_genes <- sum(x <= t_estimate)
        # final value (this is M_i)
        num_single_copy_genes / frac_x / exp_capture_rate
      },
      error = function(e) {
        message(e$message)
      }
    )
  }))

  return(total)
}

#' use gaussian kernel to calculate the mode of transcript counts
#'
#' @param x vector of numeric values
#'
#' @return most commonly occurring (log-transformed) TPM value
#'
#' @keywords internal
dmode <- function(x) {
  if (length(x) < 2) {
    return(0)
  }
  den <- stats::density(x, kernel = c("gaussian"))
  (den$x[den$y == max(den$y)])
}
