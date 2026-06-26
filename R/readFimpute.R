#' @useDynLib SNPkit, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importClassesFrom snpStats SnpMatrix
NULL

#' Read imputed genotypes from FImpute output and return SNPDataLong object
#'
#' Reads imputed genotypes and SNP information from FImpute output,
#' builds a SnpMatrix and a corresponding map, and returns an SNPDataLong object.
#'
#' @param file Character. Path to the FImpute output directory (usually "output_fimpute").
#' @param method Character. "R" (default) for vectorized R implementation, or "Rcpp" for compiled C++ implementation.
#'
#' @return An object of class \code{SNPDataLong} with three slots:
#'   \code{geno} (a \code{SnpMatrix} with individuals as rows and SNPs as
#'   columns), \code{map} (a \code{data.frame} with columns \code{Name},
#'   \code{Chromosome}, and \code{Position}), and \code{path} (the input
#'   directory).
#' @examples
#' \dontrun{
#' # Requires a directory containing FImpute output files
#' # (genotypes_imp.txt and snp_info.txt).
#' snp_long <- read.fimpute("output_fimpute", method = "R")
#' }
#'
#' @export
read.fimpute <- function(file, method = c("R", "Rcpp")) {
  method <- match.arg(method)

  genotype_file <- file.path(file, "genotypes_imp.txt")
  snp_info_file <- file.path(file, "snp_info.txt")

  # Read SNP map using fread
  map <- data.table::fread(snp_info_file, data.table = FALSE)

  # Check required columns
  required_cols <- c("SNPID", "Chr", "BPPos")
  if (!all(required_cols %in% colnames(map))) {
    stop("snp_info.txt must contain columns: SNPID, Chr, BPPos")
  }

  snp_names <- map$SNPID

  # Determine dimensions
  temp_header <- data.table::fread(genotype_file, nrows = 1, data.table = FALSE)
  nc <- nchar(temp_header[1, 3])
  nrows <- nrow(data.table::fread(genotype_file, select = 1, data.table = FALSE))  # cross-platform fix

  message("Number of individuals: ", nrows)
  message("Number of SNPs: ", nc)

  if (method == "R") {
    temp <- data.table::fread(genotype_file, data.table = FALSE)
    char_vec <- temp[, 3]
    split_list <- stringi::stri_split_boundaries(char_vec, type = "character")
    mat <- do.call(rbind, lapply(split_list, as.integer))

    mat[mat == 2] <- 3L
    mat[mat == 1] <- 2L
    mat[mat == 0] <- 1L

    rownames(mat) <- temp[, 1]
  } else if (method == "Rcpp") {
    mat <- readFImputeCpp(genotype_file, nrows, nc)
    ids <- data.table::fread(genotype_file, select = 1, data.table = FALSE)[, 1]
    rownames(mat) <- ids
  }

  colnames(mat) <- snp_names

  snp_matrix <- new("SnpMatrix", mat)

  map_out <- data.frame(
    Name = map$SNPID,
    Chromosome = map$Chr,
    Position = map$BPPos,
    stringsAsFactors = FALSE
  )

  snp_long <- new("SNPDataLong", geno = snp_matrix, map = map_out, path = file, xref_path = "")
  return(snp_long)
}

#' Import imputed FImpute results from disk
#'
#' Reads existing imputed results from a given path and returns an object of class SNPDataLong.
#'
#' @param path Character. Path to the folder containing 'output_fimpute' (e.g., "fimpute_run_nelore").
#' @param method Character. "R" (default) or "Rcpp". Passed to read.fimpute().
#'
#' @return An object of class SNPDataLong containing the imputed genotypes and SNP map.
#' @export
importFImputeResults <- function(path, method = "R") {
  if (!is.character(path) || length(path) != 1) {
    stop("'path' must be a single character string.")
  }

  output_dir <- file.path(path, "output_fimpute")

  if (!dir.exists(output_dir)) {
    stop("Output directory 'output_fimpute' does not exist at: ", output_dir)
  }

  if (!exists("read.fimpute", mode = "function")) {
    stop("The function 'read.fimpute()' must be defined and available in the current environment.")
  }

  message("Reading FImpute results from: ", output_dir)
  res <- read.fimpute(file = output_dir, method = method)

  message("Results successfully loaded as SNPDataLong object.")
  return(res)
}
