#' Summary for SNPDataLong objects
#'
#' Provides a detailed summary of an \code{SNPDataLong} object, including sample
#' and SNP counts, proportion of missing data, and SNP distribution by chromosome
#' if mapping information is available.
#'
#' @param object An object of class \code{SNPDataLong}.
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class \code{summary.SNPDataLong}, which is a list with
#'   the following elements:
#'   \describe{
#'     \item{n_individuals}{Integer. Number of individuals (rows of \code{geno}).}
#'     \item{n_snps}{Integer. Number of SNPs (columns of \code{geno}).}
#'     \item{n_missing}{Integer. Total number of missing genotype calls.}
#'     \item{prop_missing}{Numeric. Proportion of missing genotype calls.}
#'     \item{by_chromosome}{Either a table of SNP counts per chromosome (when
#'        the map provides \code{Name} and \code{Chromosome}) or \code{NULL}.}
#'     \item{missing_by_chromosome}{Either a table of SNPs with at least one
#'        missing call per chromosome, or \code{NULL}.}
#'   }
#'   The object also has a dedicated \code{print} method that displays the
#'   summary on the console.
#'
#' @export
setMethod("summary", "SNPDataLong", function(object, ...) {
  res <- list(
    n_individuals          = 0L,
    n_snps                 = 0L,
    n_missing              = NA_integer_,
    prop_missing           = NA_real_,
    by_chromosome          = NULL,
    missing_by_chromosome  = NULL,
    valid                  = FALSE,
    note                   = NULL
  )

  if (!inherits(object@geno, "SnpMatrix")) {
    res$note <- "Slot 'geno' is not a valid SnpMatrix."
    class(res) <- "summary.SNPDataLong"
    return(res)
  }

  if (!is.data.frame(object@map)) {
    res$note <- "Slot 'map' is not a data.frame."
    class(res) <- "summary.SNPDataLong"
    return(res)
  }

  res$valid <- TRUE
  res$n_individuals <- nrow(object@geno)
  res$n_snps <- ncol(object@geno)

  if (res$n_individuals == 0 || res$n_snps == 0) {
    res$note <- "Empty object: no individuals or SNPs."
    class(res) <- "summary.SNPDataLong"
    return(res)
  }

  n_total <- res$n_individuals * res$n_snps
  res$n_missing <- sum(is.na(object@geno))
  res$prop_missing <- res$n_missing / n_total

  if ("Name" %in% colnames(object@map) && "Chromosome" %in% colnames(object@map)) {
    idx <- match(colnames(object@geno), object@map$Name)
    chr_info <- object@map$Chromosome[idx]

    if (any(is.na(chr_info))) {
      res$note <- "Some SNPs were not found in the map and were ignored in chromosome counts."
    }

    chr_info_clean <- chr_info[!is.na(chr_info)]
    if (length(chr_info_clean) > 0) {
      res$by_chromosome <- table(chr_info_clean)

      snp_na_count <- colSums(is.na(object@geno))
      chr_na <- chr_info
      names(chr_na) <- colnames(object@geno)[idx]
      res$missing_by_chromosome <- tapply(snp_na_count > 0, chr_na, sum)
    }
  } else {
    res$note <- "Map does not contain expected columns ('Name', 'Chromosome'); chromosome summary omitted."
  }

  class(res) <- "summary.SNPDataLong"
  res
})

#' Print method for SNPDataLong summary
#'
#' Displays the contents of a \code{summary.SNPDataLong} object on the console.
#'
#' @param x An object of class \code{summary.SNPDataLong}.
#' @param ... Further arguments (currently unused).
#'
#' @return The input \code{x}, returned invisibly.
#' @export
print.summary.SNPDataLong <- function(x, ...) {
  cat("Summary of SNPDataLong object\n")
  cat("-----------------------------\n")

  if (!isTRUE(x$valid)) {
    cat(x$note, "\n", sep = "")
    return(invisible(x))
  }

  cat("Individuals :", x$n_individuals, "\n")
  cat("SNPs        :", x$n_snps, "\n\n")

  cat("Missing data (NA):\n")
  cat(" - Total     :", x$n_missing, "of",
      x$n_individuals * x$n_snps, "\n")
  cat(" - Proportion:", round(100 * x$prop_missing, 2), "%\n\n")

  if (!is.null(x$by_chromosome)) {
    if (!is.null(x$note)) {
      cat("Note:", x$note, "\n")
    }
    cat("Distribution of SNPs by chromosome:\n")
    print(x$by_chromosome)
    cat("\n")

    if (!is.null(x$missing_by_chromosome)) {
      cat("SNPs with missing data by chromosome:\n")
      print(x$missing_by_chromosome)
    }
  } else if (!is.null(x$note)) {
    cat(x$note, "\n")
  }

  invisible(x)
}
