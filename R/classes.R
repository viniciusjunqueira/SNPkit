# Ensure that the SnpMatrix class is registered
if (requireNamespace("snpStats", quietly = TRUE)) {
  getClass("SnpMatrix", where = asNamespace("snpStats"))
  setClassUnion("SnpMatrixOrNULL", c("SnpMatrix", "NULL"))
} else {
  warning("Package 'snpStats' not found. Using fallback 'NULL' for SnpMatrixOrNULL.")
  setClassUnion("SnpMatrixOrNULL", "NULL")
}

# Define union for map slot before using it
setClassUnion("MapDataFrameOrList", c("data.frame", "list"))


#' SNPDataLong Class
#'
#' A class to store SNP genotype data in long format, including genotype matrix, marker map, and file paths.
#'
#' @slot geno A SnpMatrix containing genotype data. Must have individuals in rows and markers in columns.
#' @slot map A data.frame or list containing marker information. The number of markers in map should match the number of columns in geno.
#' @slot path A character string with the file path or identifier. Must be of length 1.
#' @slot xref_path A character vector with additional identifiers or paths. May have one entry per individual or a single value.
#'
#' Validity function for SNPDataLong (internal)
#' @noRd
#' @keywords internal
.validSNPDataLong <- function(object) {
  errors <- character()

  if (!inherits(object@geno, "SnpMatrix")) {
    errors <- c(errors, "Slot 'geno' must be an object of class SnpMatrix.")
  }

  if (!(is.data.frame(object@map) || is.list(object@map))) {
    errors <- c(errors, "Slot 'map' must be a data.frame or list.")
  }

  if (length(object@path) != 1 || !is.character(object@path)) {
    errors <- c(errors, "Slot 'path' must be a character string of length 1.")
  }

  if (!is.character(object@xref_path)) {
    errors <- c(errors, "Slot 'xref_path' must be a character vector.")
  }

  if (inherits(object@geno, "SnpMatrix") && (is.data.frame(object@map) || is.list(object@map))) {
    n_markers_geno <- ncol(object@geno)
    if (is.data.frame(object@map)) {
      n_markers_map <- nrow(object@map)
    } else if (is.list(object@map) && !is.null(object@map$markers)) {
      n_markers_map <- length(object@map$markers)
    } else {
      n_markers_map <- NA_integer_
    }

    if (!is.na(n_markers_map) && n_markers_geno != n_markers_map) {
      errors <- c(errors, sprintf(
        "Number of SNPs in 'geno' (%d) does not match the number of SNPs in 'map' (%d).",
        n_markers_geno, n_markers_map
      ))
    }
  }

  if (length(errors) == 0) TRUE else errors
}


setClass("SNPDataLong",
         slots = c(
           geno = "SnpMatrix",
           map = "MapDataFrameOrList",
           path = "character",
           xref_path = "character"
         ),
         validity = .validSNPDataLong)



#' SNPFileConfig Class
#'
#' A class for configuring SNP file import options.
#'
#' @slot path Path to the SNP file.
#' @slot fields A list specifying column mappings or field configurations.
#' @slot codes Character vector for genotype or allele codes.
#' @slot threshold Numeric value for filtering or quality control.
#' @slot sep Character specifying the field separator.
#' @slot skip Number of lines to skip at the top of the file.
#'
#' @export
setClass("SNPFileConfig",
         slots = c(
           path = "character",
           fields = "list",
           codes = "character",
           threshold = "numeric",
           sep = "character",
           skip = "numeric"
         ))


#' SNPImportList Class
#'
#' A class for managing a list of SNP file import configurations.
#'
#' @slot configs A list of SNPFileConfig objects.
#'
#' @export
setClass("SNPImportList",
         slots = c(
           configs = "list"
         ))


#' FImputeExport Class
#'
#' A class to handle export preparation for FImpute.
#'
#' @slot geno A SnpMatrix or NULL containing genotype data.
#' @slot map A data.frame containing marker information.
#' @slot path Output file path.
#' @slot name Project or file name.
#'
#' @export
setClass("FImputeExport",
         slots = c(
           geno = "SnpMatrixOrNULL",
           map = "data.frame",
           path = "character",
           name = "character"
         ))


#' FImputeRunner Class
#'
#' A class to manage FImpute execution and results.
#'
#' @slot export An FImputeExport object.
#' @slot par_file Path to parameter file.
#' @slot exec_path Path to FImpute executable.
#' @slot results A data.frame containing results or summary information.
#'
#' @export
setClass("FImputeRunner",
         slots = c(
           export = "FImputeExport",
           par_file = "character",
           exec_path = "character",
           results = "data.frame"
         ))
