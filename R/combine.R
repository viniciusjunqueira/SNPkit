#' Safe cbind for SnpMatrix preserving dimnames
#'
#' This function performs a column-wise binding of multiple \code{SnpMatrix} objects,
#' explicitly preserving row names and column names, avoiding unexpected "object has no names" warnings.
#'
#' @param ... SnpMatrix objects to combine (must have identical row names).
#'
#' @return A single combined \code{SnpMatrix} with preserved row and column names.
#'
#' @examples
#' m1 <- methods::new("SnpMatrix",
#'                    matrix(as.raw(1:3), nrow = 3, ncol = 2,
#'                           dimnames = list(c("S1", "S2", "S3"),
#'                                           c("SNP1", "SNP2"))))
#' m2 <- methods::new("SnpMatrix",
#'                    matrix(as.raw(1:3), nrow = 3, ncol = 2,
#'                           dimnames = list(c("S1", "S2", "S3"),
#'                                           c("SNP3", "SNP4"))))
#' cbind_SnpMatrix(m1, m2)
#'
#' @export
cbind_SnpMatrix <- function(...) {
  mats <- list(...)

  # Check that all matrices have identical row names
  row_names_list <- lapply(mats, rownames)
  if (!all(sapply(row_names_list, function(x) all(x == row_names_list[[1]])))) {
    stop("All matrices must have identical row names to cbind safely.")
  }

  message("Performing safe cbind on SnpMatrix objects...")

  # Perform cbind using base
  res <- do.call(base::cbind, mats)

  # Reassign dimnames explicitly
  rownames(res) <- row_names_list[[1]]
  colnames(res) <- unlist(lapply(mats, colnames), use.names = FALSE)

  message("Row and column names preserved after cbind.")

  res
}

#' Safe rbind for SnpMatrix preserving dimnames
#'
#' This function performs a row-wise binding of multiple \code{SnpMatrix} objects,
#' explicitly preserving row names and column names, avoiding unexpected "object has no names" warnings.
#'
#' @param ... SnpMatrix objects to combine (must have identical column names).
#'
#' @return A single combined \code{SnpMatrix} with preserved row and column names.
#'
#' @examples
#' m1 <- methods::new("SnpMatrix",
#'                    matrix(as.raw(1:3), nrow = 2, ncol = 3,
#'                           dimnames = list(c("S1", "S2"),
#'                                           c("SNP1", "SNP2", "SNP3"))))
#' m2 <- methods::new("SnpMatrix",
#'                    matrix(as.raw(1:3), nrow = 2, ncol = 3,
#'                           dimnames = list(c("S3", "S4"),
#'                                           c("SNP1", "SNP2", "SNP3"))))
#' rbind_SnpMatrix(m1, m2)
#'
#' @export
rbind_SnpMatrix <- function(...) {
  mats <- list(...)

  # Check that all matrices have identical column names
  col_names_list <- lapply(mats, colnames)
  if (!all(sapply(col_names_list, function(x) all(x == col_names_list[[1]])))) {
    stop("All matrices must have identical column names to rbind safely.")
  }

  message("Performing safe rbind on SnpMatrix objects...")

  # Perform rbind using base
  res <- do.call(base::rbind, mats)

  # Reassign dimnames explicitly
  rownames(res) <- unlist(lapply(mats, rownames), use.names = FALSE)
  colnames(res) <- col_names_list[[1]]

  message("Row and column names preserved after rbind.")

  res
}

#' Combine multiple SNPDataLong objects
#'
#' This function merges a list of \code{SNPDataLong} objects, typically representing different SNP panels
#' or datasets, into a single unified \code{SNPDataLong} object. It ensures that all genotype matrices
#' have the same set of SNPs (filling missing SNPs with NA), and merges the marker map information while
#' removing duplicate SNP entries.
#'
#' @param lista A list of \code{SNPDataLong} objects to be combined.
#'
#' @return A single \code{SNPDataLong} object containing the combined genotype matrix, merged map,
#' and a concatenated path string.
#'
#' @examples
#' \donttest{
#' make_obj <- function(samples, snps) {
#'   m <- methods::new("SnpMatrix",
#'                     matrix(as.raw(1:3),
#'                            nrow = length(samples),
#'                            ncol = length(snps),
#'                            dimnames = list(samples, snps)))
#'   methods::new("SNPDataLong",
#'                geno = m,
#'                map  = data.frame(Name = snps,
#'                                  Chromosome = 1,
#'                                  Position = seq_along(snps)),
#'                path = tempfile(),
#'                xref_path = "chip1")
#' }
#' obj1 <- make_obj(c("S1", "S2"), c("SNP1", "SNP2"))
#' obj2 <- make_obj(c("S3", "S4"), c("SNP2", "SNP3"))
#' combined <- combineSNPData(list(obj1, obj2))
#' }
#'
#' @importFrom methods new as
#' @export
combineSNPData <- function(lista) {
  stopifnot(length(lista) > 0)

  message("Starting SNPDataLong combination...")

  snps_all <- Reduce(union, lapply(lista, function(x) colnames(x@geno)))
  message("Unified SNP panel with ", length(snps_all), " SNPs.")

  geno_list <- lapply(lista, function(x) {
    geno <- x@geno
    missing_snps <- setdiff(snps_all, colnames(geno))

    if (length(missing_snps) > 0) {
      message("Adding ", length(missing_snps), " missing SNPs filled with NA for one matrix...")
      na_block <- new("SnpMatrix", matrix(as.raw(0), nrow = nrow(geno), ncol = length(missing_snps)))
      colnames(na_block) <- missing_snps
      rownames(na_block) <- rownames(geno)
      geno <- cbind_SnpMatrix(geno, na_block)
    }

    geno <- geno[, snps_all, drop = FALSE]

    if (is.null(rownames(geno)) || any(rownames(geno) == "")) {
      warning("Some samples had missing rownames. Assigned default sample names.")
      rownames(geno) <- sprintf("Sample_%d", seq_len(nrow(geno)))
    }

    geno
  })

  geno_comb <- do.call(rbind_SnpMatrix, geno_list)

  map_all <- do.call(rbind, lapply(lista, function(x) x@map))
  map_all <- map_all[!duplicated(map_all$Name), , drop = FALSE]
  map_final <- map_all[match(snps_all, map_all$Name), , drop = FALSE]

  if (any(is.na(map_final$Name))) {
    missing_snps <- snps_all[is.na(match(snps_all, map_all$Name))]
    warning("Some SNPs missing in combined map: ", paste(missing_snps, collapse = ", "))
  }

  if (!inherits(geno_comb, "SnpMatrix")) {
    geno_comb <- as(geno_comb, "SnpMatrix")
  }

  message("Combination complete. Final matrix: ", nrow(geno_comb), " samples x ", ncol(geno_comb), " SNPs.")

  new("SNPDataLong",
      geno      = geno_comb,
      map       = map_final,
      path      = paste(sapply(lista, function(x) x@path), collapse = ";"),
      xref_path = paste(sapply(lista, function(x) x@xref_path), collapse = ";"))
}
