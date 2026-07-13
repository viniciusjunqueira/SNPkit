# Test helpers for combineSNPData equivalence -------------------------------
#
# These build small SNPDataLong objects and provide a faithful copy of the
# PRE-refactor combineSNPData implementation (`combineSNPData_legacy`) so the
# new implementation can be checked to return byte-identical results.

## ---- Builders -------------------------------------------------------------

# A SnpMatrix with the given sample (row) and SNP (column) names, filled with
# non-missing codes 1-3 (0x01 AA, 0x02 AB, 0x03 BB). Deterministic when `seed`
# is supplied.
make_snpmatrix <- function(samples, snps, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  vals <- as.raw(sample(1:3, length(samples) * length(snps), replace = TRUE))
  m <- matrix(vals, nrow = length(samples), ncol = length(snps),
              dimnames = list(samples, snps))
  methods::new("SnpMatrix", m)
}

# A minimal valid SNPDataLong wrapping such a matrix plus a matching map.
make_obj <- function(samples, snps, seed = NULL, chr = 1L,
                     path = "p", xref = "chip") {
  geno <- make_snpmatrix(samples, snps, seed)
  map <- data.frame(Name = snps, Chromosome = chr,
                    Position = seq_along(snps), stringsAsFactors = FALSE)
  methods::new("SNPDataLong", geno = geno, map = map,
               path = path, xref_path = xref)
}

## ---- Reference (pre-refactor) implementation ------------------------------
# Verbatim behavior of combineSNPData before the rbindSnpFlexible refactor,
# reusing the same internal primitives so it is a true golden reference. Only
# message() calls were dropped (they do not affect the returned object).
combineSNPData_legacy <- function(lista) {
  stopifnot(length(lista) > 0)

  snps_all <- Reduce(union, lapply(lista, function(x) colnames(x@geno)))

  geno_list <- lapply(lista, function(x) {
    geno <- x@geno
    missing_snps <- setdiff(snps_all, colnames(geno))

    if (length(missing_snps) > 0) {
      na_block <- methods::new("SnpMatrix", matrix(as.raw(0),
        nrow = nrow(geno),
        ncol = length(missing_snps),
        dimnames = list(rownames(geno), missing_snps)))
      geno <- SNPkit:::cbind_SnpMatrix(geno, na_block)
    }

    geno <- geno[, snps_all, drop = FALSE]

    if (is.null(rownames(geno)) || any(rownames(geno) == "")) {
      rownames(geno) <- sprintf("Sample_%d", seq_len(nrow(geno)))
    }

    geno
  })

  geno_comb <- do.call(SNPkit:::rbind_SnpMatrix, geno_list)

  map_all <- do.call(rbind, lapply(lista, function(x) x@map))
  map_all <- map_all[!duplicated(map_all$Name), , drop = FALSE]
  map_final <- map_all[match(snps_all, map_all$Name), , drop = FALSE]

  if (!inherits(geno_comb, "SnpMatrix")) {
    geno_comb <- methods::as(geno_comb, "SnpMatrix")
  }

  methods::new("SNPDataLong",
    geno      = geno_comb,
    map       = map_final,
    path      = paste(sapply(lista, function(x) x@path), collapse = ";"),
    xref_path = paste(sapply(lista, function(x) x@xref_path), collapse = ";"))
}

## ---- Quiet wrappers (silence message() noise, keep warnings) ---------------
combine_new    <- function(lista) suppressMessages(combineSNPData(lista))
combine_legacy <- function(lista) suppressMessages(combineSNPData_legacy(lista))

## ---- Comparators ----------------------------------------------------------

# Underlying raw bytes of a SnpMatrix, stripped of attributes.
snpmatrix_bytes <- function(x) {
  d <- x@.Data
  attributes(d) <- NULL
  d
}

expect_snpmatrix_equal <- function(a, b) {
  expect_true(inherits(a, "SnpMatrix"))
  expect_true(inherits(b, "SnpMatrix"))
  expect_identical(dim(a), dim(b))
  expect_identical(dimnames(a), dimnames(b))
  expect_identical(snpmatrix_bytes(a), snpmatrix_bytes(b))
}

expect_snpdata_equal <- function(a, b) {
  expect_snpmatrix_equal(a@geno, b@geno)
  expect_identical(a@map, b@map)
  expect_identical(a@path, b@path)
  expect_identical(a@xref_path, b@xref_path)
}
