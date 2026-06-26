#' Quality Control for SNPDataLong with optional criteria
#'
#' Applies flexible quality control filters on an object of class \code{SNPDataLong}.
#' Supports call rate filtering, minor allele frequency (MAF), Hardy-Weinberg equilibrium (HWE),
#' removal of monomorphic SNPs, exclusion of specific chromosomes, optionally removing SNPs without positions,
#' and optionally removing SNPs at the same genomic position (keeping the one with highest MAF).
#'
#' @param x An object of class SNPDataLong.
#' @param missing_ind Maximum allowed proportion of missing data per individual (currently not implemented).
#' @param missing_snp Maximum allowed proportion of missing data per SNP (currently not implemented).
#' @param min_snp_cr Minimum acceptable call rate for SNPs (e.g., 0.95). SNPs below this threshold are removed.
#' @param min_maf Minimum minor allele frequency allowed for SNPs (e.g., 0.05). SNPs with lower MAF are removed.
#' @param hwe p-value threshold for Hardy-Weinberg equilibrium test (e.g., 1e-6). SNPs violating this are removed.
#' @param snp_position Logical. If TRUE, removes SNPs mapped to the same position, retaining only the one with highest MAF.
#' @param no_position Logical. If TRUE, removes SNPs without defined genomic positions.
#' @param snp_mono Logical. If TRUE, removes monomorphic SNPs (with no variation).
#' @param remove_chr Character vector of chromosomes to exclude (e.g., c("X", "Y")).
#' @param action One of "report" (returns a list of removed SNPs), "filter" (returns filtered SNPDataLong), or "both" (returns both).
#' @param ... Additional optional arguments.
#'
#' @return Depending on the action argument:
#' - "report": list of SNPs removed by each filter and SNPs retained.
#' - "filter": filtered SNPDataLong object.
#' - "both": list containing the filtered object and detailed report.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' raw_mat <- matrix(as.raw(sample(1:3, 100, TRUE)), nrow = 10, ncol = 10)
#' colnames(raw_mat) <- paste0("snp", 1:10)
#' rownames(raw_mat) <- paste0("ind", 1:10)
#' geno <- methods::new("SnpMatrix", raw_mat)
#' map <- data.frame(Name = colnames(geno), Chromosome = 1, Position = 1:10)
#' x <- methods::new("SNPDataLong",
#'                   geno = geno,
#'                   map  = map,
#'                   path = tempfile(),
#'                   xref_path = "chip1")
#'
#' qcSNPs(x,
#'        min_snp_cr = 0.8,
#'        min_maf = 0.05,
#'        snp_mono = TRUE,
#'        no_position = TRUE,
#'        snp_position = TRUE,
#'        action = "filter")
#' }
#'
#' @export
setGeneric("qcSNPs", function(x, ...) standardGeneric("qcSNPs"))

#' @rdname qcSNPs
#' @export
setMethod("qcSNPs", "SNPDataLong", function(x,
                                            missing_ind = NULL,
                                            missing_snp = NULL,
                                            min_snp_cr = NULL,
                                            min_maf = NULL,
                                            hwe = NULL,
                                            snp_position = NULL,
                                            no_position = NULL,
                                            snp_mono = FALSE,
                                            remove_chr = NULL,
                                            action = c("report", "filter", "both")) {
  action <- match.arg(action)

  geno <- x@geno
  map <- x@map
  keep_snps <- colnames(geno)

  qc_header("Quality Control on SNPs")

  message("Initial number of SNPs: ", length(keep_snps))
  message("Applying quality control filters...")

  snpsum <- snpStats::col.summary(geno)

  # Initialize vectors to store SNPs to be removed
  low_callrate_snps <- low_maf <- dev.hwe <- mono <- discard_chr <- snpstoremove <- no_pos <- character()

  # Filter by call rate
  if (!is.null(min_snp_cr)) {
    low_callrate_snps <- check.call.rate(snpsum, min.call.rate = min_snp_cr)
    keep_snps <- setdiff(keep_snps, low_callrate_snps)
    message(sprintf("  - Call rate filter: %d SNP(s) removed; %d retained.",
                    length(low_callrate_snps), length(keep_snps)))
  }

  # Filter by MAF
  if (!is.null(min_maf)) {
    low_maf <- check.snp.maf(snpsum, min.maf = min_maf)
    keep_snps <- setdiff(keep_snps, low_maf)
    message(sprintf("  - MAF filter: %d SNP(s) removed; %d retained.",
                    length(low_maf), length(keep_snps)))
  }

  # Filter by HWE
  if (!is.null(hwe)) {
    dev.hwe <- check.snp.hwe.chi2(snpsum, hwe)
    keep_snps <- setdiff(keep_snps, dev.hwe)
    message(sprintf("  - HWE filter: %d SNP(s) removed; %d retained.",
                    length(dev.hwe), length(keep_snps)))
  }

  # Filter SNPs with no position
  if (!is.null(no_position) && no_position) {
    no_pos <- check.snp.no.position(map)
    keep_snps <- setdiff(keep_snps, no_pos)
    message(sprintf("  - No-position filter: %d SNP(s) removed; %d retained.",
                    length(no_pos), length(keep_snps)))
  }

  # Filter SNPs at the same genomic position
  if (!is.null(snp_position) && snp_position) {
    # Force removal of SNPs with missing position before position filter
    if (is.null(no_position) || !no_position) {
      no_pos_manual <- map$Name[is.na(map$Position)]
      if (length(no_pos_manual) > 0) {
        warning("SNPs without position removed automatically before snp_position filter.")
        keep_snps <- setdiff(keep_snps, no_pos_manual)
        map <- map[map$Name %in% keep_snps, , drop = FALSE]
        snpsum <- snpsum[keep_snps, , drop = FALSE]
      }
    }

    snp_same <- check.snp.same.position(map)
    message("  - Positions with overlapping SNPs: ", length(snp_same))
    n <- length(snp_same)
    if (n > 0) {
      for (i in seq_len(n)) {
        snpsum1 <- snpsum[snp_same[[i]], , drop = FALSE]
        snp.high.maf <- rownames(snpsum1[snpsum1[, "MAF"] == max(snpsum1[, "MAF"], na.rm = TRUE), , drop = FALSE])[1]
        snpstoremove <- union(snpstoremove, setdiff(rownames(snpsum1), snp.high.maf))
      }
      keep_snps <- setdiff(keep_snps, snpstoremove)
      message(sprintf("  - Same-position filter: %d SNP(s) removed; %d retained.",
                      length(snpstoremove), length(keep_snps)))
    }
  }

  # Filter monomorphic SNPs
  if (snp_mono) {
    mono <- check.snp.monomorf(snpsum)
    keep_snps <- setdiff(keep_snps, mono)
    message(sprintf("  - Monomorphic filter: %d SNP(s) removed; %d retained.",
                    length(mono), length(keep_snps)))
  }

  # Filter by chromosomes
  if (!is.null(remove_chr)) {
    discard_chr <- check.snp.chromo(map, remove_chr)
    keep_snps <- setdiff(keep_snps, discard_chr)
    message(sprintf("  - Chromosome filter: %d SNP(s) removed; %d retained.",
                    length(discard_chr), length(keep_snps)))
  }

  # Report-only output
  if (action == "report") {
    return(list(
      removed_by_callrate = low_callrate_snps,
      removed_by_maf = low_maf,
      removed_by_hwe = dev.hwe,
      removed_no_position = no_pos,
      removed_same_position = snpstoremove,
      removed_monomorphic = mono,
      removed_by_chr = discard_chr,
      kept_snps = keep_snps
    ))
  }

  # Apply final filtering
  filtered_geno <- geno[, keep_snps, drop = FALSE]
  filtered_map  <- map[map$Name %in% keep_snps, , drop = FALSE]
  filtered_obj  <- new("SNPDataLong",
                       geno = filtered_geno,
                       map = filtered_map,
                       path = x@path,
                       xref_path = x@xref_path)

  if (action == "filter") {
    return(filtered_obj)
  }

  if (action == "both") {
    return(list(
      filtered = filtered_obj,
      report = list(
        removed_by_callrate = low_callrate_snps,
        removed_by_maf = low_maf,
        removed_by_hwe = dev.hwe,
        removed_no_position = no_pos,
        removed_same_position = snpstoremove,
        removed_monomorphic = mono,
        removed_by_chr = discard_chr,
        kept_snps = keep_snps
      )
    ))
  }
})


#' Check SNPs mapped to the same position
#'
#' Identifies groups of SNPs that are mapped to the exact same genomic position on each chromosome.
#' Returns a list where each element corresponds to one group of overlapping SNPs.
#'
#' @param snpmap Data frame containing at least columns "Name", "Chromosome", and "Position".
#'
#' @return A list of character vectors, each with names of SNPs found at the same position.
#'
#' @export
check.snp.same.position <- function(snpmap) {
  chromo <- unique(snpmap[, "Chromosome"])
  n <- length(chromo)
  snps <- list()
  k <- 1
  for (i in seq_len(n)) {
    message("Analyzing chromosome ", chromo[i])
    snpmap.chr <- snpmap[snpmap[, "Chromosome"] == chromo[i], ]

    snpmap.chr <- snpmap.chr[!is.na(snpmap.chr[, "Position"]), , drop = FALSE]

    m <- nrow(snpmap.chr)
    if (m < 2) {
      next
    }

    sorted.snpmap.chr <- snpmap.chr[order(snpmap.chr[, "Position"]), ]

    for (j in 1:(m - 1)) {
      j1 <- j + 1
      pos_j <- sorted.snpmap.chr[j, "Position"]
      pos_j1 <- sorted.snpmap.chr[j1, "Position"]

      # Check robusto
      if (!is.na(pos_j) && !is.na(pos_j1) && pos_j == pos_j1) {
        message("SNPs in same position: ", sorted.snpmap.chr[j, "Name"], " - ", sorted.snpmap.chr[j1, "Name"])
        if (length(snps) < k) {
          snps[[k]] <- c(as.character(sorted.snpmap.chr[j, "Name"]), as.character(sorted.snpmap.chr[j1, "Name"]))
        } else {
          snps[[k]] <- c(snps[[k]], as.character(sorted.snpmap.chr[j1, "Name"]))
        }
      } else {
        if (length(snps) == k) {
          k <- k + 1
        }
      }
    }
  }
  return(snps)
}
