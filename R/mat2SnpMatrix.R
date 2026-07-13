#' Convert a genotype matrix or data.frame to snpStats::SnpMatrix
#'
#' This function converts a genotype matrix coded as 0/1/2/NA or AA/AB/BB to a
#' \code{snpStats::SnpMatrix} object. It includes checks for coding validity,
#' missing values, and duplicate sample or SNP IDs, and preserves row and column
#' names from the input.
#'
#' @details
#' The function accepts both \code{matrix} and \code{data.frame} inputs. For
#' \code{data.frame} objects, all columns are coerced to a common type using
#' \code{as.matrix()}, which preserves \code{rownames} and \code{colnames}.
#'
#' The returned \code{SnpMatrix} object stores each genotype as a single byte,
#' which is memory-efficient compared to integer storage. However, large datasets
#' still require substantial RAM. For very large genotype sets, consider using
#' on-disk formats such as \pkg{SNPRelate} (GDS) or \pkg{bigsnpr}.
#'
#' @param geno A samples x SNPs matrix or data.frame with genotypes coded as
#'        0, 1, 2, or NA. Can be numeric/integer or character. \code{rownames} =
#'        sample IDs, \code{colnames} = SNP IDs.
#' @param coding One of \code{"012"} or \code{"AAABBB"}. For character inputs only.
#'        \code{"012"} expects \code{"0"}, \code{"1"}, \code{"2"}, and \code{missing_codes}.
#'        \code{"AAABBB"} expects \code{"AA"}, \code{"AB"}, \code{"BB"}, and \code{missing_codes}.
#' @param missing_codes Character values to treat as missing (only used when
#'        \code{geno} is character), e.g., \code{c("NA","-9",".")}.
#' @param check_ids If \code{TRUE}, verifies that row and column names are unique
#'        (recommended).
#'
#' @return A \code{snpStats::SnpMatrix} with the same \code{dimnames} as \code{geno}.
#'
#' @examples
#' # Numeric 0/1/2 with NAs
#' set.seed(1)
#' geno <- matrix(sample(c(0L,1L,2L,NA), 20, replace=TRUE), nrow=5)
#' rownames(geno) <- paste0("ind", 1:5)
#' colnames(geno) <- paste0("snp", 1:4)
#' SM <- as_snpmatrix(geno)
#'
#' # Character AA/AB/BB
#' geno_c <- matrix(sample(c("AA","AB","BB","."), 20, replace=TRUE,
#'                         prob=c(.35,.3,.3,.05)), nrow=5)
#' rownames(geno_c) <- rownames(geno)
#' colnames(geno_c) <- colnames(geno)
#' SMc <- as_snpmatrix(geno_c, coding="AAABBB", missing_codes=".")
#'
#' @importFrom methods new
#' @importClassesFrom snpStats SnpMatrix
#' @export
as_snpmatrix <- function(geno,
                         coding = c("012","AAABBB"),
                         missing_codes = c("NA","-9",".",""),
                         check_ids = TRUE) {
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("Package 'snpStats' is required. Please install it.")
  }
  coding <- match.arg(coding)

  # Coerce to matrix early
  if (!is.matrix(geno)) geno <- as.matrix(geno)

  # Validate dimnames
  if (is.null(rownames(geno)) || is.null(colnames(geno))) {
    stop("geno must have rownames (samples) and colnames (SNP IDs).")
  }
  if (check_ids) {
    if (anyDuplicated(rownames(geno))) stop("Duplicate sample IDs in rownames(geno).")
    if (anyDuplicated(colnames(geno))) stop("Duplicate SNP IDs in colnames(geno).")
  }

  # Normalize to integer 0/1/2 with NA
  if (is.character(geno)) {
    if (coding == "012") {
      geno[geno %in% missing_codes] <- NA_character_
      allowed <- c("0","1","2", NA_character_)
      bad <- !(geno %in% allowed)
      if (any(bad, na.rm = TRUE)) {
        idx <- which(bad, arr.ind = TRUE)[1, , drop=TRUE]
        stop(sprintf("Invalid genotype '%s' at [%d,%d]. Expect '0','1','2' or missing.",
                     geno[idx[1], idx[2]], idx[1], idx[2]))
      }
      geno <- suppressWarnings(matrix(as.integer(geno),
                                      nrow=nrow(geno),
                                      dimnames=dimnames(geno)))
    } else if (coding == "AAABBB") {
      geno[geno %in% missing_codes] <- NA_character_
      map <- c("AA"=0L, "AB"=1L, "BA"=1L, "BB"=2L) # treat BA as AB if present
      known <- names(map)
      bad <- !(is.na(geno) | geno %in% known)
      if (any(bad, na.rm = TRUE)) {
        idx <- which(bad, arr.ind = TRUE)[1, , drop=TRUE]
        stop(sprintf("Invalid genotype '%s' at [%d,%d]. Expect 'AA','AB','BB' or missing.",
                     geno[idx[1], idx[2]], idx[1], idx[2]))
      }
      geno <- matrix(ifelse(is.na(geno), NA_integer_, unname(map[geno])),
                     nrow=nrow(geno), dimnames=dimnames(geno))
    }
  } else {
    storage.mode(geno) <- "integer"
  }

  # Validate value set
  bad_val <- !(is.na(geno) | geno %in% c(0L,1L,2L))
  if (any(bad_val)) {
    idx <- which(bad_val, arr.ind = TRUE)[1, , drop=TRUE]
    stop(sprintf("geno must contain only 0/1/2/NA. Found '%s' at [%d,%d].",
                 as.character(geno[idx[1], idx[2]]), idx[1], idx[2]))
  }

  # Encode to raw for SnpMatrix: 0x00 = missing (no call), and 0/1/2 copies of
  # the counted allele map to 0x01/0x02/0x03 respectively. This is the encoding
  # snpStats expects (and matches getGeno()); do NOT shift it.
  nr <- nrow(geno); nc <- ncol(geno)
  rawG <- matrix(as.raw(0L), nrow=nr, ncol=nc, dimnames=dimnames(geno))
  nz <- !is.na(geno)
  rawG[nz & geno == 0L] <- as.raw(1L)
  rawG[nz & geno == 1L] <- as.raw(2L)
  rawG[nz & geno == 2L] <- as.raw(3L)

  # Construct the SnpMatrix
  methods::new("SnpMatrix", rawG)
}
