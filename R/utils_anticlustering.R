if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("PC1", "PC2", "Group"))
}


# Internal: build the (optionally centered/scaled) numeric genotype matrix,
# dropping monomorphic SNPs. Kept as a matrix -- never a data.frame -- because
# for wide genotype data (hundreds of thousands of SNPs) data.frame operations
# are far slower and use much more memory.
#' @noRd
.scaledGenoMatrix <- function(object, center = FALSE, scale = FALSE) {
  if (!inherits(object, "SNPDataLong")) {
    stop("Input object must be of class SNPDataLong.")
  }

  snpsum <- snpStats::col.summary(object = object@geno)
  mono <- check.snp.monomorf(snpsum)
  if (is.null(mono)) mono <- character(0)  # in case your current checker returns NULL

  # drop monomorphic SNPs only if there are any
  if (length(mono) > 0) {
    object <- Subset(object = object, index = mono, margin = 2, keep = FALSE)
  } else {
    message("No monomorphic SNPs detected. Skipping subset.")
  }

  geno_matrix <- as(object@geno, "numeric")
  rownames(geno_matrix) <- rownames(object@geno)
  colnames(geno_matrix) <- colnames(object@geno)

  if (isTRUE(center) || isTRUE(scale) || is.numeric(center) || is.numeric(scale)) {
    message("Applying centering and/or scaling to SNP columns...")
    geno_matrix <- scale(geno_matrix, center = center, scale = scale)
    attr(geno_matrix, "scaled:center") <- NULL
    attr(geno_matrix, "scaled:scale") <- NULL
  }

  message("Genotype matrix prepared with dimensions: ",
          nrow(geno_matrix), " x ", ncol(geno_matrix))

  geno_matrix
}

#' Convert geno slot from SNPDataLong to a data.frame
#'
#' Converts the genotype matrix (geno slot) of a SNPDataLong object to a data.frame,
#' with optional centering and scaling per SNP (column).
#'
#' @param object An object of class SNPDataLong.
#' @param center Logical or numeric. If TRUE (default FALSE), center columns to mean zero.
#' @param scale Logical or numeric. If TRUE (default FALSE), scale columns to standard deviation one.
#'
#' @return A data.frame with individuals as rows and SNPs as columns (numeric 0/1/2, or centered/scaled values).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' raw_mat <- matrix(as.raw(sample(1:3, 100, TRUE)), nrow = 10, ncol = 10)
#' rownames(raw_mat) <- paste0("S", 1:10)
#' colnames(raw_mat) <- paste0("SNP", 1:10)
#' geno <- methods::new("SnpMatrix", raw_mat)
#' obj <- methods::new("SNPDataLong",
#'                     geno = geno,
#'                     map  = data.frame(Name = colnames(geno),
#'                                       Chromosome = 1,
#'                                       Position = 1:10),
#'                     path = tempfile(),
#'                     xref_path = "chip1")
#' df <- genoToDF(obj, center = TRUE, scale = TRUE)
#' head(df[, 1:5])
#' }
#' @export
genoToDF <- function(object, center = FALSE, scale = FALSE) {
  geno_matrix <- .scaledGenoMatrix(object, center = center, scale = scale)
  geno_df <- as.data.frame(geno_matrix)

  message("Genotype data converted to data.frame with dimensions: ",
          nrow(geno_df), " x ", ncol(geno_df))

  return(geno_df)
}

#' Run PCA on a SNPDataLong object
#'
#' Computes principal components of the (optionally centered/scaled) genotype
#' matrix, without any clustering. For wide data (more SNPs than individuals)
#' the PCA is obtained from the n x n Gram matrix, so the large rotation matrix
#' is never formed; when a fixed number of PCs is requested and \pkg{RSpectra}
#' is installed, only the top PCs are computed with a matrix-free solver. This
#' is the same PCA engine used by \code{\link{runAnticlusteringPCA}}, so the
#' scores are directly comparable.
#'
#' @param object An object of class \code{SNPDataLong}.
#' @param n_pcs Number of principal components to return. \code{NULL} (default)
#'   returns all PCs; a value \code{>= 1} returns that many and enables the fast
#'   \pkg{RSpectra} path; a value \code{< 1} is the proportion of variance to
#'   reach (e.g. \code{0.98}).
#' @param center Logical or numeric. Passed to \code{\link[base]{scale}}.
#'   If \code{TRUE}, center columns; if numeric, a vector of column means.
#'   Default: \code{TRUE}.
#' @param scale Logical or numeric. Passed to \code{\link[base]{scale}}.
#'   If \code{TRUE}, scale to unit variance; if numeric, a vector of column sds.
#'   Default: \code{TRUE}.
#'
#' @returns
#' A list with components:
#' \describe{
#'   \item{pca}{A \code{prcomp}-like object: \code{sdev}, \code{x} (scores) and
#'     \code{totvar} (total column variance). \code{rotation} is \code{NULL} for
#'     the wide-data paths.}
#'   \item{pcs}{Numeric matrix with the selected top principal components.}
#' }
#'
#' @examplesIf exists("nelore_imputed")
#' pr <- runPCA(nelore_imputed, n_pcs = 10)
#' head(pr$pcs[, 1:2])
#'
#' @export
#' @importFrom stats prcomp
runPCA <- function(object, n_pcs = NULL, center = TRUE, scale = TRUE) {
  if (!inherits(object, "SNPDataLong")) {
    stop("Input object must be of class SNPDataLong.")
  }
  if (!is.null(n_pcs) &&
      (!is.numeric(n_pcs) || length(n_pcs) != 1L || is.na(n_pcs))) {
    stop("`n_pcs` must be NULL or a single numeric value.")
  }

  geno_mat <- .scaledGenoMatrix(object, center = center, scale = scale)

  n <- nrow(geno_mat)
  p <- ncol(geno_mat)

  message("Running PCA...")
  # centering/scaling was already applied above, so no further centering here.
  if (p > n) {
    # Wide data (more SNPs than individuals). The scores/sdev of the SNP-space
    # PCA equal those of the eigendecomposition of the n x n Gram matrix
    # X X^T (X = geno_mat), so we never form the huge p x n rotation matrix.
    #
    # Fast path: when a fixed number of PCs is requested and RSpectra is
    # available, compute ONLY the top n_pcs eigenpairs with a matrix-free
    # operator v -> X (X^T v), avoiding the n x n Gram matrix entirely.
    fixed_k <- !is.null(n_pcs) && n_pcs >= 1
    use_rspectra <- fixed_k &&
      requireNamespace("RSpectra", quietly = TRUE) &&
      as.integer(n_pcs) <= (n - 2L)

    totvar <- sum(geno_mat^2) / max(1L, n - 1L)

    if (use_rspectra) {
      message("Using matrix-free truncated PCA via RSpectra (top ",
              as.integer(n_pcs), " PCs; p = ", p, " > n = ", n, ").")
      k    <- as.integer(n_pcs)
      Afun <- function(v, args) args %*% crossprod(args, v)  # X (X^T v)
      eig  <- RSpectra::eigs_sym(Afun, k = k, n = n, which = "LM",
                                 args = geno_mat)
      ev     <- pmax(eig$values, 0)
      sdev   <- sqrt(ev / max(1L, n - 1L))
      scores <- sweep(eig$vectors, 2, sqrt(ev), `*`)
    } else {
      message("Using Gram-matrix PCA (p = ", p, " > n = ", n, ").")
      gram <- tcrossprod(geno_mat)                 # n x n
      eig  <- eigen(gram, symmetric = TRUE)
      ev   <- pmax(eig$values, 0)                  # guard tiny negatives
      sdev   <- sqrt(ev / max(1L, n - 1L))
      scores <- sweep(eig$vectors, 2, sqrt(ev), `*`)  # U %*% diag(D) = prcomp$x
    }

    rownames(scores) <- rownames(geno_mat)
    colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
    # Mimic a prcomp object (rotation intentionally omitted). `totvar` keeps
    # variance-explained percentages correct even with only top-k sdev values.
    pca_res <- structure(
      list(sdev = sdev, rotation = NULL, center = FALSE, scale = FALSE,
           x = scores, totvar = totvar),
      class = "prcomp"
    )
  } else {
    pca_res <- stats::prcomp(geno_mat, center = FALSE, scale. = FALSE)
    pca_res$totvar <- sum(pca_res$sdev^2)
  }

  # Number of PCs to return. Denominator is the total column variance so
  # percentages stay correct even when only top-k sdev values are available.
  cum_var <- cumsum(pca_res$sdev^2 / pca_res$totvar)
  if (is.null(n_pcs)) {
    n_selected <- ncol(pca_res$x)
  } else if (n_pcs < 1) {
    n_selected <- which(cum_var >= n_pcs)[1]
    if (is.na(n_selected)) {
      stop("Could not reach requested variance proportion with available PCs.")
    }
    message("Automatically selecting ", n_selected,
            " PCs to explain at least ", round(n_pcs * 100, 1), "% variance.")
  } else {
    n_selected <- as.integer(n_pcs)
    message("Using fixed ", n_selected, " PCs.")
  }

  if (n_selected > ncol(pca_res$x)) {
    stop("Requested number of PCs (", n_selected,
         ") exceeds available PCs (", ncol(pca_res$x), ").")
  }

  list(pca = pca_res, pcs = pca_res$x[, seq_len(n_selected), drop = FALSE])
}

#' Run PCA and anticlustering on SNPDataLong
#'
#' Runs PCA on a \code{SNPDataLong} object (via \code{\link{runPCA}}) and then
#' performs anticlustering on the selected principal components.
#'
#' @param object An object of class \code{SNPDataLong}.
#' @param K Number of groups for anticlustering, or a vector of group sizes
#'   (as in \pkg{anticlust}).
#' @param n_pcs Number of top principal components to use. If \code{< 1}, it is
#'   interpreted as the proportion of variance to be explained (e.g. \code{0.8}).
#'   The fast matrix-free PCA path is used when a fixed number (\code{>= 1}) is
#'   requested and \pkg{RSpectra} is installed.
#' @param center Logical or numeric, passed to \code{\link{runPCA}}. Default
#'   \code{TRUE}.
#' @param scale Logical or numeric, passed to \code{\link{runPCA}}. Default
#'   \code{TRUE}.
#' @param anticlust_method Which \pkg{anticlust} optimiser to use.
#'   \code{"exchange"} (default) calls \code{anticlust::anticlustering};
#'   \code{"fast"} calls \code{anticlust::fast_anticlustering}, which scales to
#'   large numbers of individuals and may return different assignments.
#'
#' @returns
#' A list with components:
#' \describe{
#'   \item{groups}{Integer vector with anticlustering group assignments.}
#'   \item{pca}{The PCA result object (a \code{prcomp}-like list), as returned
#'     by \code{\link{runPCA}}.}
#'   \item{pcs}{Numeric matrix of the PCs used for anticlustering.}
#' }
#'
#' @examplesIf requireNamespace("anticlust", quietly = TRUE) && exists("nelore_imputed")
#' res <- runAnticlusteringPCA(nelore_imputed, K = 2, n_pcs = 0.8)
#' table(res$groups)
#'
#' @export
runAnticlusteringPCA <- function(object, K = 2, n_pcs = 20, center = TRUE,
                                 scale = TRUE,
                                 anticlust_method = c("exchange", "fast")) {
  anticlust_method <- match.arg(anticlust_method)

  if (!is.numeric(n_pcs) || length(n_pcs) != 1L || is.na(n_pcs)) {
    stop("`n_pcs` must be a single numeric value.")
  }

  pr      <- runPCA(object, n_pcs = n_pcs, center = center, scale = scale)
  pca_res <- pr$pca
  top_pcs <- pr$pcs
  message("Top PCs extracted.")

  if (!requireNamespace("anticlust", quietly = TRUE)) {
    stop("Package 'anticlust' is required. Please install it.")
  }

  message("Running anticlustering (", anticlust_method, ") with K = ",
          if (length(K) == 1) K else paste(K, collapse = ", "), " ...")
  if (anticlust_method == "fast") {
    groups <- anticlust::fast_anticlustering(top_pcs, K = K)
  } else {
    groups <- anticlust::anticlustering(
      x = top_pcs,
      K = K,
      standardize = TRUE
    )
  }

  message("Anticlustering completed. Groups assigned.")

  list(
    groups = groups,
    pca = pca_res,
    pcs = top_pcs
  )
}

#' Plot PCA groups from anticlustering result
#'
#' @param pca_res A prcomp object.
#' @param groups A factor or vector of group assignments.
#' @param pcs Vector of length 2 indicating which PCs to plot (default: c(1, 2)).
#' @param filename Optional. If provided, saves plot to this file (e.g., "antic.png").
#'
#' @return A ggplot object (also prints to screen).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pca_res <- stats::prcomp(matrix(rnorm(200), nrow = 20))
#' groups <- sample(1:2, 20, replace = TRUE)
#' plotPCAgroups(pca_res, groups)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal theme element_rect ggsave
#' @export
plotPCAgroups <- function(pca_res, groups, pcs = c(1, 2), filename = NULL) {
  # Use the stored total column variance when available so percentages stay
  # correct even when pca_res$sdev holds only the top-k PCs (truncated PCA).
  denom <- if (!is.null(pca_res$totvar)) pca_res$totvar else sum(pca_res$sdev^2)
  explained_var <- pca_res$sdev^2 / denom
  pc1_var <- round(100 * explained_var[pcs[1]], 2)
  pc2_var <- round(100 * explained_var[pcs[2]], 2)

  pc_df <- data.frame(
    PC1 = pca_res$x[, pcs[1]],
    PC2 = pca_res$x[, pcs[2]],
    Group = as.factor(groups)
  )

  p <- ggplot2::ggplot(pc_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = "PCA plot colored by Anticlustering Group",
      x = paste0("PC", pcs[1], " (", pc1_var, "%)"),
      y = paste0("PC", pcs[2], " (", pc2_var, "%)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))

  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 7, height = 5, dpi = 300)
    message("Plot saved to: ", filename)
  }

  return(p)
}
