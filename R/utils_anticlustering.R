if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("PC1", "PC2", "Group"))
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
  geno_df <- as.data.frame(geno_matrix)

  rownames(geno_df) <- rownames(object@geno)
  colnames(geno_df) <- colnames(object@geno)

  if (isTRUE(center) || isTRUE(scale) || is.numeric(center) || is.numeric(scale)) {
    message("Applying centering and/or scaling to SNP columns...")
    geno_df <- as.data.frame(scale(geno_df, center = center, scale = scale))
  }

  message("Genotype data converted to data.frame with dimensions: ",
          nrow(geno_df), " x ", ncol(geno_df))

  return(geno_df)
}

#' Run PCA and anticlustering on SNPDataLong
#'
#' Converts a SNPDataLong object to a data.frame, runs PCA, and performs
#' anticlustering on the selected principal components.
#'
#' @param object An object of class \code{SNPDataLong}.
#' @param K Number of groups for anticlustering, or a vector of group sizes
#'   (as in \pkg{anticlust}).
#' @param n_pcs Number of top principal components to use. If \code{< 1},
#'   it is interpreted as the proportion of variance to be explained (e.g.,
#'   \code{0.8} means PCs explaining at least 80\% variance).
#' @param center Logical or numeric. Passed to \code{\link[base]{scale}} via
#'   \code{genoToDF}. If \code{TRUE}, center columns; if numeric, a vector of
#'   column means. Default: \code{TRUE}.
#' @param scale Logical or numeric. Passed to \code{\link[base]{scale}} via
#'   \code{genoToDF}. If \code{TRUE}, scale to unit variance; if numeric,
#'   a vector of column sds. Default: \code{TRUE}.
#'
#' @returns
#' A list with components:
#' \describe{
#'   \item{groups}{Integer vector with anticlustering group assignments.}
#'   \item{pca}{The PCA result object (from \code{stats::prcomp}).}
#'   \item{pcs}{Numeric matrix of the PCs used for anticlustering.}
#' }
#'
#' @examplesIf requireNamespace("anticlust", quietly = TRUE) && exists("nelore_imputed")
#' res <- runAnticlusteringPCA(nelore_imputed, K = 2, n_pcs = 0.8)
#' table(res$groups)
#'
#' @export
#' @importFrom stats prcomp
runAnticlusteringPCA <- function(object, K = 2, n_pcs = 20, center = TRUE, scale = TRUE) {
  if (!inherits(object, "SNPDataLong")) {
    stop("Input object must be of class SNPDataLong.")
  }

  geno_df <- genoToDF(object, center = center, scale = scale)
  message("Genotype data frame created for PCA.")

  message("Running PCA...")
  # genoToDF already applied centering/scaling, so keep prcomp(center=FALSE, scale.=FALSE)
  pca_res <- stats::prcomp(geno_df, center = FALSE, scale. = FALSE)

  # Determine number of PCs to use
  var_explained <- pca_res$sdev^2
  var_explained <- var_explained / sum(var_explained)
  cum_var <- cumsum(var_explained)

  if (!is.numeric(n_pcs) || length(n_pcs) != 1L || is.na(n_pcs)) {
    stop("`n_pcs` must be a single numeric value.")
  }

  if (n_pcs < 1) {
    n_selected <- which(cum_var >= n_pcs)[1]
    if (is.na(n_selected)) {
      stop("Could not reach requested variance proportion with available PCs.")
    }
    message(
      "Automatically selecting ", n_selected,
      " PCs to explain at least ", round(n_pcs * 100, 1), "% variance."
    )
  } else {
    n_selected <- as.integer(n_pcs)
    message("Using fixed ", n_selected, " PCs.")
  }

  if (n_selected > ncol(pca_res$x)) {
    stop("Requested number of PCs (", n_selected, ") exceeds available PCs (", ncol(pca_res$x), ").")
  }

  top_pcs <- pca_res$x[, seq_len(n_selected), drop = FALSE]
  message("Top PCs extracted.")

  if (!requireNamespace("anticlust", quietly = TRUE)) {
    stop("Package 'anticlust' is required. Please install it.")
  }

  message("Running anticlustering with K = ",
          if (length(K) == 1) K else paste(K, collapse = ", "), " ...")
  groups <- anticlust::anticlustering(
    features = top_pcs,
    K = K,
    standardize = TRUE
  )

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
  explained_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
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
