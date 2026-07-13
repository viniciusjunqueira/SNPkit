# =========================================================================
# Functions from fQC package
# Original author: Roberto Higa <roberto.higa@embrapa.br>
# License: GPL-3
# Copied, integrated, and adapted into SNPkit by Vinicius Junqueira, 2025
# =========================================================================

#' Check SNP call rate
#'
#' Identifies SNPs with call rates below a minimum threshold.
#'
#' @param summary A data frame with SNP summary statistics (must contain `Call.rate` column).
#' @param min.call.rate Numeric value specifying the minimum acceptable call rate.
#'
#' @return Character vector with SNP names below threshold. Returns `NULL` if none.
#'
#' @examples
#' df <- data.frame(Call.rate = c(0.85, 0.95), row.names = c("SNP1", "SNP2"))
#' check.call.rate(df, 0.9)
#'
#' @export
check.call.rate <- function(summary, min.call.rate) {
  result <- summary$Call.rate < min.call.rate
  result[is.na(result)] <- FALSE
  names <- NULL
  if (sum(result) > 0) {
    names <- rownames(summary[result, ])
  }
  return(names)
}

#' Check Identity-By-State (IBS) for a genotype pair
#'
#' Checks IBS status for two genotypes.
#'
#' @param gen Numeric vector of length two with genotype codes.
#'
#' @return Integer: 2 if identical non-heterozygotes, 0 if opposite homozygotes, -1 otherwise.
#'
#' @examples
#' check.ibs(c(1, 1))
#' check.ibs(c(1, 3))
#'
#' @export
check.ibs <- function(gen) {
  ret <- -1
  if (gen[1] != 2 && gen[1] == gen[2]) {
    ret <- 2
  } else if (gen[1] == 1 && gen[2] == 3) {
    ret <- 0
  } else if (gen[1] == 3 && gen[2] == 1) {
    ret <- 0
  }
  return(ret)
}

#' Check identical samples based on distance
#'
#' Identifies sample pairs considered identical based on genotype distances.
#'
#' @param genotypes Genotype matrix (samples x SNPs) or SnpMatrix.
#' @param threshold Numeric distance threshold. Default 0.
#'
#' @return Data frame of identical sample pairs.
#'
#' @examples
#' mat <- matrix(sample(0:2, 20, TRUE), nrow = 5)
#' rownames(mat) <- paste0("S", 1:5)
#' check.identical.samples(mat, 0.5)
#'
#' @importFrom methods as
#' @export
check.identical.samples <- function(genotypes, threshold = 0) {
  if (inherits(genotypes, "SnpMatrix")) {
    numeric_geno <- as(genotypes, "numeric")
  } else {
    numeric_geno <- genotypes
  }

  empty <- data.frame(Sample1 = character(),
                      Sample2 = character(),
                      Distance = numeric(),
                      stringsAsFactors = FALSE)

  n <- nrow(numeric_geno)
  if (is.null(n) || n < 2) {
    return(empty)
  }

  mdistm <- as.matrix(stats::dist(numeric_geno))
  sample.names <- rownames(mdistm)

  # Vectorized extraction of the upper-triangle pairs within threshold, ordered
  # row-major (i < j) -- avoids the quadratic rbind-in-loop.
  hits <- which(upper.tri(mdistm) & mdistm <= threshold, arr.ind = TRUE)
  if (nrow(hits) == 0) {
    return(empty)
  }
  hits <- hits[order(hits[, "row"], hits[, "col"]), , drop = FALSE]

  sample.pairs <- data.frame(Sample1 = sample.names[hits[, "row"]],
                             Sample2 = sample.names[hits[, "col"]],
                             Distance = mdistm[hits],
                             stringsAsFactors = FALSE)

  for (r in seq_len(nrow(sample.pairs))) {
    warning(paste("Identical samples:", sample.pairs$Sample1[r],
                  "-", sample.pairs$Sample2[r]))
  }
  sample.pairs
}


#' Check identical samples by block
#'
#' Identifies sample pairs that stay identical (within \code{threshold}) across
#' \emph{every} SNP block, scanning the markers in blocks of \code{blcsize}.
#' Each block only re-checks the samples still in a confirmed pair, and pairs
#' that separate in any block are dropped, so the result is the intersection of
#' the per-block identical pairs.
#'
#' @param genotypes Genotype matrix (samples x SNPs) or SnpMatrix with sample
#'   names as rownames.
#' @param blcsize Block size (number of SNPs).
#' @param threshold Distance threshold. Default 0.
#'
#' @return A data.frame of identical sample pairs (columns \code{Sample1},
#'   \code{Sample2}, \code{Distance}); \code{Distance} is taken from the first
#'   block. Empty data.frame if none.
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(sample(1:3, 40, TRUE), nrow = 4)
#' rownames(mat) <- paste0("S", 1:4)
#' check.identical.samples.by.block(mat, blcsize = 5, threshold = 0)
#'
#' @export
check.identical.samples.by.block <- function(genotypes, blcsize, threshold = 0) {
  empty <- data.frame(Sample1 = character(), Sample2 = character(),
                      Distance = numeric(), stringsAsFactors = FALSE)
  n_snp <- ncol(genotypes)
  if (is.null(n_snp) || n_snp < 1 || nrow(genotypes) < 2) {
    return(empty)
  }

  # Order-independent key so {A,B} and {B,A} match across blocks.
  pair_ids <- function(df) {
    if (nrow(df) == 0) return(character(0))
    paste(pmin(df$Sample1, df$Sample2), pmax(df$Sample1, df$Sample2), sep = "\r")
  }

  confirmed <- NULL   # pairs identical in every block processed so far
  ini <- 1
  while (ini <= n_snp) {
    fin <- min(ini + blcsize - 1, n_snp)
    message("Analyzing block ", ini, "-", fin)

    # Only samples still in a confirmed pair need re-checking.
    rows <- if (is.null(confirmed)) rownames(genotypes)
            else unique(c(confirmed$Sample1, confirmed$Sample2))
    block <- suppressWarnings(
      check.identical.samples(genotypes[rows, ini:fin, drop = FALSE], threshold)
    )

    if (is.null(confirmed)) {
      confirmed <- block
    } else {
      confirmed <- confirmed[pair_ids(confirmed) %in% pair_ids(block), , drop = FALSE]
    }
    if (nrow(confirmed) == 0) {
      return(empty)
    }
    ini <- fin + 1
  }
  rownames(confirmed) <- NULL
  confirmed
}

#' Check Mendelian inconsistencies
#'
#' Identifies Mendelian inconsistencies between father-child pairs.
#'
#' @param genotypes Genotype matrix.
#' @param father Vector of father sample IDs.
#' @param child Vector of child sample IDs.
#'
#' @return Data frame summarizing inconsistencies per pair.
#'
#' @examples
#' set.seed(1)
#' genotypes <- matrix(sample(1:3, 30, TRUE), nrow = 3,
#'                     dimnames = list(c("F1", "C1", "C2"), NULL))
#' check.mendelian.inconsistencies(genotypes,
#'                                 father = "F1",
#'                                 child  = c("C1", "C2"))
#'
#' @export
check.mendelian.inconsistencies <- function(genotypes, father, child) {
  m <- length(child)
  n <- length(father)
  if (n == 0 || m == 0) {
    return(data.frame(Father = character(), Child = character(),
                      N = numeric(), Total = numeric(), Rate = numeric(),
                      stringsAsFactors = FALSE))
  }
  sample1 <- NULL
  sample2 <- NULL
  n.inconsist <- NULL
  t.inconsist <- NULL
  tx.inconsist <- NULL
  for (i in seq_len(n)) {
    g1 <- genotypes[father[i], ]
    nam1 <- father[i]
    for (j in seq_len(m)) {
      nam2 <- child[j]
      if (nam1 != nam2) {
        sample1 <- c(sample1, paste(nam1))
        g2 <- genotypes[child[j], ]
        sample2 <- c(sample2, nam2)
        counts <- check.mendelian.inconsistencies.pair(g1, g2)
        n.inconsist <- c(n.inconsist, counts[1])
        t.inconsist <- c(t.inconsist, counts[2])
        tx.i <- counts[1] / counts[2]
        tx.inconsist <- c(tx.inconsist, tx.i)
        message(nam1, " - ", nam2, " = ", counts[1], " ", counts[2], " ", tx.i)
      }
    }
  }
  result <- data.frame(sample1 = sample1, sample2 = sample2, inconsist = n.inconsist, total = t.inconsist, rate = tx.inconsist)
  colnames(result) <- c("Father", "Child", "N", "Total", "Rate")
  return(result)
}

#' Check Mendelian inconsistencies for a pair
#'
#' Calculates number of inconsistencies and total comparable SNPs for a parent-child pair.
#'
#' @param g1 Genotype vector for parent.
#' @param g2 Genotype vector for child.
#'
#' @return Numeric vector: [# inconsistencies, # comparable SNPs].
#'
#' @examples
#' g1 <- c(1, 1, 3, 3, 2)
#' g2 <- c(3, 1, 1, 3, 2)
#' check.mendelian.inconsistencies.pair(g1, g2)
#'
#' @export
check.mendelian.inconsistencies.pair <- function(g1, g2) {
  inconsist <- (g1 == 1 & g2 == 3) | (g1 == 3 & g2 == 1)
  homoz2 <- (g1 == 1 | g1 == 3) & (g2 == 1 | g2 == 3)
  ret <- c(sum(inconsist), sum(homoz2))
  return(ret)
}

#' Check sample heterozygosity
#'
#' Identifies samples with heterozygosity values deviating beyond a specified threshold.
#'
#' @param sample.summary Data frame containing sample summary (must have `Heterozygosity` column).
#' @param max.dev Maximum number of standard deviations allowed from mean.
#'
#' @return Character vector with sample names considered outliers. Returns `NULL` if none.
#'
#' @examples
#' ss <- data.frame(Heterozygosity = c(0.2, 0.5, 0.7))
#' rownames(ss) <- c("Ind1", "Ind2", "Ind3")
#' check.sample.heterozygosity(ss, 1)
#'
#' @importFrom stats sd
#' @export
check.sample.heterozygosity <- function(sample.summary, max.dev) {
  m <- mean(sample.summary[, "Heterozygosity"], na.rm = TRUE)
  s <- sd(sample.summary[, "Heterozygosity"], na.rm = TRUE)

  below <- sample.summary[, "Heterozygosity"] < m - max.dev * s
  below[is.na(below)] <- FALSE

  above <- sample.summary[, "Heterozygosity"] > m + max.dev * s
  above[is.na(above)] <- FALSE

  smps <- union(rownames(sample.summary)[below], rownames(sample.summary)[above])

  if (length(smps) == 0) {
    return(NULL)
  }

  return(smps)
}


#' Check SNP by chromosome
#'
#' Filters SNP names belonging to specified chromosomes.
#'
#' @param snpmap Data frame with SNP map info (must contain columns `Chromosome` and `Name`).
#' @param chromosomes Vector of chromosome identifiers to filter.
#'
#' @return Character vector with SNP names.
#'
#' @examples
#' snpmap <- data.frame(Chromosome = c(1, 1, 2), Name = c("SNP1", "SNP2", "SNP3"))
#' check.snp.chromo(snpmap, 1)
#'
#' @export
check.snp.chromo <- function(snpmap, chromosomes) {
  snps <- snpmap[snpmap$Chromosome %in% chromosomes, "Name"]
  if (length(snps) == 0) {
    snps <- NULL
  }
  return(as.character(snps))
}

#' Check SNP Hardy-Weinberg equilibrium deviation
#'
#' Identifies SNPs deviating from HWE beyond a z-score threshold.
#'
#' @param snp.summary Data frame with SNP summary (must contain `z.HWE` column).
#' @param max.dev Maximum z-score allowed.
#'
#' @return Character vector with SNP names deviating from HWE. Returns `NULL` if none.
#'
#' @examples
#' df <- data.frame(z.HWE = c(2, 5), row.names = c("SNP1", "SNP2"))
#' check.snp.hwe(df, 3)
#'
#' @export
check.snp.hwe <- function(snp.summary, max.dev) {
  result <- snp.summary$z.HWE^2 >= max.dev^2
  result[is.na(result)] <- FALSE
  snps <- NULL
  if (sum(result) > 0) {
    snps <- rownames(snp.summary[result, ])
  }
  return(snps)
}

#' Check SNP minor allele frequency
#'
#' Identifies SNPs with minor allele frequency below a minimum threshold.
#'
#' @param snp.summary Data frame with SNP summary (must contain `MAF` column).
#' @param min.maf Minimum MAF allowed.
#'
#' @return Character vector with SNP names below threshold. Returns `NULL` if none.
#'
#' @examples
#' df <- data.frame(MAF = c(0.01, 0.2), row.names = c("SNP1", "SNP2"))
#' check.snp.maf(df, 0.05)
#'
#' @export
check.snp.maf <- function(snp.summary, min.maf) {
  result <- snp.summary$MAF < min.maf
  result[is.na(result)] <- FALSE
  snps <- NULL
  if (sum(result) > 0) {
    snps <- rownames(snp.summary[result, ])
  }
  return(snps)
}

#' Check SNP missing genotype frequencies
#'
#' Identifies SNPs with genotype frequencies below a minimum threshold.
#'
#' @param snp.summary Data frame with columns `P.AA`, `P.AB`, `P.BB`.
#' @param min.mgf Minimum genotype frequency allowed.
#'
#' @return Character vector with SNP names below threshold. Returns `NULL` if none.
#'
#' @examples
#' df <- data.frame(P.AA = c(0.01, 0.5), P.AB = c(0.02, 0.4), P.BB = c(0.01, 0.1))
#' rownames(df) <- c("SNP1", "SNP2")
#' check.snp.mgf(df, 0.05)
#'
#' @export
check.snp.mgf <- function(snp.summary, min.mgf) {
  result <- snp.summary$P.AA < min.mgf | snp.summary$P.AB < min.mgf | snp.summary$P.BB < min.mgf
  result[is.na(result)] <- FALSE
  mgf <- NULL
  if (sum(result) > 0) {
    mgf <- rownames(snp.summary[result, ])
  }
  return(mgf)
}

#' Check SNP monomorphic status
#'
#' Identifies SNPs considered monomorphic.
#'
#' @param snp.summary Data frame with columns `P.AA`, `P.AB`, `P.BB`.
#'
#' @return Character vector with monomorphic SNP names. Returns `NULL` if none.
#'
#' @examples
#' df <- data.frame(P.AA = c(1, 0.5), P.AB = c(0, 0.5), P.BB = c(0, 0))
#' rownames(df) <- c("SNP1", "SNP2")
#' check.snp.monomorf(df)
#'
#' @export
check.snp.monomorf <- function(snp.summary) {
  result <- snp.summary$P.AA == 1 | snp.summary$P.AB == 1 | snp.summary$P.BB == 1
  result[is.na(result)] <- FALSE
  snps <- NULL
  if (sum(result) > 0) {
    snps <- rownames(snp.summary[result, ])
  }
  return(snps)
}

#' Check SNP no position
#'
#' Identifies SNPs without a usable genomic position, i.e. whose position is
#' missing (`NA`), blank, non-numeric, or zero. The `Position` column may be
#' numeric or character (`getGeno()` reads maps as character), so it is coerced
#' to numeric first.
#'
#' @param snpmap Data frame with columns `Position` and `Name`.
#'
#' @return Character vector with SNP names without position. Returns `NULL` if none.
#'
#' @examples
#' df <- data.frame(Position = c(0, 100, NA), Name = c("SNP1", "SNP2", "SNP3"))
#' check.snp.no.position(df)  # SNP1 (zero) and SNP3 (missing)
#'
#' @export
check.snp.no.position <- function(snpmap) {
  # Coerce first so numeric and character maps behave the same; blank/non-numeric
  # values become NA and count as "no position", as does an explicit zero.
  pos_num <- suppressWarnings(as.numeric(snpmap[["Position"]]))
  no_pos  <- is.na(pos_num) | pos_num == 0
  snps <- as.character(snpmap[["Name"]])[no_pos]
  if (length(snps) == 0) {
    return(NULL)
  }
  snps
}

#' Check SNPs with same position
#'
#' Identifies SNPs that share the same position on the same chromosome.
#'
#' @param snpmap Data frame with columns `Chromosome`, `Position`, and `Name`.
#'
#' @return List of SNP groups sharing positions.
#'
#' @examples
#' df <- data.frame(Chromosome = c(1, 1, 2),
#'                  Position = c(100, 100, 200),
#'                  Name = c("SNP1", "SNP2", "SNP3"))
#' check.snp.same.position(df)
#'
#' @export
check.snp.same.position <- function(snpmap) {
  pos <- snpmap[["Position"]]
  ok  <- !is.na(pos)
  if (!any(ok)) {
    return(list())
  }
  # Group SNP names by chromosome + position; groups with more than one SNP
  # share a locus. Vectorized and free of adjacency/index bookkeeping, so it is
  # safe for single-SNP chromosomes and missing positions.
  key    <- paste(snpmap[["Chromosome"]][ok], pos[ok], sep = ":")
  snp_id <- as.character(snpmap[["Name"]])[ok]
  groups <- split(snp_id, key)
  unname(groups[lengths(groups) > 1])
}


#' IBS pair statistics
#'
#' Calculates IBS mean and standard deviation between two samples.
#'
#' @param g1 Genotype vector for first sample.
#' @param g2 Genotype vector for second sample.
#'
#' @return Numeric vector: [mean IBS, standard deviation].
#'
#' @examples
#' g1 <- sample(0:2, 10, TRUE)
#' g2 <- sample(0:2, 10, TRUE)
#' ibs.pair(g1, g2)
#'
#' @export
ibs.pair <- function(g1, g2) {
  mat <- rbind(g1, g2)
  vet <- apply(mat, 2, check.ibs)
  vet[vet < 0] <- mean(vet[vet >= 0])
  m <- mean(vet)
  s <- sd(vet)
  return(c(m, s))
}

#' Convert pairs to sets
#'
#' Groups sample pairs into sets of related samples.
#'
#' @param pairs Matrix or list of sample pairs.
#'
#' @return List of sets of samples.
#'
#' @examples
#' pairs <- matrix(c("A", "B", "B", "C", "D", "E"), ncol = 2, byrow = TRUE)
#' pairs2sets(pairs)
#'
#' @export
pairs2sets <- function(pairs) {
  if (length(pairs) > 0) {
    sample.pairs <- matrix(pairs[, 1:2], ncol = 2)
    idx <- 1:dim(sample.pairs)[1]
    n <- length(idx)
    k <- 1
    sample.ident <- list()
    while (n > 0) {
      toremove <- idx[1]
      sample.ident[[k]] <- as.character(sample.pairs[idx[1], ])
      pivot <- sample.ident[[k]]
      if (n >= 2) {
        for (i in 2:n) {
          settest <- as.character(sample.pairs[idx[i], ])
          if (length(intersect(pivot, settest)) > 0) {
            pivot <- union(pivot, settest)
            toremove <- c(toremove, idx[i])
          }
        }
      }
      sample.ident[[k]] <- pivot
      k <- k + 1
      idx <- setdiff(idx, toremove)
      n <- length(idx)
    }
    return(sample.ident)
  }
}

#' Do genome relationship matrix PCA (deprecated)
#'
#' @description
#' \strong{Deprecated}. Performs PCA using the genome relationship matrix (GRM)
#' on a raw \code{SnpMatrix}. Use \code{\link{runPCA}} instead, which operates on
#' a \code{SNPDataLong} object, standardises SNPs, and returns scores directly
#' comparable to \code{\link{runAnticlusteringPCA}}.
#'
#' @param genotypes Genotype matrix.
#'
#' @return List containing `pcs` (principal components) and `eigen` (eigenvalues).
#'
#' @seealso \code{\link{runPCA}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' mat <- matrix(sample(as.raw(1:3), 200, TRUE), nrow = 10, ncol = 20)
#' geno <- methods::new("SnpMatrix", mat)
#' rownames(geno) <- paste0("S", 1:10)
#' colnames(geno) <- paste0("SNP", 1:20)
#' res <- suppressWarnings(doPCA(geno))  # doPCA is deprecated; use runPCA()
#' str(res)
#' }
#'
#' @export
doPCA <- function(genotypes) {
  .Deprecated("runPCA")
  xxmat <- snpStats::xxt(genotypes, correct.for.missing = FALSE)
  evv <- eigen(xxmat, symmetric = TRUE)
  pcs <- evv$vectors
  evals <- evv$values
  message("Eigenvalues near zero set to zero (|eigenvalue| < 1e-3)")
  evals[abs(evals) < 0.001] <- 0
  btr <- snpStats::snp.pre.multiply(genotypes, diag(1/sqrt(evals)) %*% t(pcs))
  pcs <- snpStats::snp.post.multiply(genotypes, t(btr))
  return(list(pcs = pcs, eigen = evals))
}

#' Exploratory plots for SNP and sample summary
#'
#' Generates exploratory plots: MAF histograms, HWE plots, heterozygosity scatter, MDS, and dendrogram.
#'
#' @param snp.summary Data frame with SNP summary.
#' @param snps.plot Filename for SNP histogram plot.
#' @param sample.summary Data frame with sample summary.
#' @param samples.plot Filename for heterozygosity plot.
#' @param distm Distance matrix for samples.
#' @param glabels Sample labels for plots.
#' @param mds.plot Filename for MDS plot.
#' @param hierq.plot Filename for hierarchical cluster plot.
#'
#' @return NULL (plots are saved as JPEG files)
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".jpg")
#' snp.summary <- data.frame(
#'   MAF   = runif(20),
#'   z.HWE = rnorm(20),
#'   Calls = rep(100, 20),
#'   P.AA  = runif(20, 0, 0.5),
#'   P.AB  = runif(20, 0, 0.5),
#'   P.BB  = runif(20, 0, 0.5)
#' )
#' sample.summary <- data.frame(
#'   Call.rate      = runif(5, 0.9, 1),
#'   Heterozygosity = runif(5, 0.2, 0.4),
#'   row.names = paste0("S", 1:5)
#' )
#' distm <- stats::dist(matrix(rnorm(25), nrow = 5))
#' exploratory.plots(snp.summary,
#'                   snps.plot      = tempfile(fileext = ".jpg"),
#'                   sample.summary = sample.summary,
#'                   samples.plot   = tempfile(fileext = ".jpg"),
#'                   distm          = distm,
#'                   glabels        = paste0("S", 1:5),
#'                   mds.plot       = tempfile(fileext = ".jpg"),
#'                   hierq.plot     = tempfile(fileext = ".jpg"))
#' }
#'
#' @importFrom grDevices jpeg dev.off
#' @importFrom graphics par hist text plot
#' @importFrom MASS isoMDS
#' @importFrom stats hclust
#' @export
exploratory.plots <- function(snp.summary, snps.plot, sample.summary, samples.plot, distm, glabels, mds.plot, hierq.plot) {
  # Helper to safely open and close JPEG devices and restore par() afterwards
  safe_jpeg <- function(file, expr, ...) {
    grDevices::jpeg(file, ...)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit({
      graphics::par(oldpar)
      grDevices::dev.off()
    }, add = TRUE)
    force(expr)
  }

  # SNP histograms (MAF and HWE)
  safe_jpeg(snps.plot, {
    graphics::par(mfrow = c(1, 2))
    graphics::hist(snp.summary$MAF, main = "Histogram of MAF", xlab = "MAF", col = "grey")
    graphics::hist(snp.summary$z.HWE, main = "Histogram of HWE (z-score)", xlab = "HWE z-score", col = "grey")
  }, width = 900, height = 600)

  # Histogram of HWE p-values (chi2)
  base_sample <- tools::file_path_sans_ext(samples.plot)
  safe_jpeg(paste0(base_sample, ".chi2.jpg"), {
    pvchi2 <- get.hwe.chi2(snp.summary)
    graphics::hist(pvchi2, main = "Histogram of HWE (Chi2 p-values)", xlab = "HWE p-value", col = "grey")
  }, width = 800, height = 600)

  # Call rate vs heterozygosity
  safe_jpeg(samples.plot, {
    graphics::par(mfrow = c(1, 1))
    graphics::plot(sample.summary$Call.rate, sample.summary$Heterozygosity,
                   xlab = "Call rate", ylab = "Heterozygosity", main = "Call rate vs Heterozygosity",
                   pch = 19, col = "blue")
  }, width = 800, height = 600)

  # Same plot with labels
  safe_jpeg(paste0(base_sample, ".1.jpg"), {
    graphics::plot(sample.summary$Call.rate, sample.summary$Heterozygosity,
                   xlab = "Call rate", ylab = "Heterozygosity", main = "Call rate vs Heterozygosity (labeled)",
                   pch = 19, col = "blue")
    graphics::text(sample.summary$Call.rate, sample.summary$Heterozygosity, labels = rownames(sample.summary), pos = 3, cex = 0.7)
  }, width = 1000, height = 800)

  # MDS plot
  iso <- MASS::isoMDS(distm, tol = 1e-10, maxit = 500)

  safe_jpeg(mds.plot, {
    graphics::plot(iso$points[, 1], iso$points[, 2],
                   xlab = "Dim 1", ylab = "Dim 2", main = "Samples MDS",
                   pch = 19, col = "darkgreen")
  }, width = 800, height = 600)

  # MDS with labels
  base_mds <- tools::file_path_sans_ext(mds.plot)
  safe_jpeg(paste0(base_mds, ".1.jpg"), {
    graphics::plot(iso$points[, 1], iso$points[, 2],
                   xlab = "Dim 1", ylab = "Dim 2", main = "Samples MDS (labeled)",
                   pch = 19, col = "darkgreen")
    graphics::text(iso$points[, 1], iso$points[, 2], labels = glabels, pos = 3, cex = 0.7)
  }, width = 1000, height = 800)

  # Hierarchical clustering
  safe_jpeg(hierq.plot, {
    hcl <- stats::hclust(distm, method = "single")
    graphics::plot(hcl, main = "Hierarchical Cluster", xlab = "Samples", ylab = "Distances")
  }, width = 800, height = 600)

  invisible(NULL)
}


#' Get correlation (fc method)
#'
#' Calculates genotype correlation using a fast check (fc) method.
#'
#' @param g1 Genotype vector.
#' @param g2 Genotype vector.
#'
#' @return Numeric value of correlation.
#'
#' @examples
#' g1 <- sample(0:2, 10, TRUE)
#' g2 <- sample(0:2, 10, TRUE)
#' get.correl.fc(g1, g2)
#'
#' @export
get.correl.fc <- function(g1, g2) {
  g1 <- as.raw(g1)
  g2 <- as.raw(g2)
  # av marks positions called in both samples (non-zero). Because it already
  # excludes zeros, the former `== 0` concordance term was always zero and is
  # dropped here without changing the result.
  av <- as.logical(g1) & as.logical(g2)
  t1 <- sum(av)
  t2 <- sum(g1[av] == 1 & g2[av] == 1) + sum(g1[av] == 2 & g2[av] == 2)
  return(ifelse(t1, t2 / t1, 0))
}

#' Get gender based on heterozygosity
#'
#' Infers gender using heterozygosity thresholds.
#'
#' @param sample.summary Data frame with `Heterozygosity` column.
#' @param threshM Numeric threshold for males.
#' @param threshF Numeric threshold for females.
#'
#' @return Data frame with columns `heterozygosity` and `sex`.
#'
#' @examples
#' df <- data.frame(Heterozygosity = c(0.1, 0.3, 0.6))
#' rownames(df) <- c("A", "B", "C")
#' get.gender(df, 0.2, 0.5)
#'
#' @export
get.gender <- function(sample.summary, threshM, threshF) {
  if (threshM > threshF | threshM <= 0 | threshF <= 0) {
    stop("Invalid thresholds.")
  }
  h <- sample.summary$Heterozygosity
  sex <- rep("I", length(h))
  sex[h < threshM] <- "M"
  sex[h >= threshF] <- "F"
  ret <- data.frame(heterozygosity = h, sex = sex)
  rownames(ret) <- rownames(sample.summary)
  return(ret)
}

#' Get HWE chi-square p-values
#'
#' Calculates Hardy-Weinberg equilibrium chi-square p-values for SNPs.
#'
#' @param snp.summary Data frame with columns `Calls`, `P.AA`, `P.AB`, `P.BB`.
#'
#' @return Numeric vector with p-values.
#'
#' @examples
#' df <- data.frame(Calls = c(100, 100), P.AA = c(0.6, 0.4), P.AB = c(0.3, 0.4), P.BB = c(0.1, 0.2))
#' get.hwe.chi2(df)
#'
#' @importFrom stats pchisq
#' @export
get.hwe.chi2 <- function(snp.summary) {
  # Observed counts
  ObsCountAA <- snp.summary$Calls * snp.summary$P.AA
  ObsCountAB <- snp.summary$Calls * snp.summary$P.AB
  ObsCountBB <- snp.summary$Calls * snp.summary$P.BB

  # Allele frequency A
  freqA <- (2 * snp.summary$P.AA + snp.summary$P.AB) / 2

  # Expected counts under HWE
  ExpCountAA <- snp.summary$Calls * freqA^2
  ExpCountAB <- 2 * snp.summary$Calls * freqA * (1 - freqA)
  ExpCountBB <- snp.summary$Calls * (1 - freqA)^2

  # Avoid division by zero
  ExpCountAA[ExpCountAA == 0] <- 1e-6
  ExpCountAB[ExpCountAB == 0] <- 1e-6
  ExpCountBB[ExpCountBB == 0] <- 1e-6

  # Chi-square statistic
  chi2stat <- (ObsCountAA - ExpCountAA)^2 / ExpCountAA +
    (ObsCountAB - ExpCountAB)^2 / ExpCountAB +
    (ObsCountBB - ExpCountBB)^2 / ExpCountBB

  # P-values
  pvalues <- pchisq(chi2stat, df = 1, lower.tail = FALSE)
  return(pvalues)
}


#' Check SNPs for Hardy-Weinberg equilibrium deviation using chi-square p-values
#'
#' This function identifies SNP markers whose Hardy-Weinberg equilibrium (HWE) chi-square p-values
#' indicate significant deviation beyond a specified threshold. It uses the p-values computed by
#' \code{get.hwe.chi2} on the input summary data frame.
#'
#' @param snp.summary A data frame or matrix containing summary statistics for SNP markers.
#'        The row names should correspond to SNP identifiers. It must be compatible with
#'        the function \code{get.hwe.chi2}.
#' @param max.dev A numeric value specifying the maximum acceptable p-value threshold.
#'        SNPs with p-values below this threshold are considered as deviating from HWE.
#'
#' @return A character vector of SNP identifiers (rownames) that fail the HWE test (p-value < \code{max.dev}).
#'         If no SNPs fail, an empty vector is returned.
#'
#' @details Any SNP with missing p-value (NA) is treated as not failing (returned as FALSE).
#'
#' @seealso \code{\link{get.hwe.chi2}}
#'
#' @examples
#' snp.summary <- data.frame(
#'   Calls = c(100, 100),
#'   P.AA  = c(0.25, 0.7),
#'   P.AB  = c(0.50, 0.05),
#'   P.BB  = c(0.25, 0.25),
#'   row.names = c("SNP1", "SNP2")
#' )
#' check.snp.hwe.chi2(snp.summary, max.dev = 0.05)
#'
#' @export
check.snp.hwe.chi2 <- function (snp.summary, max.dev)
{
    pvalues <- get.hwe.chi2(snp.summary)
    result <- pvalues < max.dev
    result[is.na(result)] <- FALSE
    snps <- NULL
    if (sum(result) > 0) {
        snps <- rownames(snp.summary[result, ])
    }
    return(snps)
}
