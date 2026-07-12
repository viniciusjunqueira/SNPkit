#' Save SNPDataLong object to PLINK format
#'
#' Saves genotype and map data from an SNPDataLong object in PLINK format (.ped/.map and optionally binary files).
#'
#' @param object An object of class SNPDataLong.
#' @param path Character. Directory where files will be saved. Must be supplied
#'   by the caller (e.g. a folder inside \code{tempdir()} for examples).
#' @param name Character. Base name for PLINK output files.
#' @param run_plink Logical. If TRUE (default), runs PLINK1 to convert to binary files. If FALSE, only .ped and .map files are saved.
#' @param chunk_size Integer. Number of individuals per chunk for writing .ped file (default: 1000).
#' @param extra_args Character vector. Extra arguments appended verbatim to the
#'   PLINK command line when \code{run_plink = TRUE}. Useful for non-human data;
#'   e.g. \code{"--chr-set 29"} makes PLINK treat autosomes 1-29 correctly for
#'   cattle instead of reading 23-26 as human X/Y/XY/MT. Default is NULL.
#'
#' @return No return value, called for side effects. Files (\code{.ped}/\code{.map},
#'   and \code{.bed}/\code{.bim}/\code{.fam} when \code{run_plink = TRUE}) are
#'   written under \code{path}.
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
#' savePlink(obj, path = tempdir(), name = "demo",
#'           run_plink = FALSE, chunk_size = 5)
#' }
#' @importFrom utils write.table
#' @export
savePlink <- function(object, path, name = "plink_data", run_plink = TRUE, chunk_size = 1000, extra_args = NULL) {
  if (!inherits(object, "SNPDataLong")) {
    stop("Input object must be of class SNPDataLong.")
  }

  if (missing(path) || !is.character(path) || length(path) != 1) {
    stop("'path' must be a single character string indicating the output directory.")
  }

  qc_header("Saving Files in Plink Format")

  geno <- object@geno
  map <- object@map
  n_ind <- nrow(geno)

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created output directory: ", path)
  }

  ## ----- PED file -----
  ped_file <- file.path(path, paste0(name, ".ped"))
  smp <- rownames(geno)

  message("Writing .ped file in chunks...")
  con <- file(ped_file, "wt")

  for (start in seq(1, n_ind, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n_ind)
    idx <- start:end

    # Convert block of individuals
    geno_chr <- as(geno[idx, , drop = FALSE], "character")
    geno_chr <- gsub("/", " ", geno_chr)
    geno_chr[is.na(geno_chr)] <- "0 0"

    # Build lines for this chunk
    lines <- vapply(seq_along(idx), function(i) {
      paste("NA", smp[idx[i]], "0 0 -9 -9", paste(geno_chr[i, ], collapse = " "))
    }, character(1L))

    writeLines(lines, con)
    message(sprintf(" Wrote individuals %d to %d", start, end))
  }

  close(con)
  message(".ped file written: ", ped_file)

  ## ----- MAP file -----
  message("Writing .map file...")
  map_file <- file.path(path, paste0(name, ".map"))
  map_out <- data.frame(
    Chromosome = map$Chromosome,
    SNP_ID = map$Name,
    Position = map$Position,
    stringsAsFactors = FALSE
  )
  utils::write.table(map_out, map_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  message(".map file written: ", map_file)

  ## ----- Optionally run PLINK -----
  if (run_plink) {
    message("Running PLINK to generate binary files...")
    extra <- if (!is.null(extra_args)) paste(extra_args, collapse = " ") else ""
    cmd <- paste("cd", shQuote(path), "&& plink1 --file", shQuote(name), "--map3 --out", shQuote(name), "--make-bed --noweb", extra)
    status <- system(cmd)

    if (status == 0) {
      message("PLINK binary files created successfully.")
    } else {
      warning("PLINK execution failed. Please check your installation and logs.")
    }
  } else {
    message("Skipping PLINK binary conversion as requested.")
  }

  message("All done.")
  invisible(NULL)
}
