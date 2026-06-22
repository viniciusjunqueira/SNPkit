#' @title Save genotype and map files in FImpute format
#'
#' @description
#' S4 method to export genotype (.gen), map (.map), and parameter (fimpute.par) files compatible with [FImpute](https://www.aps.uoguelph.ca/~msargol/fimpute/).
#'
#' @param object An object of class `FImputeExport` or `SNPDataLong`.
#' @param path Output directory. Must be supplied by the caller (e.g. a path
#'   inside \code{tempdir()} for examples).
#' @param ... Additional arguments passed to methods.
#'
#' @return No return value, called for side effects. The function writes the
#'   files \code{data.gen}, \code{data.map}, and \code{fimpute.par} to the
#'   directory \code{path} and returns \code{NULL} invisibly.
#'
#' @export
setGeneric("saveFImpute", function(object, ...) standardGeneric("saveFImpute"))

#' @rdname saveFImpute
#' @export
setMethod("saveFImpute", "FImputeExport", function(object) {
  save_fimpute_raw(object@geno, object@map, object@path)
})

#' @rdname saveFImpute
#' @export
setMethod("saveFImpute", "SNPDataLong", function(object, path) {
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("The 'snpStats' package is required. Install it with install.packages('snpStats').")
  }

  if (!inherits(object@geno, "SnpMatrix")) {
    stop("The 'geno' slot must be of class 'SnpMatrix'.")
  }

  if (!is.data.frame(object@map)) {
    stop("The 'map' slot must be a data.frame.")
  }

  if (missing(path) || !is.character(path) || length(path) != 1) {
    stop("'path' must be a single character string indicating the output directory.")
  }

  save_fimpute_raw(object@geno, object@map, path, xref = object@xref_path)
})

#' @title Export genotypes and map using basic arguments
#'
#' @description
#' Convenience function to export FImpute files directly from a `SnpMatrix` and map `data.frame`.
#'
#' @param geno A `SnpMatrix` object.
#' @param map A data.frame with columns 'Name', 'Chromosome', 'Position', and 'SourcePath'.
#' @param path Output directory.
#' @param xref Optional vector of identifiers per individual (used to assign numeric chip IDs).
#'
#' @return No return value, called for side effects. The function writes
#'   three files (\code{data.gen}, \code{data.map}, and \code{fimpute.par}) to
#'   the directory specified by \code{path} and returns \code{NULL} invisibly.
#'
#' @importFrom utils write.table
#' @export
saveFImputeRaw <- function(geno, map, path, xref = NULL) {
  save_fimpute_raw(geno, map, path, xref = xref)
}

#' Internal function: writes files in FImpute format (.gen, .map, fimpute.par)
#'
#' @noRd
save_fimpute_raw <- function(genotype, map, path, xref = NULL) {
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("The 'snpStats' package is required. Install it with install.packages('snpStats').")
  }

  if (!inherits(genotype, "SnpMatrix")) {
    stop("The 'genotype' object must be of class 'SnpMatrix'.")
  }

  if (!is.data.frame(map)) {
    stop("The 'map' argument must be a data.frame.")
  }

  qc_header("Saving FImpute Files")

  required_cols <- c("Name", "Chromosome", "Position", "SourcePath")
  missing_cols <- setdiff(required_cols, colnames(map))
  if (length(missing_cols) > 0) {
    stop("The map is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.character(path) || length(path) != 1) {
    stop("'path' must be a single character string indicating the output directory.")
  }

  if (!dir.exists(path)) {
    message("Creating output directory: ", path)
    dir.create(path, recursive = TRUE)
  } else {
    existing_files <- list.files(path, pattern = "data\\.(gen|map|par)$", full.names = TRUE)
    if (length(existing_files) > 0) {
      warning("The following files will be overwritten:\n  ",
              paste(basename(existing_files), collapse = "\n  "))
    }
  }

  snp_order <- colnames(genotype)

  map_final <- map[map$Name %in% snp_order, , drop = FALSE]
  map_final <- map_final[match(snp_order, map_final$Name), ]

  stopifnot(identical(map_final$Name, snp_order))

  chips <- unique(map_final$SourcePath)
  chip_names <- paste0("Chip", seq_along(chips))

  # if (!is.null(xref)) {
  #   if (!all(xref %in% chips)) {
  #     stop("The values in 'xref' must exactly match 'SourcePath' values in the map.")
  #   }
  #   xref_to_num <- setNames(seq_along(chips), chips)
  #   chip_ids <- as.integer(xref_to_num[xref])
  # } else {
  #   chip_ids <- rep(1L, nrow(genotype))
  # }

  # Temporary solutions
  chip_ids <- rep(1L, nrow(genotype))

  ## Write .gen file
  gen_file <- file.path(path, "data.gen")
  con <- file(gen_file, "wt")

  smp <- rownames(genotype)
  if (is.null(smp)) stop("The 'genotype' object must have row names (individual IDs).")

  nC <- max(nchar(smp), na.rm = TRUE)
  writeLines("ID    Chip                   Call...", con)

  pb <- utils::txtProgressBar(min = 0, max = nrow(genotype), style = 3)

  for (i in seq_len(nrow(genotype))) {
    temp <- as(genotype[i, ], "character")
    temp[temp == "A/A"] <- "0"
    temp[temp == "A/B"] <- "1"
    temp[temp == "B/B"] <- "2"
    temp[is.na(temp) | temp == "NA"] <- "5"
    linha <- paste(format(smp[i], width = nC, justify = "left"),
                   chip_ids[i],
                   paste(temp, collapse = ""))
    writeLines(linha, con)

    utils::setTxtProgressBar(pb, i)
  }

  close(pb)
  close(con)
  message("File written: ", gen_file)

  ## Write .map file
  map_file <- file.path(path, "data.map")

  # Chips
  # chips <- unique(map_final$SourcePath)
  # chip_names <- paste0("Chip", seq_along(chips))

  # header <- c("SNP_ID", "Chr", "Pos", chip_names)
  # header_line <- paste(header, collapse = " ")
  header <- c("SNP_ID", "Chr", "Pos", "Chip1")
  header_line <- paste(header, collapse = " ")

  # Escrever
  write(header_line, file = map_file)

  # for (chip in chip_names) {
  #   map_final[[chip]] <- 0
  # }

  # chip_cols <- chip_names[match(map_final$SourcePath, chips)]

  # for (chip in chip_names) {
  #   idx <- which(chip_cols == chip)
  #   if (length(idx) > 0) {
  #     map_final[idx, chip] <- seq_along(idx)
  #   }
  # }

  map_out <- data.frame(
    SNP_ID = map_final$Name,
    Chr = map_final$Chromosome,
    Pos = map_final$Position,
    # map_final[, chip_names, drop = FALSE],
    1:nrow(map_final),
    stringsAsFactors = FALSE
  )

  utils::write.table(
    map_out,
    file = map_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = " ",
    append = TRUE
  )
  message("File written: ", map_file)

  ## Write fimpute.par file
  par_file <- file.path(path, "fimpute.par")
  writeLines(c(
    'title="FImpute imputation";',
    'genotype_file="data.gen";',
    'snp_info_file="data.map";',
    'parentage_test /ert_mm=0.02 /remove_conflict /find_match_cnflt /find_match_mp /find_match_ugp;',
    'output_folder="output_fimpute";',
    'save_genotype;',
    'njob=24;'
  ), con = par_file)
  message("File written: ", par_file)

  invisible(NULL)
}
