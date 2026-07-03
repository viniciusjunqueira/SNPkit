if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Name"))
}

#' Flexible and efficient genotype file reading with autodetection using fread
#'
#' Allows flexible import of SNP genotype data from Illumina FinalReport files,
#' using fast initial column detection via \code{data.table::fread}, followed by
#' full genotype matrix construction with \code{snpStats::read.snps.long}.
#'
#' @param path Path to the directory containing \code{FinalReport.txt}
#' @param fields List specifying column indices (sample, snp, allele1, allele2, confidence)
#' @param codes Allele codes (e.g., \code{c("A", "B")})
#' @param threshold Confidence threshold
#' @param sep Field separator
#' @param skip Lines to skip
#' @param verbose Logical; show progress
#' @param every Frequency for progress updates
#' @param ... Additional optional arguments.
#'
#' @return An \code{SNPDataLong} object
#'
#' @importFrom utils read.table
#' @importFrom dplyr %>% select
#' @export
setGeneric("getGeno", function(...) standardGeneric("getGeno"))

#' @rdname getGeno
#' @export
setMethod("getGeno", signature(),
          function(path,
                   fields = list(sample = 2, snp = 1, allele1 = 7, allele2 = 8, confidence = 9),
                   codes = c("A", "B"),
                   threshold = 0.15,
                   sep = "\t",
                   skip = 0,
                   verbose = TRUE,
                   every = NULL) {

            if (!requireNamespace("data.table", quietly = TRUE)) {
              stop("The 'data.table' package is required. Install it using install.packages('data.table').")
            }

            if (!file.exists(file.path(path, "FinalReport.txt"))) {
              warning("File FinalReport.txt not found at: ", path)
              return(NULL)
            }

            # Fast initial read using fread
            cols_to_read <- unique(unlist(fields[c("sample", "snp")]))
            col_select <- as.integer(cols_to_read)

            fread_result <- tryCatch({
              data.table::fread(
                file = file.path(path, "FinalReport.txt"),
                sep = sep,
                skip = skip,
                header = TRUE,
                select = col_select,
                data.table = FALSE,
                colClasses = "character"
              )
            }, error = function(e) {
              warning("Failed to read FinalReport.txt with fread: ", e$message)
              return(NULL)
            })

            if (is.null(fread_result)) return(NULL)

            sample.col <- match(fields$sample, col_select)
            snp.col    <- match(fields$snp, col_select)

            sample.id <- unique(fread_result[[sample.col]])
            snp.id    <- unique(fread_result[[snp.col]])

            if (length(sample.id) == 0) {
              warning("Sample IDs could not be determined correctly; setting default numeric row names.")
            }
            if (length(snp.id) == 0) {
              warning("SNP IDs could not be determined correctly; setting default numeric column names.")
            }

            if (is.null(every)) every <- length(snp.id)

            # Full genotype matrix
            read_snps_long <- function(flds) {
              snpStats::read.snps.long(
                file      = file.path(path, "FinalReport.txt"),
                sample.id = sample.id,
                snp.id    = snp.id,
                fields    = flds,
                codes     = codes,
                threshold = threshold,
                sep       = sep,
                skip      = skip,
                verbose   = verbose,
                every     = every
              )
            }

            data <- tryCatch({
              read_snps_long(fields)
            }, error = function(e) {
              if (grepl("confidence", e$message, ignore.case = TRUE) &&
                  !is.null(fields$confidence) && fields$confidence > 0) {
                warning(
                  "Confidence score reading failed (", e$message, "). ",
                  "Removing malformed lines and retrying on a clean temporary file."
                )
                orig_file <- file.path(path, "FinalReport.txt")
                tmp_file  <- tempfile(fileext = ".txt")
                con_in    <- NULL
                con_out   <- NULL
                on.exit({
                  if (!is.null(con_in))  try(close(con_in),  silent = TRUE)
                  if (!is.null(con_out)) try(close(con_out), silent = TRUE)
                  unlink(tmp_file)
                }, add = TRUE)

                # --- Strategy 1: awk (fast, no RAM overhead) ---
                cleaned <- FALSE
                if (nchar(Sys.which("awk")) > 0) {
                  # Keep header lines and data lines where confidence field
                  # starts with a digit (valid numeric GC Score).
                  awk_prog  <- sprintf(
                    "NR<=%d || $%d~/^[0-9]/",
                    skip + 1L, fields$confidence
                  )
                  exit_code <- system2(
                    "awk",
                    args   = c("-F\t", awk_prog, orig_file),
                    stdout = tmp_file,
                    stderr = FALSE
                  )
                  cleaned <- (exit_code == 0 &&
                                file.exists(tmp_file) &&
                                file.info(tmp_file)$size > 0)
                }

                # --- Strategy 2: readLines chunk fallback ---
                if (!cleaned) {
                  conf_df <- tryCatch(
                    data.table::fread(
                      orig_file, sep = sep, skip = skip, header = TRUE,
                      select = as.integer(fields$confidence),
                      data.table = FALSE, colClasses = "character"
                    ),
                    error = function(ef) {
                      warning("fread failed reading confidence column: ", ef$message)
                      NULL
                    }
                  )
                  if (is.null(conf_df)) {
                    warning("Cannot preprocess file; skipping this path.")
                    return(NULL)
                  }
                  bad_rows  <- which(
                    suppressWarnings(is.na(as.numeric(conf_df[[1]])))
                  )
                  bad_lines <- bad_rows + skip + 1L
                  tryCatch({
                    con_in  <- file(orig_file, "r", blocking = FALSE)
                    con_out <- file(tmp_file, "w")
                    lnum <- 0L
                    repeat {
                      ch <- readLines(con_in, n = 50000L, warn = FALSE)
                      if (!length(ch)) break
                      idx  <- seq_along(ch) + lnum
                      keep <- !(idx %in% bad_lines)
                      if (any(keep)) writeLines(ch[keep], con_out)
                      lnum <- lnum + length(ch)
                    }
                    close(con_in);  con_in  <- NULL
                    close(con_out); con_out <- NULL
                    cleaned <- file.exists(tmp_file) && file.info(tmp_file)$size > 0
                  }, error = function(ec) {
                    warning("Error during file preprocessing: ", ec$message)
                  })
                }

                if (!cleaned) {
                  warning("Failed to create preprocessed file; skipping this path.")
                  return(NULL)
                }

                tryCatch(
                  snpStats::read.snps.long(
                    file      = tmp_file,
                    sample.id = sample.id,
                    snp.id    = snp.id,
                    fields    = fields,
                    codes     = codes,
                    threshold = threshold,
                    sep       = sep,
                    skip      = skip,
                    verbose   = verbose,
                    every     = every
                  ),
                  error = function(e2) {
                    warning("Error while running read.snps.long: ", e2$message)
                    NULL
                  }
                )
              } else {
                warning("Error while running read.snps.long: ", e$message)
                NULL
              }
            })

            if (is.null(data)) return(NULL)

            # Force row and column names
            if (is.null(rownames(data))) {
              rownames(data) <- sample.id
            }
            if (is.null(colnames(data))) {
              colnames(data) <- snp.id
            }

            # Read SNP map
            map_file <- file.path(path, "SNP_Map.txt")
            if (!file.exists(map_file)) {
              warning("SNP_Map.txt file not found at: ", path)
              return(NULL)
            }

            map <- tryCatch({
              utils::read.table(map_file, colClasses = "character", sep = sep, header = TRUE)
            }, error = function(e) {
              warning("Error reading SNP_Map.txt: ", e$message)
              return(NULL)
            })

            if (is.null(map)) return(NULL)

            # Chromosome column name
            possible_chr_cols <- c("Chromosome", "Chr")
            chr_name <- intersect(possible_chr_cols, colnames(map))

            if (length(chr_name) == 0) {
              stop("No chromosome column found (expected 'Chromosome' or 'Chr') on map file.")
            } else {
              chr_name <- chr_name[1]
            }

            # Select columns
            map <- map %>%
              dplyr::select(Name, !!chr_name, Position)

            new("SNPDataLong",
                geno = data,
                map  = map,
                path = path)
          }
)

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("combineSNPData"))
}

#' Import and combine multiple genotype configurations
#'
#' Imports genotype data from multiple configurations defined in an
#' \code{SNPImportList} object and combines them into a unified \code{SNPDataLong} object.
#'
#' @param object An \code{SNPImportList} object.
#'
#' @return A combined \code{SNPDataLong} object.
#'
#' @export
setGeneric("importAllGenos", function(object) standardGeneric("importAllGenos"))

#' @rdname importAllGenos
#' @export
setMethod("importAllGenos", "SNPImportList", function(object) {
  if (!inherits(object, "SNPImportList")) {
    stop("Input must be of class SNPImportList.")
  }

  all_genos <- lapply(object@configs, function(cfg) {
    tryCatch(
      getGeno(cfg),
      error = function(e) {
        warning("Error in getGeno(): ", e$message)
        NULL
      }
    )
  })

  all_genos <- Filter(Negate(is.null), all_genos)

  if (length(all_genos) == 0) {
    stop("No genotype data was successfully imported.")
  }

  combined <- combineSNPData(all_genos)

  return(combined)
})
