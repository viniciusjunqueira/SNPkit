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

            orig_file <- file.path(path, "FinalReport.txt")

            # Fast initial read using fread. Include the confidence column so we
            # can detect malformed lines BEFORE the (single) read.snps.long call:
            # once read.snps.long fails on a bad line it leaves the snpStats C
            # state corrupted, so a second call in the same session silently
            # skips most of the data. We therefore never let read.snps.long fail
            # -- we repair the confidence field into a temporary file if needed.
            # The original file on disk is never modified.
            has_conf <- !is.null(fields$confidence) && fields$confidence > 0
            cols_to_read <- unique(unlist(
              fields[c("sample", "snp", if (has_conf) "confidence")]
            ))
            col_select <- as.integer(cols_to_read)

            fread_result <- tryCatch({
              data.table::fread(
                file = orig_file,
                sep = sep,
                skip = skip,
                header = TRUE,
                select = col_select,
                fill = TRUE,
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

            # Detect data rows whose confidence field cannot be parsed by
            # read.snps.long. An empty field (two consecutive separators) makes
            # read.snps.long abort with "Failure to read confidence score".
            # Note that literal "NaN"/"Inf" ARE tolerated by read.snps.long
            # (handled as no-calls), so we deliberately do NOT flag those, to
            # avoid needlessly rewriting files that would read fine as-is.
            read_file <- orig_file
            tmp_file  <- NULL
            if (has_conf) {
              conf.col  <- match(fields$confidence, col_select)
              conf_chr  <- trimws(fread_result[[conf.col]])
              conf_num  <- suppressWarnings(as.numeric(conf_chr))
              tolerated <- conf_chr %in% c("NaN", "nan", "NAN",
                                           "Inf", "-Inf", "inf", "-inf")
              bad_rows  <- which(is.na(conf_num) & !tolerated)
              if (length(bad_rows) > 0L) {
                message("Repairing ", length(bad_rows),
                        " line(s) with an empty/unreadable confidence score ",
                        "(set to 0 = no call) before reading genotypes.")
                # Data row i is at file line (skip + 1 + i): the first `skip`
                # lines plus the column-header line precede the data.
                bad_lines <- bad_rows + skip + 1L
                tmp_file  <- tempfile(fileext = ".txt")
                # Read the raw text and, ONLY on the affected lines, set the
                # confidence field to "0" (lowest confidence -> rejected by the
                # threshold -> treated as a no call, exactly like the well-formed
                # no-call rows that already carry a GC Score of 0). Every other
                # line is written verbatim, preserving the original formatting
                # and the sample/SNP identifiers. A single blocking readLines
                # avoids the line-splitting a chunked non-blocking read can cause.
                cf <- fields$confidence
                ok <- tryCatch({
                  all_lines <- readLines(orig_file, warn = FALSE)
                  bad_lines <- bad_lines[bad_lines <= length(all_lines)]
                  fixed <- vapply(all_lines[bad_lines], function(ln) {
                    v <- strsplit(ln, sep, fixed = TRUE)[[1]]
                    if (length(v) >= cf) v[cf] <- "0"
                    paste(v, collapse = sep)
                  }, character(1), USE.NAMES = FALSE)
                  all_lines[bad_lines] <- fixed
                  writeLines(all_lines, tmp_file)
                  file.exists(tmp_file) && file.info(tmp_file)$size > 0L
                }, error = function(ec) {
                  warning("Error writing preprocessed temp file: ", ec$message)
                  FALSE
                })
                if (ok) {
                  read_file <- tmp_file
                } else {
                  warning("Preprocessing failed; attempting to read the original file.")
                  unlink(tmp_file); tmp_file <- NULL
                }
              }
            }
            if (!is.null(tmp_file)) on.exit(unlink(tmp_file), add = TRUE)

            # Full genotype matrix (single read.snps.long call on a clean file)
            data <- tryCatch({
              snpStats::read.snps.long(
                file      = read_file,
                sample.id = sample.id,
                snp.id    = snp.id,
                fields    = fields,
                codes     = codes,
                threshold = threshold,
                sep       = sep,
                skip      = skip,
                verbose   = verbose,
                every     = every
              )
            }, error = function(e) {
              warning("Error while running read.snps.long: ", e$message)
              NULL
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
