if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Name"))
}

#' Flexible and efficient genotype file reading using fread
#'
#' Imports SNP genotype data from Illumina FinalReport files using
#' \code{data.table::fread} and builds the \code{SnpMatrix} directly from the
#' long-format calls. This is reliable even for very large files (millions of
#' lines, hundreds of samples), where \code{snpStats::read.snps.long} may fail
#' to read all samples. Empty or unreadable confidence values are treated as no
#' calls. The original file on disk is never modified.
#'
#' @param path Path to the directory containing \code{FinalReport.txt}
#' @param fields List specifying column indices (sample, snp, allele1, allele2, confidence)
#' @param codes Allele codes (e.g., \code{c("A", "B")}); a genotype is coded as
#'   the count of \code{codes[2]} alleles (homozygous \code{codes[1]},
#'   heterozygous, homozygous \code{codes[2]}).
#' @param threshold Confidence threshold; calls below it are set to missing
#' @param sep Field separator
#' @param skip Lines to skip
#' @param verbose Logical; show progress
#' @param every Deprecated; kept for backward compatibility and ignored.
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

            # Read the required columns with fread and build the SnpMatrix
            # directly. We deliberately do NOT use snpStats::read.snps.long: its
            # internal search fails to scale to very large long-format files
            # (millions of lines / hundreds of samples), silently reading only
            # the first sample and reporting the rest as "not found". fread reads
            # the file reliably regardless of size, and empty/unreadable
            # confidence fields become NA here and are simply treated as no
            # calls -- no temporary file or line repair is needed.
            has_conf <- !is.null(fields$confidence) && fields$confidence > 0
            req <- c(sample  = fields$sample,
                     snp     = fields$snp,
                     allele1 = fields$allele1,
                     allele2 = fields$allele2)
            if (has_conf) req <- c(req, confidence = fields$confidence)
            sel <- as.integer(req)

            # Look up each required column's ORIGINAL header name from the
            # column-header line (the line right after `skip`). We reference
            # fread's output by these names rather than by position, because
            # different data.table versions return `select`ed columns either in
            # file order or in `select` order -- referencing by the preserved
            # header name is correct regardless.
            header_line <- readLines(orig_file, n = skip + 1L)
            col_names   <- strsplit(header_line[length(header_line)],
                                    sep, fixed = TRUE)[[1]]
            role_name <- vapply(req, function(i) col_names[i], character(1))
            if (anyNA(role_name)) {
              warning("Field index out of range for FinalReport columns at: ", path)
              return(NULL)
            }

            dt <- tryCatch({
              data.table::fread(
                file = orig_file,
                sep = sep,
                skip = skip,
                header = TRUE,
                select = sel,
                fill = TRUE,
                data.table = FALSE,
                colClasses = "character"
              )
            }, error = function(e) {
              warning("Failed to read FinalReport.txt with fread: ", e$message)
              return(NULL)
            })

            if (is.null(dt)) return(NULL)

            # Reference columns by their original header name (robust to the
            # order in which fread returns selected columns).
            dt_sample <- dt[[role_name[["sample"]]]]
            dt_snp    <- dt[[role_name[["snp"]]]]
            dt_a1     <- dt[[role_name[["allele1"]]]]
            dt_a2     <- dt[[role_name[["allele2"]]]]
            dt_conf   <- if (has_conf) dt[[role_name[["confidence"]]]] else NULL

            sample.id <- unique(dt_sample)
            snp.id    <- unique(dt_snp)

            if (length(sample.id) == 0) {
              warning("Sample IDs could not be determined correctly; setting default numeric row names.")
            }
            if (length(snp.id) == 0) {
              warning("SNP IDs could not be determined correctly; setting default numeric column names.")
            }

            # Encode each call into the SnpMatrix byte scheme:
            #   0x00 = missing/no call, 0x01 = codes[1]/codes[1],
            #   0x02 = heterozygous, 0x03 = codes[2]/codes[2].
            # A call is missing when either allele is a no call ("-") or, when a
            # confidence field is present, when it is below the threshold (this
            # also captures empty/unreadable confidence values, which parse to
            # NA).
            a1 <- dt_a1
            a2 <- dt_a2
            miss <- is.na(a1) | is.na(a2) | a1 == "-" | a2 == "-"
            if (has_conf) {
              conf_num <- suppressWarnings(as.numeric(dt_conf))
              miss <- miss | is.na(conf_num) | conf_num < threshold
            }
            n2   <- (a1 == codes[2]) + (a2 == codes[2])   # count of second allele
            code <- raw(nrow(dt))                         # all missing (0x00)
            keep <- !miss
            code[keep] <- as.raw(n2[keep] + 1L)

            ri <- match(dt_sample, sample.id)
            ci <- match(dt_snp,    snp.id)
            geno_mat <- matrix(as.raw(0),
                               nrow = length(sample.id),
                               ncol = length(snp.id),
                               dimnames = list(sample.id, snp.id))
            geno_mat[cbind(ri, ci)] <- code
            data <- new("SnpMatrix", geno_mat)

            if (verbose) {
              message("Built SnpMatrix with ", length(sample.id),
                      " samples x ", length(snp.id), " SNPs from ",
                      format(nrow(dt), big.mark = ","), " calls.")
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

            # Select columns and normalise the chromosome column name to
            # "Chromosome" so maps from panels using "Chr" and panels using
            # "Chromosome" can be row-bound together in combineSNPData().
            map <- map %>%
              dplyr::select(Name, !!chr_name, Position)
            names(map)[names(map) == chr_name] <- "Chromosome"

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
