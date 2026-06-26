#' Run ADMIXTURE analysis
#'
#' This function runs the ADMIXTURE program on a set of PLINK files (.bed/.bim/.fam)
#' located in a specified directory, using a given file prefix. It supports both unsupervised
#' and supervised analyses, optional cross-validation, and custom output file prefixes to avoid overwriting results.
#'
#' @param path Character. Path to the folder containing PLINK files.
#' @param prefix Character. File prefix (without extension). The function will look for `<prefix>.bed`, `<prefix>.bim`, and `<prefix>.fam` in `path`.
#' @param admixture_path Character. Path to the ADMIXTURE executable, or "admixture" if in system PATH. Default is "admixture".
#' @param K Integer. Number of ancestral populations to estimate.
#' @param supervised Logical. If TRUE, runs ADMIXTURE in supervised mode (requires \code{pop_assignments}). Default is FALSE.
#' @param pop_assignments Character vector. Population assignments for each individual (length equal to number of individuals in `.fam`). Use \code{NA} or "-" for missing. Required if \code{supervised = TRUE}.
#' @param extra_args Character vector. Additional arguments to pass to ADMIXTURE (e.g., other flags). Default is NULL.
#' @param out_prefix Character. Optional prefix for renaming output files (.Q, .P, .log) after the run completes. Default is NULL.
#' @param cv Integer. Number of folds for cross-validation (e.g., 5 or 10). If provided, adds \code{--cv=cv}. Default is NULL.
#'
#' @return No value returned. Runs ADMIXTURE as a side effect. Generates output files in the specified directory. Messages indicate progress and output file names.
#'
#' @details
#' When \code{supervised = TRUE}, a `.pop` file is automatically created in the specified directory.
#' Each line in this file corresponds to one individual, containing the population name or "-" for missing assignments.
#' 
#' If \code{out_prefix} is provided, the function renames the standard ADMIXTURE output files 
#' (e.g., `<prefix>.3.Q`) to use this prefix (e.g., `myrun.Q`).
#'
#' The function only works on Linux or macOS systems.
#'
#' @examples
#' \dontrun{
#' # Requires the external ADMIXTURE binary and PLINK files prepared beforehand.
#' work_dir <- file.path(tempdir(), "admixture_demo")
#' run_admixture(
#'   path = work_dir,
#'   prefix = "plink_data",
#'   admixture_path = "admixture",
#'   K = 3,
#'   out_prefix = "run1_k3"
#' )
#'
#' pop_vec <- c("A", "A", "B", "B", "-", "-", "A", "B", "A", "-")
#' run_admixture(
#'   path = work_dir,
#'   prefix = "plink_data",
#'   admixture_path = "admixture",
#'   K = 3,
#'   supervised = TRUE,
#'   pop_assignments = pop_vec,
#'   cv = 10,
#'   out_prefix = "supervised_k3_cv10"
#' )
#' }
#'
#' @export
run_admixture <- function(path, prefix, admixture_path = "admixture", K,
                          supervised = FALSE, pop_assignments = NULL,
                          extra_args = NULL, out_prefix = NULL, cv = NULL) {
  # Check OS
  sys <- Sys.info()[["sysname"]]
  if (!sys %in% c("Linux", "Darwin")) {
    stop("ADMIXTURE can only be run on Linux or macOS systems.")
  }
  
  # Check if admixture executable is available
  admix_exec <- Sys.which(admixture_path)
  if (admix_exec == "") {
    stop("ADMIXTURE executable not found. Please ensure it is installed and available in your PATH or provide full path.")
  }
  
  # Build file paths
  bed_file <- file.path(path, paste0(prefix, ".bed"))
  bim_file <- file.path(path, paste0(prefix, ".bim"))
  fam_file <- file.path(path, paste0(prefix, ".fam"))
  
  if (!file.exists(bed_file)) stop("BED file not found: ", bed_file)
  if (!file.exists(bim_file)) stop("BIM file not found: ", bim_file)
  if (!file.exists(fam_file)) stop("FAM file not found: ", fam_file)
  
  # Read FAM
  fam_data <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  n_individuals <- nrow(fam_data)
  
  if (supervised) {
    if (is.null(pop_assignments)) {
      stop("When supervised = TRUE, you must provide pop_assignments vector.")
    }
    if (length(pop_assignments) != n_individuals) {
      stop("Length of pop_assignments does not match number of individuals in .fam file.")
    }
    
    # Write .pop file
    pop_file <- file.path(path, paste0(prefix, ".pop"))
    pop_values <- ifelse(is.na(pop_assignments), "-", pop_assignments)
    writeLines(pop_values, con = pop_file)
    message(".pop file created at: ", pop_file)
  }
  
  # Change to analysis folder first
  oldwd <- getwd()
  setwd(path)
  on.exit(setwd(oldwd), add = TRUE)
  
  # Local file name only
  bed_file_name <- paste0(prefix, ".bed")
  cmd_args <- character()
  if (supervised) {
    cmd_args <- c(cmd_args, "--supervised")
  }
  if (!is.null(cv)) {
    if (!is.numeric(cv) || cv <= 1) {
      stop("cv must be a numeric value > 1 (e.g., 5 or 10).")
    }
    cmd_args <- c(cmd_args, paste0("--cv=", cv))
  }
  if (!is.null(extra_args)) {
    cmd_args <- c(cmd_args, extra_args)
  }
  cmd_args <- c(cmd_args, bed_file_name, as.character(K))
  
  # Run
  message("Running ADMIXTURE: ", admix_exec, " ", paste(cmd_args, collapse = " "))
  res <- system2(admix_exec, args = cmd_args)
  
  if (res != 0) {
    warning("ADMIXTURE did not finish successfully. Check logs.")
  } else {
    message("ADMIXTURE run completed successfully.")
    
    if (!is.null(out_prefix)) {
      q_file <- paste0(prefix, ".", K, ".Q")
      p_file <- paste0(prefix, ".", K, ".P")
      log_file <- paste0(prefix, ".log")
      
      new_q_file <- file.path(path, paste0(out_prefix, ".Q"))
      new_p_file <- file.path(path, paste0(out_prefix, ".P"))
      new_log_file <- file.path(path, paste0(out_prefix, ".log"))
      
      if (file.exists(q_file)) file.rename(q_file, new_q_file)
      if (file.exists(p_file)) file.rename(p_file, new_p_file)
      if (file.exists(log_file)) file.rename(log_file, new_log_file)
      
      message("Output files renamed with prefix: ", out_prefix)
    }
  }
  
  invisible(NULL)
}
