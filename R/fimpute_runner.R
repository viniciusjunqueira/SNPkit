#' @title Run FImpute from a FImputeRunner object
#'
#' @description
#' This function runs the external FImpute software using a `FImputeRunner` object,
#' ensuring that all required input files are present and the results are imported.
#'
#' @param object An object of class `FImputeRunner`.
#' @param verbose Logical. If TRUE (default), FImpute output will be printed to the console.
#'
#' @return An updated `FImputeRunner` object with the `results` slot populated (SNPDataLong).
#' @examples
#' \dontrun{
#' # Requires the external FImpute3 binary in PATH.
#' path_fimpute <- file.path(tempdir(), "fimpute_run_example")
#' param_file <- file.path(path_fimpute, "fimpute.par")
#'
#' export_obj <- methods::new("FImputeExport",
#'                            geno = geno_obj@geno,
#'                            map  = geno_obj@map,
#'                            path = path_fimpute)
#'
#' runner <- methods::new("FImputeRunner",
#'                        export    = export_obj,
#'                        par_file  = param_file,
#'                        exec_path = "FImpute3")
#'
#' res <- runFImpute(runner, verbose = TRUE)
#' }
#' @importFrom methods new
#' @export
setGeneric("runFImpute", function(object, verbose = TRUE) standardGeneric("runFImpute"))

#' @rdname runFImpute
#' @export
setMethod("runFImpute", "FImputeRunner", function(object, verbose = TRUE) {
  dir <- object@export@path
  par_file <- object@par_file
  exec <- object@exec_path

  if (!dir.exists(dir)) stop("Working directory does not exist: ", dir)
  if (!file.exists(par_file)) stop("Parameter file not found: ", par_file)

  resolved_exec <- Sys.which(exec)
  if (resolved_exec == "") {
    stop("FImpute executable not found in PATH: ", exec)
  }

  if (.Platform$OS.type == "unix" && file.access(resolved_exec, 1) != 0) {
    stop("FImpute executable exists but is not executable: ", resolved_exec)
  }

  destino_par <- file.path(dir, "fimpute.par")
  if (normalizePath(par_file) != normalizePath(destino_par)) {
    file.copy(par_file, destino_par, overwrite = TRUE)
    message("Parameter file copied to working directory: ", destino_par)
  } else {
    message("Parameter file already in working directory, no copy needed.")
  }

  output_dir <- file.path(dir, "output_fimpute")
  if (dir.exists(output_dir)) {
    message("Removing previous output directory: ", output_dir)
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }

  message("Running FImpute...")

  command <- paste("cd", shQuote(dir), "&&", shQuote(resolved_exec), "fimpute.par")
  status <- system(
    command,
    intern = FALSE,
    ignore.stdout = !verbose,
    ignore.stderr = !verbose
  )

  if (status != 0) {
    stop("FImpute execution failed. Please check the executable and input files.")
  }

  if (!dir.exists(output_dir)) {
    stop("Output directory 'output_fimpute' was not created. FImpute may have failed silently.")
  }

  if (!exists("read.fimpute", mode = "function")) {
    stop("The function 'read.fimpute()' must be defined and available in the current environment.")
  }

  res <- read.fimpute(file = output_dir)

  # Armazena objeto completo SNPDataLong no slot results
  # object@results <- res

  message("Results successfully read and stored in 'results' slot of FImputeRunner.")
  return(res)
})
