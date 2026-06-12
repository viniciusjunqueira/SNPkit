# Build the package
{
  unlink("NAMESPACE")
  devtools::document()

  unlink("vignettes/*.html", recursive = TRUE)
  unlink("vignettes/*.rds", recursive = TRUE)
  devtools::clean_vignettes()
  devtools::build_vignettes()

  devtools::document()
  devtools::clean_dll()
  devtools::install(upgrade = "never")
  pkgdown::build_site()

  pkg_dir <- getwd()

  desc <- read.dcf(file.path(pkg_dir, "DESCRIPTION"))
  versao <- desc[1, "Version"]
  pdf_out <- file.path(pkg_dir, paste0("SNPkit_", versao, ".pdf"))

  if (file.exists(pdf_out)) file.remove(pdf_out)

  cmd <- paste("R CMD Rd2pdf", shQuote(pkg_dir), "-o", shQuote(pdf_out))
  system(cmd)

  if (file.exists(pdf_out)) {
    message("PDF manual successfully generated: ", pdf_out)
    # Uncomment to open automatically:
    browseURL(pdf_out)
  } else {
    warning("PDF manual generation failed!")
  }

  system(paste("R CMD build", shQuote(pkg_dir)))
}

#
# Check if tar.gz is fine as this is what is needed for CRAN
{
  system("R CMD check SNPkit_0.1.0.tar.gz")
}

#
# Check if all good
{
  # unlink("NAMESPACE")
  unlink("man/*.Rd")
  devtools::document()
  # devtools::build()
  # devtools::check()
  devtools::check(cran = TRUE)
  # devtools::check("SNPkit_0.1.0.tar.gz")
  # system("R CMD check SNPkit_0.1.0.tar.gz")
}


{
  tools::showNonASCIIfile("R/fimpute_runner.R")
  # Get all R source files in the R/ folder
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)

  # Check each file for non-ASCII characters
  for (f in r_files) {
    message("Checking: ", f)
    tools::showNonASCIIfile(f)
  }
}
