# run_admixture() renames each run's <prefix>.<K>.Q/.P to <out_prefix>.* so that
# successive same-K runs don't overwrite each other. The rename must not depend
# on the working directory (a relative output path used to double-nest and fail
# silently). These lock the fixed .rename_admixture_outputs() helper.

test_that("renames Q/P/log to the out_prefix inside an absolute dir", {
  d <- normalizePath(tempfile("admx"), mustWork = FALSE)
  dir.create(d)
  writeLines("q", file.path(d, "plink_data.3.Q"))
  writeLines("p", file.path(d, "plink_data.3.P"))
  writeLines("log", file.path(d, "plink_data.log"))

  ok <- suppressMessages(
    SNPkit:::.rename_admixture_outputs(d, "plink_data", 3, "unsupervised_run1_k3"))

  expect_true(ok)
  expect_true(file.exists(file.path(d, "unsupervised_run1_k3.Q")))
  expect_true(file.exists(file.path(d, "unsupervised_run1_k3.P")))
  expect_true(file.exists(file.path(d, "unsupervised_run1_k3.log")))
  expect_false(file.exists(file.path(d, "plink_data.3.Q")))  # moved, not copied
})

test_that("works when dir is '.' relative to the working directory (pre5 case)", {
  # Reproduces the failure: caller has setwd()'d into the output folder and the
  # rename target used to be re-resolved against it (double-nested).
  old <- getwd()
  d <- normalizePath(tempfile("admx"), mustWork = FALSE)
  dir.create(d)
  setwd(d)
  on.exit(setwd(old), add = TRUE)

  writeLines("q", "plink_data.2.Q")
  writeLines("p", "plink_data.2.P")

  ok <- suppressMessages(
    SNPkit:::.rename_admixture_outputs(".", "plink_data", 2, "supervised_run1_k2"))

  expect_true(ok)
  expect_true(file.exists("supervised_run1_k2.Q"))
  expect_false(file.exists("plink_data.2.Q"))
})

test_that("successive same-K runs do not overwrite once renamed", {
  d <- normalizePath(tempfile("admx"), mustWork = FALSE)
  dir.create(d)

  writeLines("unsup", file.path(d, "plink_data.3.Q"))
  writeLines("unsupP", file.path(d, "plink_data.3.P"))
  suppressMessages(
    SNPkit:::.rename_admixture_outputs(d, "plink_data", 3, "unsupervised_run1_k3"))

  # A second K=3 run writes plink_data.3.Q again; it must not clobber the first.
  writeLines("sup", file.path(d, "plink_data.3.Q"))
  writeLines("supP", file.path(d, "plink_data.3.P"))
  suppressMessages(
    SNPkit:::.rename_admixture_outputs(d, "plink_data", 3, "supervised_run1_k3"))

  expect_equal(readLines(file.path(d, "unsupervised_run1_k3.Q")), "unsup")
  expect_equal(readLines(file.path(d, "supervised_run1_k3.Q")), "sup")
})

test_that("warns and returns FALSE when the .Q output is missing", {
  d <- normalizePath(tempfile("admx"), mustWork = FALSE)
  dir.create(d)
  expect_warning(
    ok <- SNPkit:::.rename_admixture_outputs(d, "plink_data", 3, "x"),
    "not found")
  expect_false(ok)
})
