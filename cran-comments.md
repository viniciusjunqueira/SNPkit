## Submission

This is an update from the current CRAN version 0.1.0 to 0.1.2. The intervening
0.1.1 was prepared but not released on CRAN, so NEWS.md documents both the
0.1.1 and 0.1.2 sections. Highlights of the changes since 0.1.0 (see NEWS.md
for the full list):

New features and improvements

* New `runPCA()`: genotype PCA without clustering, reused internally by
  `runAnticlusteringPCA()`. For wide data it avoids forming the large rotation
  matrix and can use a matrix-free truncated solver when the optional
  'RSpectra' package (Suggests) is installed.
* `savePlink()` gains an `extra_args` argument to pass flags through to PLINK
  (e.g. `--chr-set`) for non-human data.
* `getGeno()` now builds the genotype matrix directly from
  `data.table::fread`, which is reliable for very large long-format
  FinalReport files; it also handles empty confidence fields and normalises
  the map chromosome column name so panels can be combined.

Bug fixes

* `as_snpmatrix()` used a shifted byte encoding (genotype 0 was stored as
  missing and missing as BB); it now uses the canonical `SnpMatrix` encoding.
* `run_admixture()` renamed its output files to a doubly-nested path when the
  output directory was relative, so the rename failed silently; the path is
  now resolved to an absolute path and the rename is verified.
* Correctness and input-guarding fixes in `combineSNPData()`, `qcSNPs()`,
  `check.snp.no.position()`, `check.snp.same.position()` (previously defined
  twice), `check.identical.samples()`, `pairs2sets()`,
  `check.mendelian.inconsistencies()`, and the `SNPDataLong` validity check.

Deprecations and interface changes

* `doPCA()` is deprecated in favour of `runPCA()` (still works, emits a
  message).
* The low-level helpers `rbind_SnpMatrix()`, `cbind_SnpMatrix()` and
  `rbindSnpFlexible()`, exported in 0.1.0, are no longer exported; they are
  internal implementation details of `combineSNPData()`. There are no reverse
  dependencies, so no packages are affected.
* `qcSNPs()` drops the never-implemented `missing_ind` and `missing_snp`
  arguments.

Testing

* Added a testthat test suite.

## Test environments

* local macOS, R release: `R CMD check --as-cran`
* GitHub Actions: ubuntu-latest (release, oldrel-1, devel), windows-latest
  (release), macos-latest (release)

## R CMD check results

0 errors | 0 warnings | 0 notes attributable to the package.

On the local macOS check the only messages are environmental and do not
originate in package code:

* one compiler warning from R's own header
  (`R_ext/Boolean.h`: `-Wfixed-enum-extension`), emitted by a very recent Apple
  clang and not reproducible on the Linux/Windows checks;
* a note "unable to verify current time" (no network access to the time
  server);
* a note that the local HTML Tidy is too old to validate the HTML manual.

## Reverse dependencies

There are no reverse dependencies.
