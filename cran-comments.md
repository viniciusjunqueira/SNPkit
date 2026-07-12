## Submission

This is a minor update (0.1.0 -> 0.1.1) with bug fixes and one new exported
function. Summary of changes (see NEWS.md for details):

* New `runPCA()`: run the genotype PCA without clustering; `runAnticlusteringPCA()`
  now reuses it internally.
* `doPCA()` is deprecated in favour of `runPCA()` (still works, emits a message).
* `getGeno()` bug fixes: builds the genotype matrix directly from
  `data.table::fread` (reliable for very large long-format FinalReport files),
  handles empty confidence fields, and normalises the map chromosome column
  name so panels can be combined.
* Fixes in `combineSNPData()`, `qcSNPs()`, and the `SNPDataLong` validity check.
* Performance: faster, lower-memory PCA for wide data; optional `RSpectra`
  (Suggests) enables a matrix-free truncated PCA.

## Test environments

* local macOS, R release
* GitHub Actions: ubuntu-latest (release, oldrel-1, devel), windows-latest
  (release), macos-latest (release)

## R CMD check results

0 errors | 0 warnings | 0 notes on the CI Linux/Windows checks.

On macOS with the current Apple clang / SDK, `R CMD check` reports one
compiler note from R's own header
(`R_ext/Boolean.h`: `-Wfixed-enum-extension`). It originates in R itself, not
in package code, and does not appear on the Linux/Windows checks. The package
compiled and passed the previous submission (0.1.0) unchanged in this respect.

## Reverse dependencies

There are no reverse dependencies.
