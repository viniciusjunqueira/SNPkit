# SNPkit 0.1.2

## Bug fixes

* `run_admixture()` renamed its outputs to a doubly-nested, non-existent path
  when `path` was relative (it `setwd()`s into `path`, then built the rename
  target with `file.path(path, ...)` again). The rename failed silently, so
  successive same-K runs overwrote each other's `.Q`/`.P`. `path` is now
  resolved to an absolute path and the rename is verified (warns on failure).
* `as_snpmatrix()` used a shifted byte encoding, so genotype `0` was stored as
  missing, missing as `BB`, and heterozygotes/alt-homozygotes were off by one.
  Any snpStats statistic (call rate, MAF, HWE, GRM/PCA) computed on its output
  was therefore wrong. It now uses the canonical `SnpMatrix` encoding
  (`0x00 = missing`, `0/1/2 -> 0x01/0x02/0x03`), matching `getGeno()`.
* `check.snp.same.position()` was defined twice; the shadowing copy lacked
  guards. It is now a single vectorized definition that groups SNPs by
  chromosome + position, safe for single-SNP chromosomes and missing positions.
* `check.identical.samples()`, `pairs2sets()` and
  `check.mendelian.inconsistencies()` no longer error on empty or single-element
  inputs (`1:n` / `2:n` off-by-one guards). `check.identical.samples()` also
  extracts pairs vectorially instead of growing a data.frame in a nested loop.
* `check.snp.no.position()` now treats a SNP as unmapped when its position is
  missing (`NA`), blank, non-numeric, or zero (previously only `== 0`, which
  also returned spurious `NA` names when positions were missing). `qcSNPs()`
  uses the same check consistently before the same-position filter.
* `check.identical.samples.by.block()` now returns the pairs identical across
  *every* block (data.frame), instead of only the last block's pairs, and is
  guarded against fewer than two samples.
* Removed the always-zero concordance term in `get.correl.fc()` (result
  unchanged).

## Removed

* `qcSNPs()` drops the `missing_ind` and `missing_snp` arguments, which were
  never implemented. Use `min_snp_cr` for per-SNP call rate and
  `qcSamples(smp_cr = ...)` for per-individual call rate.

## Internal changes

* `combineSNPData()` now assembles the combined genotype matrix through
  `rbindSnpFlexible()` instead of a duplicated inline implementation. The
  returned object is byte-identical to before (covered by a new equivalence
  test suite); only the intermediate "Adding N missing SNPs" message is no
  longer emitted.
* `rbindSnpFlexible()`, `rbind_SnpMatrix()` and `cbind_SnpMatrix()` are now
  internal helpers (no longer exported). Use `combineSNPData()` for the
  supported, high-level combining workflow.

## Tests

* Added a `testthat` suite, including a broad characterization battery that
  checks `combineSNPData()` against a reference copy of the pre-refactor
  implementation across overlapping/partial/disjoint/multi-object/randomized
  inputs and edge cases.


# SNPkit 0.1.1

## New features

* New exported function `runPCA()` runs the genotype PCA on a `SNPDataLong`
  object without any clustering, returning a `prcomp`-like object and the
  selected top principal components. It uses the same PCA engine as
  `runAnticlusteringPCA()` (standardised SNPs, Gram-matrix eigendecomposition,
  optional matrix-free `RSpectra` fast path), so scores are directly
  comparable. `runAnticlusteringPCA()` now calls `runPCA()` internally.

## Deprecations

* `doPCA()` is deprecated in favour of `runPCA()`. It still works but emits a
  deprecation message. `doPCA()` used an unscaled GRM (`snpStats::xxt`) on a raw
  `SnpMatrix`; `runPCA()` is the recommended, standardised PCA on a
  `SNPDataLong`.

## Performance

* `runAnticlusteringPCA()` is now dramatically faster and lighter on memory for
  wide genotype data. When there are more SNPs than individuals (the usual
  case), PCA is computed from the small n x n Gram matrix instead of a full SVD
  on the n x p matrix, avoiding the construction of the huge p x n rotation
  matrix (which could be several GB and was never used). The genotype matrix is
  also no longer converted to a `data.frame` before PCA. Scores and standard
  deviations are unchanged. `runAnticlusteringPCA()` no longer uses
  `anticlust`'s removed `features` argument.

* `runAnticlusteringPCA()`: when a fixed number of PCs is requested and the
  optional `RSpectra` package is installed, only the top `n_pcs` components are
  computed with a matrix-free solver, so the n x n Gram matrix is never formed
  at all. This is the dominant remaining cost for very wide data. Falls back to
  the Gram-matrix `eigen` decomposition when `RSpectra` is absent or a
  proportion of variance is requested (`n_pcs < 1`).

* `runAnticlusteringPCA()` gains an `anticlust_method` argument. The default
  `"exchange"` preserves current behaviour; `"fast"` uses
  `anticlust::fast_anticlustering`, which scales to large numbers of
  individuals.

## Bug fixes

* `getGeno()`: the chromosome column of the returned map is now always named
  `Chromosome` (previously it kept the original header name, which could be
  `Chr` or `Chromosome`). This fixes a `match.names ... names do not match`
  error in `combineSNPData()` / `import_geno_list()` when combining datasets
  whose `SNP_Map.txt` files used different chromosome column names.

* `getGeno()` now builds the `SnpMatrix` directly from the `data.table::fread`
  output instead of calling `snpStats::read.snps.long`. The latter's internal
  search does not scale to very large long-format `FinalReport.txt` files
  (millions of lines / hundreds of samples): it silently read only the first
  sample and reported all remaining rows as "not found", producing a matrix in
  which every sample but one was empty. Building the matrix from `fread` is
  reliable regardless of file size and is also robust to malformed lines --
  empty or unreadable confidence (GC Score) fields simply parse to `NA` and are
  treated as no calls, so no line repair or temporary file is required and the
  original file on disk is never modified.

* `combineSNPData()`: fixed spurious `"object has no names"` warning from
  `snpStats` when filling missing SNPs with NA. The `SnpMatrix` block is now
  constructed with `dimnames` set at creation time.

* `SNPDataLong`: relaxed validation of the `xref_path` slot to accept character
  vectors of any length (one entry per individual), resolving an error in
  `import_geno_list()` when datasets contained more than one individual.

* `qcSNPs()`: fixed a warning (`max` returning `-Inf`) and potential incorrect
  removal of all SNPs in a same-position group when all MAF values were `NA`.
  The first SNP in the group is now kept in that case.

# SNPkit 0.1.0

* Initial CRAN release.
