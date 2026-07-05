# SNPkit 0.1.1

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
