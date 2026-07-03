# SNPkit 0.1.1

## Bug fixes

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
