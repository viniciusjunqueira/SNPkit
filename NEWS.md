# SNPkit 0.1.1

## Bug fixes

* `getGeno()`: when `snpStats::read.snps.long` fails with a confidence score
  reading error (e.g. a malformed line with an empty confidence field), the
  function now removes the offending lines via `fread`, writes a clean
  temporary file, and retries `read.snps.long` on that file. The previous
  retry with `confidence = 0` was unreliable due to persistent internal state
  in `snpStats` after a failed call, causing most genotypes to be silently
  lost.

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
