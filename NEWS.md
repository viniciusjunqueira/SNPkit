# SNPkit 0.1.1

## Bug fixes

* `getGeno()`: when `snpStats::read.snps.long` fails with a confidence score
  reading error (e.g. malformed lines in the input file), the function now
  automatically retries without confidence filtering instead of returning
  `NULL`. A warning is emitted to inform the user that the threshold was
  not applied.

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
