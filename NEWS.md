# SNPkit 0.1.1

## Bug fixes

* `getGeno()`: malformed lines in `FinalReport.txt` (an empty/unreadable
  confidence field, or a structurally incomplete line with fewer fields than
  expected) are now detected *before* reading genotypes. The initial `fread`
  scan (with `fill = TRUE`) flags rows whose confidence value is missing or
  non-numeric; if any are found, those lines are removed from the raw text and
  `read.snps.long` reads a clean temporary file instead. This matters because
  once `read.snps.long` fails on a bad line it leaves the `snpStats` C state
  corrupted, so any retry in the same session silently skips most of the data.
  By cleaning the file up front, `read.snps.long` is only ever called once, on
  valid data, and all genotypes are recovered.

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
