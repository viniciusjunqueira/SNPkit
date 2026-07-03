# SNPkit 0.1.1

## Bug fixes

* `getGeno()`: when `snpStats::read.snps.long` fails because the input file
  contains malformed lines (an empty/unreadable confidence field, or a
  structurally incomplete line with fewer fields than expected), the function
  now re-reads the file with `data.table::fread(fill = TRUE)`, drops every row
  whose confidence value is missing or non-numeric, writes a clean temporary
  file, and retries `read.snps.long` on it. This recovers all valid genotypes
  instead of returning `NULL` and skipping the whole dataset.

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
