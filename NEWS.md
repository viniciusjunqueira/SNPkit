# SNPkit 0.1.1

## Bug fixes

* `getGeno()`: lines in `FinalReport.txt` whose confidence (GC Score) field was
  exported empty are now handled gracefully. Such empty fields make
  `snpStats::read.snps.long` abort with "Failure to read confidence score", and
  once it fails it leaves the `snpStats` C state corrupted so that any retry in
  the same session silently skips most of the data. `getGeno()` now detects
  these rows up front (during the initial `fread` scan) and writes a temporary
  copy in which only the affected confidence fields are set to `0` -- the lowest
  confidence, so they are rejected by the threshold and treated as no calls,
  exactly like the well-formed rows that already carry a GC Score of `0`. The
  original file on disk is never modified, `read.snps.long` is called only once
  on valid data, and no genotypes are lost. Literal `NaN`/`Inf` confidence
  values are left untouched, since `read.snps.long` already tolerates them.

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
