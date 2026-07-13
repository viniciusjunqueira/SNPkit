# combineSNPData must return exactly the same object as the pre-refactor
# implementation across a wide range of inputs.

test_that("fully overlapping SNPs (identical colnames)", {
  o1 <- make_obj(c("A1", "A2", "A3"), c("s1", "s2", "s3"), seed = 1)
  o2 <- make_obj(c("B1", "B2"), c("s1", "s2", "s3"), seed = 2)
  expect_snpdata_equal(combine_new(list(o1, o2)), combine_legacy(list(o1, o2)))
})

test_that("partially overlapping SNPs", {
  o1 <- make_obj(c("A1", "A2"), c("s1", "s2", "s3"), seed = 3)
  o2 <- make_obj(c("B1", "B2", "B3"), c("s2", "s3", "s4", "s5"), seed = 4)
  expect_snpdata_equal(combine_new(list(o1, o2)), combine_legacy(list(o1, o2)))
})

test_that("disjoint SNP sets", {
  o1 <- make_obj(c("A1", "A2"), c("s1", "s2"), seed = 5)
  o2 <- make_obj(c("B1", "B2"), c("s3", "s4"), seed = 6)
  expect_snpdata_equal(combine_new(list(o1, o2)), combine_legacy(list(o1, o2)))
})

test_that("five objects with mixed overlaps (pre5-like)", {
  objs <- list(
    make_obj(paste0("ang", 1:6),  c("s1", "s2", "s3", "s4"),       seed = 11),
    make_obj(paste0("bra", 1:3),  c("s2", "s3", "s5"),             seed = 12),
    make_obj(paste0("brg", 1:5),  c("s1", "s3", "s4", "s6"),       seed = 13),
    make_obj(paste0("nel", 1:2),  c("s4", "s5", "s6", "s7"),       seed = 14),
    make_obj(paste0("ult", 1:4),  c("s1", "s2", "s6", "s8"),       seed = 15)
  )
  expect_snpdata_equal(combine_new(objs), combine_legacy(objs))
})

test_that("SNPs in different column orders across objects", {
  o1 <- make_obj(c("A1", "A2"), c("s3", "s1", "s2"), seed = 21)
  o2 <- make_obj(c("B1", "B2"), c("s2", "s4", "s1"), seed = 22)
  expect_snpdata_equal(combine_new(list(o1, o2)), combine_legacy(list(o1, o2)))
})

test_that("single-object list does not error and matches legacy", {
  o1 <- make_obj(c("A1", "A2", "A3"), c("s1", "s2", "s3"), seed = 31)
  res_new <- combine_new(list(o1))
  res_leg <- combine_legacy(list(o1))
  expect_snpdata_equal(res_new, res_leg)
  # And its genotype/map are just the input's, unchanged.
  expect_snpmatrix_equal(res_new@geno, o1@geno)
  expect_identical(res_new@map, o1@map)
})

test_that("missing rownames trigger the same default-name fallback", {
  o1 <- make_obj(c("A1", "A2"), c("s1", "s2"), seed = 41)
  rownames(o1@geno) <- NULL
  o2 <- make_obj(c("B1", "B2"), c("s2", "s3"), seed = 42)

  expect_warning(
    res_new <- suppressMessages(combineSNPData(list(o1, o2))),
    "missing rownames"
  )
  res_leg <- suppressWarnings(combine_legacy(list(o1, o2)))
  expect_snpdata_equal(res_new, res_leg)
  # First two rows were renamed Sample_1 / Sample_2.
  expect_identical(rownames(res_new@geno)[1:2], c("Sample_1", "Sample_2"))
})

test_that("edge sizes: single SNP and single sample per object", {
  o1 <- make_obj(c("A1", "A2"), "s1", seed = 51)
  o2 <- make_obj("B1", c("s1", "s2"), seed = 52)
  expect_snpdata_equal(combine_new(list(o1, o2)), combine_legacy(list(o1, o2)))
})

test_that("duplicate SNP Name across objects dedupes identically", {
  # s2 shared between both objects but with different Position/Chromosome in the
  # map; combine keeps the first occurrence for both implementations.
  o1 <- make_obj(c("A1", "A2"), c("s1", "s2"), seed = 61, chr = 1L)
  o2 <- make_obj(c("B1", "B2"), c("s2", "s3"), seed = 62, chr = 2L)
  o2@map$Position <- c(999L, 1000L)
  expect_snpdata_equal(combine_new(list(o1, o2)), combine_legacy(list(o1, o2)))
})

test_that("path and xref_path are concatenated identically", {
  o1 <- make_obj(c("A1"), c("s1", "s2"), seed = 71, path = "pathA", xref = "chipA")
  o2 <- make_obj(c("B1"), c("s2", "s3"), seed = 72, path = "pathB", xref = "chipB")
  res <- combine_new(list(o1, o2))
  expect_identical(res@path, "pathA;pathB")
  expect_identical(res@xref_path, "chipA;chipB")
  expect_snpdata_equal(res, combine_legacy(list(o1, o2)))
})

test_that("randomized inputs: new == legacy across many seeds", {
  pool <- paste0("s", 1:12)
  for (seed in 1:25) {
    set.seed(seed)
    n_obj <- sample(2:5, 1)
    objs <- lapply(seq_len(n_obj), function(i) {
      n_samp <- sample(1:6, 1)
      n_snp  <- sample(2:length(pool), 1)
      snps   <- sample(pool, n_snp)                       # random subset & order
      samples <- paste0("o", i, "_s", seq_len(n_samp))
      make_obj(samples, snps, seed = seed * 100 + i)
    })
    expect_snpdata_equal(combine_new(objs), combine_legacy(objs))
  }
})

test_that("returned object passes full identical() (strong check)", {
  o1 <- make_obj(c("A1", "A2"), c("s1", "s2", "s3"), seed = 81)
  o2 <- make_obj(c("B1", "B2", "B3"), c("s2", "s3", "s4"), seed = 82)
  expect_identical(combine_new(list(o1, o2)), combine_legacy(list(o1, o2)))
})
