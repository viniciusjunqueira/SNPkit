# Robustness of the fQC helper checks: single/empty inputs must not error, and
# check.snp.same.position (previously duplicated) must group correctly.

test_that("check.snp.same.position groups SNPs by chromosome + position", {
  map <- data.frame(
    Name       = c("a", "b", "c", "d", "e", "f"),
    Chromosome = c(1, 1, 1, 2, 2, 3),
    Position   = c(100, 100, 200, 50, 50, 999),
    stringsAsFactors = FALSE)
  res <- check.snp.same.position(map)
  expect_length(res, 2)                               # {a,b} and {d,e}
  expect_true(any(vapply(res, setequal, logical(1), c("a", "b"))))
  expect_true(any(vapply(res, setequal, logical(1), c("d", "e"))))
})

test_that("same position on different chromosomes is not grouped", {
  map <- data.frame(Name = c("a", "b"), Chromosome = c(1, 2), Position = c(100, 100))
  expect_length(check.snp.same.position(map), 0)
})

test_that("single-SNP chromosome and NA positions are safe (no error)", {
  map <- data.frame(Name = c("a", "b", "c"),
                    Chromosome = c(1, 2, 2),
                    Position = c(100, NA, 300))
  expect_length(check.snp.same.position(map), 0)
})

test_that("three or more SNPs at one position form a single group", {
  map <- data.frame(Name = c("a", "b", "c"), Chromosome = 1, Position = 100)
  res <- check.snp.same.position(map)
  expect_length(res, 1)
  expect_setequal(res[[1]], c("a", "b", "c"))
})

test_that("check.identical.samples handles <2 samples and finds duplicates", {
  one <- matrix(c(0, 1, 2), nrow = 1, dimnames = list("s1", NULL))
  expect_equal(nrow(check.identical.samples(one)), 0)

  m <- rbind(s1 = c(0, 1, 2, 0), s2 = c(0, 1, 2, 0), s3 = c(2, 2, 0, 1))
  res <- suppressWarnings(check.identical.samples(m, threshold = 0))
  expect_equal(nrow(res), 1)
  expect_setequal(c(res$Sample1, res$Sample2), c("s1", "s2"))
})

test_that("check.identical.samples returns empty when nothing within threshold", {
  m <- rbind(s1 = c(0, 0, 0), s2 = c(2, 2, 2))
  expect_equal(nrow(suppressWarnings(check.identical.samples(m, threshold = 0))), 0)
})

test_that("pairs2sets merges chained pairs and handles a single pair", {
  single <- matrix(c("A", "B"), ncol = 2, byrow = TRUE)
  expect_equal(pairs2sets(single), list(c("A", "B")))

  chained <- matrix(c("A", "B", "B", "C", "D", "E"), ncol = 2, byrow = TRUE)
  sets <- pairs2sets(chained)
  expect_length(sets, 2)
  expect_true(any(vapply(sets, setequal, logical(1), c("A", "B", "C"))))
  expect_true(any(vapply(sets, setequal, logical(1), c("D", "E"))))
})

test_that("check.mendelian.inconsistencies: empty child is safe; counts a pair", {
  geno <- rbind(F1 = c(1, 1, 3, 3, 2), C1 = c(3, 1, 1, 3, 2))
  expect_equal(
    nrow(check.mendelian.inconsistencies(geno, father = "F1", child = character(0))),
    0
  )
  res <- suppressMessages(
    check.mendelian.inconsistencies(geno, father = "F1", child = "C1")
  )
  expect_equal(res$N, 2)      # positions 1 and 3 are 1<->3 inconsistencies
  expect_equal(res$Total, 4)  # 4 comparable homozygous positions
})
