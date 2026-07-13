# as_snpmatrix() must produce the canonical snpStats byte encoding
# (0x00 = missing, 0x01/0x02/0x03 = 0/1/2 copies) so that round-tripping and
# downstream statistics are correct. Regression guard for the encoding bug.

dn <- function(nr, nc) list(paste0("s", seq_len(nr)), paste0("m", seq_len(nc)))

test_that("numeric 0/1/2/NA round-trips through as_snpmatrix unchanged", {
  g <- matrix(c(0, 1, 2, NA,
                2, 1, 0, NA),
              nrow = 2, byrow = TRUE, dimnames = dn(2, 4))
  sm <- as_snpmatrix(g)
  expect_equal(as(sm, "numeric"), g)
})

test_that("character coding '012' round-trips (missing preserved as NA)", {
  gc <- matrix(c("0", "1", "2", "NA",
                 "2", "1", "0", "."),
               nrow = 2, byrow = TRUE, dimnames = dn(2, 4))
  sm <- as_snpmatrix(gc, coding = "012", missing_codes = c("NA", "."))
  expect_equal(
    as(sm, "numeric"),
    matrix(c(0, 1, 2, NA, 2, 1, 0, NA), nrow = 2, byrow = TRUE, dimnames = dn(2, 4))
  )
})

test_that("character coding 'AAABBB' maps AA/AB/BB -> 0/1/2 and missing -> NA", {
  ga <- matrix(c("AA", "AB", "BB", ".",
                 "BB", "AB", "AA", "."),
               nrow = 2, byrow = TRUE, dimnames = dn(2, 4))
  sm <- as_snpmatrix(ga, coding = "AAABBB", missing_codes = ".")
  expect_equal(
    as(sm, "numeric"),
    matrix(c(0, 1, 2, NA, 2, 1, 0, NA), nrow = 2, byrow = TRUE, dimnames = dn(2, 4))
  )
})

test_that("raw bytes are canonical: 0/1/2 -> 1/2/3, missing -> 0", {
  g <- matrix(c(0L, 1L, 2L), nrow = 1, dimnames = dn(1, 3))
  expect_identical(as.integer(as_snpmatrix(g)@.Data), c(1L, 2L, 3L))

  g2 <- matrix(c(0L, NA, 2L), nrow = 1, dimnames = dn(1, 3))
  expect_identical(as.integer(as_snpmatrix(g2)@.Data), c(1L, 0L, 3L))

  # Same scheme getGeno uses: non-missing byte == genotype + 1.
  g3 <- matrix(c(0L, 1L, 2L, 1L), nrow = 2, dimnames = dn(2, 2))
  expect_identical(as_snpmatrix(g3)@.Data,
                   matrix(as.raw(g3 + 1L), nrow = 2, dimnames = dn(2, 2)))
})

test_that("col.summary sees the right calls, allele frequency and call rate", {
  # One SNP, B-allele counts 0,0,1,2 -> 3 B alleles of 8 -> RAF 0.375, all called.
  g <- matrix(c(0L, 0L, 1L, 2L), ncol = 1, dimnames = list(paste0("s", 1:4), "snp1"))
  cs <- snpStats::col.summary(as_snpmatrix(g))
  expect_equal(cs$Calls, 4)
  expect_equal(cs$Call.rate, 1)
  expect_equal(cs$RAF, 0.375)
})

test_that("input validation still rejects bad genotypes and ids", {
  expect_error(
    as_snpmatrix(matrix(c(0L, 1L, 5L, 2L), nrow = 2, dimnames = dn(2, 2))),
    "0/1/2"
  )
  expect_error(
    as_snpmatrix(matrix(0L, 2, 2, dimnames = list(c("a", "a"), c("x", "y")))),
    "Duplicate sample"
  )
  expect_error(as_snpmatrix(matrix(0L, 2, 2)), "rownames")
})
