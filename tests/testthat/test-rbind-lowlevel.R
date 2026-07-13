# Low-level (now internal) SnpMatrix binding primitives.

test_that("rbindSnpFlexible fills absent SNPs with missing (raw 0)", {
  m1 <- make_snpmatrix(c("A1", "A2"), c("s1", "s2"), seed = 1)
  m2 <- make_snpmatrix(c("B1", "B2"), c("s2", "s3"), seed = 2)
  res <- SNPkit:::rbindSnpFlexible(m1, m2)

  expect_identical(colnames(res), c("s1", "s2", "s3"))
  expect_identical(rownames(res), c("A1", "A2", "B1", "B2"))
  # s3 absent in m1 -> missing for its rows; s1 absent in m2 -> missing there.
  expect_true(all(res[c("A1", "A2"), "s3"]@.Data == as.raw(0)))
  expect_true(all(res[c("B1", "B2"), "s1"]@.Data == as.raw(0)))
  # Present entries are copied verbatim.
  expect_identical(res[c("A1", "A2"), "s1"]@.Data, m1[, "s1"]@.Data)
  expect_identical(res[c("B1", "B2"), "s3"]@.Data, m2[, "s3"]@.Data)
})

test_that("rbindSnpFlexible rejects fewer than two matrices", {
  m1 <- make_snpmatrix("A1", c("s1", "s2"), seed = 3)
  expect_error(SNPkit:::rbindSnpFlexible(m1), "two SnpMatrix")
})

test_that("rbindSnpFlexible rejects non-SnpMatrix input", {
  m1 <- make_snpmatrix(c("A1", "A2"), c("s1", "s2"), seed = 4)
  expect_error(
    SNPkit:::rbindSnpFlexible(m1, matrix(as.raw(0), 2, 2)),
    "must be SnpMatrix"
  )
})

test_that("union column ordering matches Reduce(union): unique(unlist()) == Reduce(union)", {
  pool <- paste0("s", 1:10)
  for (seed in 1:30) {
    set.seed(seed)
    cols <- lapply(1:sample(2:5, 1), function(i) sample(pool, sample(2:8, 1)))
    expect_identical(unique(unlist(cols)), Reduce(union, cols))
  }
})

test_that("engines agree when columns are identical (flexible == strict)", {
  m1 <- make_snpmatrix(c("A1", "A2"), c("s1", "s2", "s3"), seed = 5)
  m2 <- make_snpmatrix(c("B1", "B2"), c("s1", "s2", "s3"), seed = 6)
  flex   <- SNPkit:::rbindSnpFlexible(m1, m2)
  strict <- SNPkit:::rbind_SnpMatrix(m1, m2)
  expect_identical(dimnames(flex), dimnames(strict))
  expect_identical(snpmatrix_bytes(flex), snpmatrix_bytes(strict))
})

test_that("rbind_SnpMatrix is strict about column names and keeps dimnames", {
  m1 <- make_snpmatrix(c("A1", "A2"), c("s1", "s2"), seed = 7)
  m2 <- make_snpmatrix(c("B1", "B2"), c("s1", "s2"), seed = 8)
  res <- SNPkit:::rbind_SnpMatrix(m1, m2)
  expect_identical(rownames(res), c("A1", "A2", "B1", "B2"))
  expect_identical(colnames(res), c("s1", "s2"))

  m3 <- make_snpmatrix(c("C1", "C2"), c("s1", "s9"), seed = 9)
  expect_error(SNPkit:::rbind_SnpMatrix(m1, m3), "identical column names")
})

test_that("cbind_SnpMatrix is strict about row names and keeps dimnames", {
  m1 <- make_snpmatrix(c("A1", "A2"), c("s1", "s2"), seed = 10)
  m2 <- make_snpmatrix(c("A1", "A2"), c("s3", "s4"), seed = 11)
  res <- SNPkit:::cbind_SnpMatrix(m1, m2)
  expect_identical(rownames(res), c("A1", "A2"))
  expect_identical(colnames(res), c("s1", "s2", "s3", "s4"))

  m3 <- make_snpmatrix(c("X1", "X2"), c("s5", "s6"), seed = 12)
  expect_error(SNPkit:::cbind_SnpMatrix(m1, m3), "identical row names")
})
