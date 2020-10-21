test_that("calculateMaxEntScanScore works for example SA", {
  x <- calculateMaxEntScanScore('TTCCAAACGAACTTTTGTAGGGA',3)
  expect_equal(x, "2.89")
})

test_that("calculateMaxEntScanScore works for example SD", {
  x <- calculateMaxEntScanScore('GAGGTAAGT',5)
  expect_equal(x, "11.08")
})


