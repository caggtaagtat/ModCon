test_that("createCodonMatrix works for example sequence", {
  x <- dim(createCodonMatrix('ATGACCGATCGAATCCGG'))
  expect_equal(x, c(2, 6))
})

