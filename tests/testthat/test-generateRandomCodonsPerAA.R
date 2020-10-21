test_that("generateRandomCodonsPerAA works for example AA seq", {
  x <- generateRandomCodonsPerAA('Lys')
  x <- x %in% c("AAA", "AAG")
  expect_equal(x, TRUE)
})


test_that("calculateHZEIint works for example sequence", {
  x <- generateRandomCodonsPerAA(c('Lys','Lys'))
  expect_is(x, "character")
})

