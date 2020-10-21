test_that("calculateHZEIint works for example sequence", {
  x <- round(calculateHZEIint('ATACCAGCCAGCTATTACATTT'),3)
  expect_equal(x, 24.593)
})


test_that("calculateHZEIint works for example sequence", {
  x <- calculateHZEIint('ATACCAGCCAGCTATTACATTT')
  expect_is(x, "numeric")
})

