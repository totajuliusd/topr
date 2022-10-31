test_that("basic manhattan works", {
  expect_error(manhattan(CD_UKBB), NA)
})


test_that("adcvanced manhattan works", {
  expect_no_error(manhattan(CD_UKBB))
})