context("test-dummy.R")

test_that("dummy works", {
  dummy_var <- 1
  expect_is(dummy_var, "numeric")
  expect_true(dummy_var == 1)
})
