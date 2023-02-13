test_that("basic regionplot works", {
  expect_error(regionplot(CD_UKBB, gene="IL23R"), NA)
})
