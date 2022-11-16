test_that("basic regionplot works", {
  expect_error(regionplot(CD_UKBB, gene="IL23R"), NA)
})

test_that("advanecd regionplot works", {
  expect_no_error(regionplot(CD_UKBB, gene="IL23R"))
})
