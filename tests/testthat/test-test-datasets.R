test_that("test CD_UKBB data", {
  dat <- head(CD_UKBB)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 6)
})
