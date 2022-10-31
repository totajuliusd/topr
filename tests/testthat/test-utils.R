test_that("get_lead_snps works", {
  out <- get_lead_snps(CD_UKBB)
  #expect_equal(2 * 2, 4)
})


test_that("annotate_with_nearest_gene works", {
  # missing test
})


test_that("get_gene_coords works", {
  # get_gene_coords("IL23R")
})


test_that("get_snps_within_region works", {
  # get_snps_within_region(CD_UKBB, region = "chr1:67138906-67259979")
})


test_that("get_top_snp works", {
  # get_top_snp(CD_UKBB, chr="chr1")
})
