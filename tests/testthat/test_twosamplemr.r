context("TwoSampleMR")


test_that("gwasvcf_to_TwoSampleMR", {
  skip_if_not_installed("S4Vectors")
	fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
	vcf1 <- VariantAnnotation::readVcf(fn)
	exposure_dat <- gwasvcf_to_TwoSampleMR(vcf1)
	expect_true(nrow(exposure_dat) == nrow(vcf1))
})


test_that("ieugwasr_to_TwoSampleMR", {
	a <- ieugwasr::tophits("ieu-a-2")
	b <- ieugwasr::associations(a$rsid, "ieu-a-7")
	exposure_dat <- ieugwasr_to_TwoSampleMR(a)
	outcome_dat <- ieugwasr_to_TwoSampleMR(b, type="outcome")
	dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat)
	out <- TwoSampleMR::mr(dat)
	expect_true(nrow(exposure_dat) == nrow(dat))
	expect_true(nrow(out) > 3)
})

