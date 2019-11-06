context("TwoSampleMR")
library(gwasglue)


test_that("gwasvcf_to_TwoSampleMR", {
	fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
	vcf1 <- VariantAnnotation::readVcf(fn)
	exposure_dat <- gwasvcf_to_TwoSampleMR(vcf1)
	expect_true(nrow(exposure_dat) == nrow(vcf1))
})


test_that("ieugwasr_to_TwoSampleMR", {
	a <- ieugwasr::tophits("IEU-a-2")
	exposure_dat <- ieugwasr_to_TwoSampleMR(a)
	expect_true(nrow(exposure_dat) == nrow(a))
})

