context("coloc")
library(gwasglue)

fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf1 <- VariantAnnotation::readVcf(fn)
vcf2 <- VariantAnnotation::readVcf(fn)

test_that("coloc", {
	a <- gwasvcf_to_coloc(vcf1, vcf2, "1:1-100000000")
	expect_true(is.list(a))

	b <- coloc::coloc.abf(a$dataset1, a$dataset2)
	expect_true(is.list(a))
})

