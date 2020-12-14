context("coloc")
library(gwasglue)

fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf1 <- readVcf(fn)
vcf2 <- readVcf(fn)

test_that("coloc", {
	a <- gwasvcf_to_coloc(vcf1, vcf2, "1:1-100000000")
	expect_true(is.list(a))

	b <- expect_warning(coloc::coloc.abf(a$dataset1, a$dataset2))
	expect_true(is.list(b))
})


test_that("colic ieugwasr", {
	chrpos <- "1:109724880-109904880"
	out <- ieugwasr_to_coloc(id1='ieu-a-300', id2='ieu-a-7', chrompos=chrpos)
	b <- expect_warning(coloc::coloc.abf(out$dataset1, out$dataset2))
	expect_true(is.list(b))
})

