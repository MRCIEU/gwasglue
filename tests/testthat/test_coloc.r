context("coloc")
library(gwasglue)

skip_if_not_installed("VariantAnnotation")
library(VariantAnnotation)

fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf1 <- readVcf(fn)
vcf2 <- readVcf(fn)

test_that("coloc vcf", {
	a <- gwasvcf_to_coloc(vcf1, vcf2, "1:1-100000000")
	expect_true(is.list(a))

	b <- expect_warning(coloc::coloc.abf(a$dataset1, a$dataset2))
	expect_true(is.list(b))
})


test_that("coloc ieugwasr", {
	chrpos <- "1:109724880-109904880"
	out <- expect_warning(ieugwasr_to_coloc(id1='ieu-a-300', id2='ieu-a-7', chrompos=chrpos))
	b <- expect_warning(coloc::coloc.abf(out$dataset1, out$dataset2))
	expect_true(is.list(b))
})


test_that("arth bbj", {
	chrpos <- "1:38228579-38328579"
	out <- ieugwasr_to_coloc(id1='bbj-a-73', id2='bbj-a-73', chrompos=chrpos, type1 = "cc", type2 = "cc")
	res <- coloc::coloc.abf(out$dataset1, out$dataset2)
	expect_true(res$summary[6] > 0.8)
})


test_that("coloc ieugwasr 2", {
	chrpos <- "1:47634677-47734677"
	out <- expect_warning(ieugwasr_to_coloc("ieu-a-2", "eqtl-a-ENSG00000162366", chrpos))
	res <- expect_warning(coloc::coloc.abf(out$dataset1, out$dataset2))
})



