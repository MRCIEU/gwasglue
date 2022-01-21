context("gassocplot")

radius <- 70000
a <- ieugwasr::tophits("ieu-a-2")
b <- ieugwasr::variants_rsid(a$rsid)
chrpos <- paste0(b$chr[1], ":", b$pos[1]-radius, "-", b$pos[1]+radius)

test_that("ieugwasr1", {

	a <- ieugwasr_to_gassocplot(chrpos, "ieu-a-2")
	expect_true(class(a) == "list")
	expect_true(all(c("data", "corr") %in% names(a)))

	skip_if_not_installed("gassocplot")
	library(gassocplot)
	b <- do.call(assoc_plot, a)
	expect_true("gtable" %in% class(b))
})



test_that("ieugwasr2", {

	a <- ieugwasr_to_gassocplot(chrpos, c("ieu-a-2", "ieu-a-7"))
	expect_true(class(a) == "list")
	expect_true(all(c("markers", "z", "corr") %in% names(a)))

	skip_if_not_installed("gassocplot")
	library(gassocplot)
	b <- do.call(stack_assoc_plot, a)
	expect_true("gtable" %in% class(b))

})

