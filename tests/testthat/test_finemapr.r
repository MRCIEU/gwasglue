context("finemapr")
library(ieugwasr)


test_that("ieugwasr_to_finemapr", {
	v <- ieugwasr::variants_rsid("rs7528419")
	r <- paste0(v[["chr"]], ":", v[["pos"]]-100000, "-", v[["pos"]]+100000)
	a <- ieugwasr_to_finemapr(r, c("ieu-a-7", "ieu-a-2"))
	expect_true(length(a) == 2)
	expect_true(class(a) == "FinemaprList")

	# options(finemapr_caviar = "/Users/gh13047/bin/caviar")
	# library(dplyr)
	# finemapr::run_caviar(a[["IEU-a-7"]]$z, a[["IEU-a-7"]]$ld, args = "-c 3")
})

