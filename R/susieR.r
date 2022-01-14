#' Perform fine mapping pipeline using susieR
#'
#' Clumps data, then maps those to LD representative regions. Within each 
#' detected LD representative region, fine mapping is performed.
#'
#' @param vcffile Path to vcf file
#' @param bfile Path to ld reference panel
#' @param plink_bin Path to plink
#' @param pop EUR, ASN or AFR
#' @param threads Numeric indicating number of threads. Default `1`
#' @param clump_kb Numeric indicating kilobases used in clumping. 
#' Default `1000`
#' @param clump_r2 Numeric indicating R-squared used in clumping. Default `0.001`
#' @param clump_p Numeric indicating clumping p-value. Default `5e-8`
#' @param ... Optional arguments to be passed to susie_rss
#'
#' @export
#' @return List
susieR_pipeline <- function(vcffile, bfile, plink_bin, pop, threads=1, clump_kb=1000, clump_r2=0.001, clump_p=5e-8, ...)
{
	message("Performing clumping")
	clumped <- clump_gwasvcf(vcffile, plink_bin=plink_bin, bfile=bfile)

	message("Map clumps to regions")
	regions <- map_variants_to_regions(clumped$chrpos, pop=pop)
	regions$rsid <- clumped$rsid[match(regions$variant, clumped$chrpos)]

	message("Obtain LD matrices for each region")
	m <- gwasvcf_to_finemapr(region=regions$region, vcf=vcffile, bfile=ldref, plink_bin=plink, threads=threads)

	message("Perform susieR finemapping in each region")
	m2 <- parallel::mclapply(1:length(m), function(i) {
		message(i)
		res <- susieR::susie_rss(
			m[[i]]$z$zscore, 
			m[[i]]$ld, 
			...
		)
		res$fmset <- sapply(m[[i]]$susieR$sets$cs, function(x){
			m[[i]]$z$snp[x[which.max(m[[i]]$susieR$pip[x])]]
		})
		return(res)
	}, mc.cores=threads)

	for(i in 1:length(m))
	{
		m[[i]]$susieR <- m2[[i]]
	}
	out <- list(clumped=clumped, regions=regions, res=m)
	return(out)
}
