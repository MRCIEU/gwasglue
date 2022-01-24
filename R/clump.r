#' Perform LD clumping
#'
#' @param vcf VCF file or VCF object
#' @param clump_kb Clumping kb window. Default is very strict, `10000`
#' @param clump_r2 Clumping r2 threshold. Default is very strict, `0.001`
#' @param clump_p Clumping sig level for index variants. Default = `1` 
#' (i.e. no threshold)
#' @param pop Super-population to use as reference panel. Default = `"EUR"`. 
#' Options are `"EUR"`, `"SAS"`, `"EAS"`, `"AFR"`, `"AMR"`. 
#' `'legacy'` also available - which is a previously used verison of the EUR 
#' panel with a slightly different set of markers
#' @param bfile If this is provided then will use the API. Default = `NULL`
#' @param plink_bin If `NULL` and bfile is not `NULL` then will detect packaged 
#' plink binary for specific OS. Otherwise specify path to plink binary. 
#' Default = `NULL`
#' @param access_token Google OAuth2 access token. Used to authenticate level 
#' of access to data
#'
#' @export
#' @return data frame of clumped results
clump_gwasvcf <- function(vcf, clump_kb=1000, clump_r2=0.001, clump_p=5e-8, 
                          pop=NULL, bfile=NULL, plink_bin=NULL, 
                          access_token=NULL)
{
	message("Applying threshold to vcf")
	sig <- gwasvcf::query_gwas(vcf, pval=clump_p)

	if(is.null(bfile))
	{
		message("Using API. ",
		        "Note that this could be slow. ", 
            "To reduce server disruption it is recommended to use local LD ",
            "reference files")
		message("See gwasglue vignette on how to do this")

		fn <- function(dat)
		{
			ieugwasr::ld_clump(dat, pop=pop, clump_kb=clump_kb, clump_r2=clump_r2, 
			                   clump_p=clump_p, 
			                   access_token=ieugwasr::check_access_token())
		}
	} else {
		fn <- function(dat)
		{
			ieugwasr::ld_clump(dat, clump_kb=clump_kb, clump_r2=clump_r2, 
			                   clump_p=clump_p, bfile=bfile, plink_bin=plink_bin)
		}
	}

	# Clump
	if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
	  stop(
	    "Package \"SummarizedExperiment\" must be installed to use this function.",
	    call. = FALSE
	  )
	}
	clumped <- sig %>%
		gwasvcf::vcf_to_tibble() %>%
		dplyr::mutate(pval=10^{-LP}) %>%
		dplyr::select(rsid, pval) %>%
		fn(.) %>%
		{.$rsid} %>%
		{sig[names(sig) %in% .]} %>%
		SummarizedExperiment::rowRanges() %>%
	  {dplyr::tibble(rsid = names(.), 
		               chrpos=paste0(SummarizedExperiment::seqnames(.), ":", 
		                             SummarizedExperiment::ranges(.)@start))}
	return(clumped)
}
