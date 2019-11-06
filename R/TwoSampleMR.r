#' Convert output from query to TwoSampleMR format
#'
#' @param x Output from ieugwasr query e.g. associations, tophits, phewas
#' @param type "exposure" (default) or "outcome"
#'
#' @export
#' @return data frame
ieugwasr_to_TwoSampleMR <- function(x, type="exposure")
{
	stopifnot(type %in% c("exposure", "outcome"))
	stopifnot(is.data.frame(x))
	names(x) <- paste0(names(x), ".", type)
	nom <- names(x)
	names(x)[nom == paste0("name.", type)] <- "SNP"
	names(x)[nom == paste0("ea.", type)] <- paste0("effect_allele.", type)
	names(x)[nom == paste0("nea.", type)] <- paste0("other_allele.", type)
	names(x)[nom == paste0("eaf.", type)] <- paste0("eaf.", type)
	names(x)[nom == paste0("p.", type)] <- paste0("pval.", type)
	names(x)[nom == paste0("n.", type)] <- paste0("samplesize.", type)
	names(x)[nom == paste0("trait.", type)]	<- type

	x[[paste0("mr_keep.", type)]] <- !is.na(x[[paste0("beta.", type)]]) & !is.na(x[[paste0("se.", type)]]) & !is.na(x[[paste0("effect_allele.", type)]]) & !is.na(x[["SNP"]])
	return(x)
}


#' Create exposure or outcome data format for TwoSampleMR from vcf
#'
#' @param vcf VCF object
#' @param type ="exposure" or "outcome"
#'
#' @export
#' @return data frame
gwasvcf_to_TwoSampleMR <- function(vcf, type="exposure")
{
	a <- vcf %>% gwasvcf::vcf_to_granges()
	S4Vectors::mcols(a)[["SNP"]] <- names(a)
	a <- dplyr::as_tibble(a)
	if(!"ES" %in% names(a)) a[["ES"]] <- NA
	if(!"SE" %in% names(a)) a[["SE"]] <- NA
	if(!"LP" %in% names(a)) a[["LP"]] <- NA
	if(!"SS" %in% names(a)) a[["SS"]] <- NA
	if(!"NC" %in% names(a)) a[["NC"]] <- NA
	if(!"id" %in% names(a)) a[["id"]] <- NA
	a[["LP"]] <- 10^-a[["LP"]]
	a[["NCONT"]] <- a[["SS"]] - a[["NC"]]
	TwoSampleMR::format_data(
		a, type=type,
		snp_col="SNP",
		effect_allele_col="ALT",
		other_allele_col="REF",
		eaf_col="AF",
		chr_col="CHROM",
		pos_col="POS",
		beta_col="ES",
		se_col="SE",
		pval_col="LP",
		samplesize_col="SS",
		ncase_col="NC",
		ncontrol_col="NCONT",
		phenotype_col="id"
	)
}
