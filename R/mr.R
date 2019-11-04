#' Convert output from query to TwoSampleMR format
#'
#' @param x Output from ieugwasr query e.g. associations, tophits, phewas
#' @param type "exposure" (default) or "outcome"
#'
#' @export
#' @return data frame
ieu_to_TwoSampleMR <- function(x, type="exposure")
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
