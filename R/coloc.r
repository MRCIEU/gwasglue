#' Perform coloc
#'
#' Perform coloc between two vcfs given a particular locus
#'
#' @param vcf1 VCF object
#' @param vcf2 VCF object
#' @param chrompos chrompos object (see input to parse_chrompos)
#'
#' @export
#' @return Result from coloc
gwasvcf_to_coloc <- function(vcf1, vcf2, chrompos)
{
	## TODO: binary or quantitative traits
	## TODO: multiallelic variants

	o <- gwasvcf::vcflist_overlaps(list(vcf1, vcf2), chrompos)
	vcf1 <- o[[1]]
	vcf2 <- o[[2]]

	stopifnot(length(vcf1) == length(vcf2))
	tab1 <- vcf1 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
	tab2 <- vcf2 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
	index <- as.character(tab1$REF) == as.character(tab2$REF) &
			as.character(tab1$ALT) == as.character(tab2$ALT) &
			as.character(tab1$seqnames) == as.character(tab2$seqnames) &
			tab1$start == tab2$start
	stopifnot(sum(index) > 0)

	tab1 <- tab1[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = "quant", snp = names(vcf2)[index])}
	tab2 <- tab2[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = "quant", snp = names(vcf2)[index])}

	return(list(dataset1=tab1, dataset2=tab2))
}



## TODO: add ieugwasr version
