#' Generate coloc dataset from vcf files
#'
#' @param vcf1 VCF object or path to vcf file
#' @param vcf2 VCF object or path to vcf file
#' @param chrompos Character of chr:pos1-pos2
#'
#' @export
#' @return List of datasets to feed into coloc
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

	type1 <- ifelse(VariantAnnotation::header(vcf1) %>%
		VariantAnnotation::meta() %>%
		{.[["SAMPLE"]][["StudyType"]]} == "Continuous", "quant", "cc")

	type2 <- ifelse(VariantAnnotation::header(vcf2) %>%
		VariantAnnotation::meta() %>%
		{.[["SAMPLE"]][["StudyType"]]} == "Continuous", "quant", "cc")

	out1 <- tab1[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = type1, snp = names(vcf2)[index], z = .$ES / .$SE, chr = .$seqnames, pos = .$start, id = VariantAnnotation::samples(VariantAnnotation::header(vcf1))[1])}
	out2 <- tab2[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = type2, snp = names(vcf2)[index], z = .$ES / .$SE, chr = .$seqnames, pos = .$start, id = VariantAnnotation::samples(VariantAnnotation::header(vcf2))[1])}

	if(type1 == "cc")
	{
		out1$s <- mean(tab1$NC / tab1$SS, na.rm=TRUE)
	} else 

	if(type2 == "cc")
	{
		out2$s <- mean(tab1$NC / tab1$SS, na.rm=TRUE)
	}

	return(list(dataset1=out1, dataset2=out2))
}


#' Generate coloc dataset from the IEU GWAS database
#'
#' @param id1 ID for trait 1
#' @param id2 ID for trait 2
#' @param chrompos Character of chr:pos1-pos2
#'
#' @export
#' @return List of datasets to feed into coloc
ieugwasr_to_coloc <- function(id1, id2, chrompos)
{
	tab1 <- ieugwasr::associations(id=id1, variants=chrompos) %>% subset(., !duplicated(rsid))
	tab2 <- ieugwasr::associations(id=id2, variants=chrompos) %>% subset(., !duplicated(rsid))
	commonsnps <- tab1$rsid[tab1$rsid %in% tab2$rsid]
	tab1 <- tab1[tab1$rsid %in% commonsnps, ] %>% arrange(rsid)
	tab2 <- tab2[tab2$rsid %in% commonsnps, ] %>% arrange(rsid)
	stopifnot(all(tab1$rsid == tab2$rsid))

	index <- as.character(tab1$ea) == as.character(tab2$ea) &
			as.character(tab1$nea) == as.character(tab2$nea) &
			as.character(tab1$rsid) == as.character(tab2$rsid) &
			tab1$position == tab2$position
	stopifnot(sum(index) > 0)
	tab1$eaf <- as.numeric(tab1$eaf)
	tab2$eaf <- as.numeric(tab2$eaf)
	tab1$eaf[which(tab1$eaf > 0.5)] <- 1 - tab1$eaf[which(tab1$eaf > 0.5)]
	tab2$eaf[which(tab2$eaf > 0.5)] <- 1 - tab2$eaf[which(tab2$eaf > 0.5)]

	info1 <- ieugwasr::gwasinfo(id1)
	type1 <- info1 %>% {ifelse(.$unit == "logOR", "cc", "quant")}
	info2 <- ieugwasr::gwasinfo(id2)
	type2 <- info2 %>% {ifelse(.$unit == "logOR", "cc", "quant")}

	tab1 <- tab1[index,] %>% {list(pvalues = .$p, N = .$n, MAF = .$eaf, beta = .$beta, varbeta = .$se^2, type = type1, snp = .$rsid, z = .$beta / .$se, chr = .$chr, pos = .$position, id = id1)}
	tab2 <- tab2[index,] %>% {list(pvalues = .$p, N = .$n, MAF = .$eaf, beta = .$beta, varbeta = .$se^2, type = type2, snp = .$rsid, z = .$beta / .$se, chr = .$chr, pos = .$position, id = id2)}

	if(type1 == "cc")
	{
		tab1$s <- info1$ncase / info1$sample_size
	} else 

	if(type2 == "cc")
	{
		tab2$s <- info2$ncase / info2$sample_size
	}

	return(list(dataset1=tab1, dataset2=tab2))
}


