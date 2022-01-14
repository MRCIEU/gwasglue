#' Write files for PWCoCo where data are read from two VCF objects or files.
#'
#' @param vcf1 VCF object or path to VCF file
#' @param vcf2 VCF object or path to VCF file
#' @param chrompos Character of the format chr:pos1-pos2
#' @param type1 How to treat vcffile1 for coloc, either `"quant"` or `"cc"`
#' @param type2 How to treat vcffile2 for coloc, either `"quant"` or `"cc"`
#' @param outfile Path to output files, without file ending
#'
#' @return `0` if success, `1` if there was a problem
#' @importFrom utils write.table
#' @importFrom gwasvcf vcflist_overlaps
gwasvcf_to_pwcoco <- function(vcf1, vcf2, chrompos, type1=NULL, type2=NULL, 
                              outfile)
{
	overlap <- gwasvcf::vcflist_overlaps(list(vcf1, vcf2), chrompos)
	vcf1 <- overlap[[1]]
	vcf2 <- overlap[[2]]
	
	if (length(vcf1) == 0 || length(vcf2) == 0)
	{
		message("No overlaps for the given chrompos in ", ifelse(length(vcf1) == 0, "vcf1", "vcf2"), ".")
		return(1)
	}
	
	# vcf1
	tib1 <- vcf1 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble() %>%
		dplyr::select(rsid, ALT, REF, AF, ES, SE, LP, SS, NC) %>%
		dplyr::rename(
			SNP = rsid,
			A1 = ALT,
			A2 = REF,
			freq = AF,
			b = ES,
			se = SE,
			p = LP,
			N = ss,
			N_case = NC
		)
	tib1$p <- 10^(-tib1$p)
	
	# Coloc type -- if study type is continuous then do not need the case column
	if (type1 == "quant" || VariantAnnotation::header(vcf1) %>% VariantAnnotation::meta() %>% {.[["SAMPLE"]][["StudyType"]]} == "Continuous")
	{
		tib1 <- tib1[c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")]
	}
	
	# vcf2
	tib2 <- vcf2 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble() %>%
		dplyr::select(rsid, ALT, REF, AF, ES, SE, LP, SS, NC) %>%
		dplyr::rename(
			SNP = rsid,
			A1 = ALT,
			A2 = REF,
			freq = AF,
			b = ES,
			se = SE,
			p = LP,
			N = ss,
			N_case = NC
		)
	tib2$p <- 10^(-tib2$p)
	
	# Coloc type -- if study type is continuous then do not need the case column
	if (type2 == "quant" || VariantAnnotation::header(vcf2) %>% VariantAnnotation::meta() %>% {.[["SAMPLE"]][["StudyType"]]} == "Continuous")
	{
		tib2 <- tib2[c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")]
	}
	
	write.table(tib1, file=paste0(outfile, "1.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
	write.table(tib1, file=paste0(outfile, "2.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
	return(0)
}

#' Write files for PWCoCo where data are read from the OpenGWAS DB.
#'
#' @param id1 ID for trait 1
#' @param id2 ID for trait 2
#' @param chrompos Character of the format chr:pos1-pos2
#' @param type1 How to treat vcffile1 for coloc, either "quant" or "cc"
#' @param type2 How to treat vcffile2 for coloc, either "quant" or "cc"
#' @param outfile Path to output files, without file ending
#'
#' @return 0 if success, 1 if there was a problem
#' @importFrom utils write.table
ieugwasr_to_pwcoco <- function(id1, id2, chrompos, type1=NULL, type2=NULL, outfile)
{
	tib1 <- ieugwasr::associations(id=id1, variants=chrompos) %>% subset(., !duplicated(rsid))
	tib2 <- ieugwasr::associations(id=id2, variants=chrompos) %>% subset(., !duplicated(rsid))
	
	if (length(tib1) < 1 || length(tib2) < 1)
	{
		message("Data could not be read using the ieugwasr package for id1 = ", id1, " and id2 = ", id2, ".")
		return(1)
	}
	
	# Matching the files is quicker for PWCoCo, so best to off-load to that?
	# Save data -- PWCoCo handles the matching and cleaning mostly by itself
	tib1 %<>% dplyr::select(rsid, ea, nea, eaf, beta, se, p, n) %>%
		dplyr::rename(
			SNP = rsid,
			A1 = ea,
			A2 = nea,
			freq = eaf,
			b = beta,
			se = se,
			p = p,
			N = n
		)
	# Need to determine whether there are cases
	info1 <- ieugwasr::gwasinfo(id1)
	if ("ncase" %in% colnames(info1))
	{
		tib1$N_case <- info1$ncase
	}
	
	tib2 %<>% dplyr::select(rsid, ea, nea, eaf, beta, se, p, n) %>%
		dplyr::rename(
			SNP = rsid,
			A1 = ea,
			A2 = nea,
			freq = eaf,
			b = beta,
			se = se,
			p = p,
			N = n
		)
	info2 <- ieugwasr::gwasinfo(id2)
	if ("ncase" %in% colnames(info2))
	{
		tib2$N_case <- info2$ncase
	}
	
	write.table(tib1, file=paste0(outfile, "1.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
	write.table(tib2, file=paste0(outfile, "2.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
	return(0)
}

#' Perform pair-wise conditional and colocalisation analysis using PWCoCo
#'
#' For a list of fine-mapped rsids, will assign to regions and generate colocalisation results for conditionally independent summary stats for each rsid
#'
#' @param id1 Path to vcffile or ID for trait1
#' @param id2 Path to vcffile2 or ID for trait2
#' @param bfile LD reference panel in Plink format (.bed, .bim, .fam)
#' @param chrompos Chromosome position (format: chr:pos1-pos2) region of interest
#' @param pwcoco Path to pwcoco binary
#' @param type1 How to treat vcffile1 for coloc, either "quant" or "cc"
#' @param type2 How to treat vcffile2 for coloc, either "quant" or "cc"
#' @param workdir Location to store files, default=tempdir()
#'
#' @export
#' @return List of colocalisation results
pwcoco <- function(id1, id2, bfile, chrompos, pwcoco, type1=NULL, type2=NULL, workdir=tempdir())
{
	if (file.exists(id1) && file.exists(id2))
	{
		message("Reading two VCF files for PWCoCo.")
		
		stopifnot(gwasvcf_to_pwcoco(id1, id2, chrompos, type1, type2, outfile=file.path(workdir, "sum_stats")) == 0)
	} else if (!file.exists(id1) && !file.exists(id2))
	{
		message("Reading two IDs from OpenGWAS.")
		
		stopifnot(ieugwasr_to_pwcoco(id1, id2, chrompos, type1, type2, outfile=file.path(workdir, "sum_stats")) == 0)
	}
	# else; mixed
	
	# PWCoCo itself is multi-threaded; is it a good idea to multi-thread this function call too?
	chr <- as.integer(strsplit(chrompos, ":")[[1]][1])
	cmd <- glue::glue("{pwcoco} --bfile {bfile} --sum_stats1 {file.path(workdir, 'sum_stats1.txt')} --sum_stats2 {file.path(workdir, 'sum_stats2.txt')} --out {file.path(workdir, 'out')} --chr {chr}")
	system(cmd)
	res <- data.table::fread(file.path(workdir, 'out.coloc'))
	
	return(res)
}
