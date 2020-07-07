cojo_sumstat_file <- function(vcffile, outfile)
{
	vcf <- VariantAnnotation::readVcf(vcffile)
	tib <- gwasvcf::vcf_to_tibble(vcf)
	tib$LP <- 10^(-tib$LP)
	tib <- tib %>% dplyr::select(rsid, ALT, REF, AF, ES, SE, LP, SS)
	names(tib) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
	write.table(tib, file=outfile, row=F, col=T, qu=F)
	return(vcf)
}

map_variants_to_regions <- function(chrpos, pop)
{
	regionfile <- system.file("extdata", "ldetect", paste0(pop, ".bed"), package="gwasglue")
	regions <- data.table::fread(regionfile, header=TRUE) %>%
		dplyr::mutate(
			chr=as.numeric(gsub("chr", "", chr)),
			start=as.numeric(start),
			stop=as.numeric(stop)
		) %>% dplyr::as_tibble()
	
	reg <- lapply(chrpos, function(i)
	{
		x <- strsplit(i, split=":")[[1]] %>% as.numeric()
		subset(regions, chr == x[1] & x[2] >= start & x[2] <= stop)[1,] %>%
			dplyr::mutate(
				variant=i,
				region=paste0(chr, ":", start, "-", stop)
			) %>%
			dplyr::select(variant, region)
	}) %>% dplyr::bind_rows()
	return(reg)
}


cojo_cond <- function(vcffile, bfile, snplist, pop, plink=genetics.binaRies::get_plink_binary(), bcftools=genetics.binaRies::get_bcftools_binary(), gcta=genetics.binaRies::get_gcta_binary(), wd=tempdir())
{
	vcf <- cojo_sumstat_file(vcffile, file.path(wd, "sum.txt"))
	rr <- SummarizedExperiment::rowRanges(vcf) %>% dplyr::as_tibble()

	ext <- vcf[names(vcf) %in% snplist] %>%
		SummarizedExperiment::rowRanges() 
	chrpos <- paste0(SummarizedExperiment::seqnames(ext), ":", SummarizedExperiment::ranges(ext)@start)

	# Get regions
	regions <- map_variants_to_regions(chrpos, pop)
	regions$rsid <- names(ext)[match(regions$variant, chrpos)]
	dup_reg <- regions$region[duplicated(regions$region)]

	l <- list()
	for(i in dup_reg)
	{
		x <- subset(regions, region == i)
		m <- list()
		extract_list <- names(query_gwas(vcf, chrompos=i))
		write.table(extract_list, file=file.path(wd, "extract.txt"), row=F, col=F, qu=F)
		for(j in x$variant)
		{
			condsnps <- subset(x, variant != j)$rsid
			write.table(condsnps, file=file.path(wd, "cond.txt"), row=F, col=F, qu=F)
			cmd <- glue::glue("{gcta} --bfile {bfile} --extract {file.path(wd, 'extract.txt')} --cojo-file {file.path(wd, 'sum.txt')} --cojo-cond {file.path(wd, 'cond.txt')} --out {file.path(wd, 'out')}")
			system(cmd)
			m[[j]] <- fread(file.path(wd, 'out.cma'))
		}
		l[[i]] <- m
	}
	return(l)
}
