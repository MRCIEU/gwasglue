#' Write vcf file to cojo sumstat file
#'
#' @param vcffile Path to vcf file
#' @param outfile Path to output file
#'
#' @export
#' @return vcf object
#' @importFrom utils write.table
cojo_sumstat_file <- function(vcffile, outfile)
{
	vcf <- VariantAnnotation::readVcf(vcffile)
	tib <- gwasvcf::vcf_to_tibble(vcf)
	tib$LP <- 10^(-tib$LP)
	tib <- tib %>% dplyr::select(rsid, ALT, REF, AF, ES, SE, LP, SS)
	names(tib) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
	write.table(tib, file=outfile, row.names = FALSE, col.names = TRUE, 
	            quote = FALSE)
	return(vcf)
}

#' For a set of variants map to LD regions
#'
#' LD regions defined here 
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731402/}
#'
#' @param chrpos Array of chr:pos
#' @param pop EUR, AFR or ASN
#'
#' @export
#' @return
map_variants_to_regions <- function(chrpos, pop)
{
	regionfile <- system.file("extdata", "ldetect", paste0(pop, ".bed"), 
	                          package="gwasglue")
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


#' Perform conditional analysis using GCTA COJO
#'
#' For a list of fine-mapped rsids, will assign to regions and generate 
#' conditionally independent summary stats for each rsid.
#'
#' @param vcffile Path to vcffile
#' @param bfile LD reference panel
#' @param snplist List of rsids
#' @param pop EUR, ASN or AFR
#' @param gcta Path to gcta binary. \cr
#' For convenience can use default=`genetics.binaRies::get_gcta_binary()`
#' @param workdir Location to store temporary files. Default=`tempdir()`
#' @param threads Number of parallel threads. Default=`1`
#'
#' @export
#' @return List of independent summary stats
#' @importFrom utils write.table
cojo_cond <- function(vcffile, bfile, snplist, pop, 
                      gcta=genetics.binaRies::get_gcta_binary(), 
                      workdir=tempdir(), threads=1)
{
	message("Formatting sumstats")
	vcf <- cojo_sumstat_file(vcffile, file.path(workdir, "sum.txt"))

	ext <- vcf[names(vcf) %in% snplist] %>%
		SummarizedExperiment::rowRanges() 
	chrpos <- paste0(SummarizedExperiment::seqnames(ext), ":", 
	                 SummarizedExperiment::ranges(ext)@start)

	message("Organising regions")
	regions <- map_variants_to_regions(chrpos, pop)
	regions$rsid <- names(ext)[match(regions$variant, chrpos)]
	dup_reg <- unique(regions$region[duplicated(regions$region)])

	message(length(dup_reg), " out of ", nrow(regions), 
	        " regions have multiple variants")

	l <- parallel::mclapply(dup_reg, function(i)
	{
		message(i)
		x <- subset(regions, region == i)
		m <- list()
		y <- gwasvcf::query_gwas(vcf, chrompos=i)
		extract_list <- names(y)
		write.table(extract_list, file=file.path(workdir, "extract.txt"), 
		            row.names = FALSE, col.names = FALSE, quote = FALSE)
		for(j in x$variant)
		{
			message(j)
			condsnps <- subset(x, variant != j)$rsid
			write.table(condsnps, file=file.path(workdir, "cond.txt"), 
			            row.names = FALSE, col.names = FALSE, quote = FALSE)
			cmd <- glue::glue("{gcta} --bfile {bfile} --extract {file.path(workdir, 'extract.txt')} --cojo-file {file.path(workdir, 'sum.txt')} --cojo-cond {file.path(workdir, 'cond.txt')} --out {file.path(workdir, 'out')}")
			system(cmd)
			res <- data.table::fread(file.path(workdir, 'out.cma.cojo'))
			m[[j]] <- dplyr::select(res, rsid=SNP, chr=Chr, pos=bp, alt=refA, ES=bC, 
			                        SE=bC_se, pval=pC, n=n)
		}
		return(m)
	}, mc.cores=threads)

	message("Adding in remaining regions in the same format")
	single_reg <- regions$region[!regions$region %in% dup_reg]

	for(i in single_reg)
	{
		message(i)
		x <- subset(regions, region == i)
		j <- x$variant
		y <- gwasvcf::query_gwas(vcf, chrompos=i) %>% gwasvcf::vcf_to_tibble() %>%
			dplyr::mutate(pval=10^{-LP})
		l[[i]][[j]] <- dplyr::select(y, rsid=rsid, chr=seqnames, pos=start, 
		                             alt=ALT, ES=ES, SE=SE, pval=pval, n=SS)
	}

	return(l)
}
