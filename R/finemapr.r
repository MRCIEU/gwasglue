#' Generate data for analysis in various finemapping methods
#'
#' Uses the finemapr package \url{https://github.com/variani/finemapr}.
#'
#' @param region Region of the genome to extract eg 1:109317192-110317192"
#' @param id Array of GWAS studies to query. See [`ieugwasr::gwasinfo`] 
#' for available studies
#' @param bfile If this is provided then will use the API. Default = `NULL`
#' @param plink_bin If `NULL` and bfile is not `NULL` then will detect packaged 
#' plink binary for specific OS. Otherwise specify path to plink binary. 
#' Default = `NULL`
#'
#' @export
#' @return Each id will be a list of z score data, ld matrix, and sample size
ieugwasr_to_finemapr <- function(region, id, bfile=NULL, plink_bin=NULL)
{
	id <- unique(id)
	message("Getting rsids in region")
	rsid <- ieugwasr::variants_to_rsid(region)
	message("Extracting rsids from data")
	as <- ieugwasr::associations(rsid, id, proxies=0)
	rsid_avail <- unique(as$rsid)
	message("Calculating LD for ", length(rsid_avail), " variants")
	ld <- suppressWarnings(ieugwasr::ld_matrix(rsid_avail, bfile, plink_bin, 
	                                           with_alleles=FALSE)) %>% 
	  greedy_remove()
	rsid_avail <- rownames(ld)
	message("Data available for ", length(rsid_avail), " variants")
	as <- subset(as, rsid %in% rsid_avail)
	out <- list()
	for(i in 1:length(unique(id)))
	{
		dat <- list()
		x <- as[as[["rsid"]] %in% rsid_avail & as[["id"]] == id[i],]
		dat[["z"]] <- dplyr::tibble(snp = x[["rsid"]], 
		                            zscore = x[["beta"]] / x[["se"]])
		index <- match(x[["rsid"]], rsid_avail)
		dat[["ld"]] <- ld[index, index]
		stopifnot(all(x[["rsid"]] == rownames(dat[["ld"]])))

		n <- x[["n"]]
		if(all(is.na(n)))
		{
			g <- ieugwasr::gwasinfo(id[i])
			n <- g[["sample_size"]]
		}
		dat[["n"]] <- n
		out[[id[i]]] <- dat
	}
	class(out) <- "FinemaprList"
	return(out)
}

print.FinemaprList <- function(x)
{
	utils::str(x)
}


#' Generate data for fine mapping analysis
#'
#' For a given region and VCF file, extracts the variants in the region along 
#' with LD matrix from a reference panel
#'
#' @param region Region of the genome to extract eg `1:109317192-110317192`  
#' Can be array
#' @param vcf Path to VCF file or VCF object
#' @param bfile LD reference panel
#' @param plink_bin Path to plink. Default = [`genetics.binaRies::get_plink_binary()`]
#' @param threads Number of threads to run in parallel. Default=1
#'
#' @export
#' @return List of datasets for finemapping
gwasvcf_to_finemapr <- function(region, vcf, bfile, plink_bin=genetics.binaRies::get_plink_binary(), threads=1)
{
	message("Extracting data from vcf")
	ext <- gwasvcf::query_gwas(vcf=vcf, chrompos=region)
	out <- parallel::mclapply(unique(region), function(i){
		message(i)
		m <- list()
		temp <- gwasvcf::query_gwas(vcf=ext, chrompos=i)
		m[["ld"]] <- ieugwasr::ld_matrix(names(temp), bfile=bfile, plink_bin=plink_bin, with_alleles=FALSE) %>%
			greedy_remove()
		tib <- gwasvcf::vcf_to_tibble(temp)
		m[["z"]] <- tib %>%
		subset(rsid %in% rownames(m[["ld"]])) %>%
		dplyr::mutate(z=ES/SE) %>%
		dplyr::select(snp=rsid, zscore=z)
		m[["n"]] <- tib[["SS"]]
		out[[i]] <- m
		return(out)
	}, mc.cores=threads)
	class(out) <- "FinemaprList"
	return(out)
}



greedy_remove <- function(ld)
{
	ind <- which(!is.finite(ld), arr.ind=TRUE)
	if(length(ind) == 0)
	{
		return(ld)
	}
	tab <- table(ind) %>% sort(decreasing=TRUE) %>% as.data.frame(stringsAsFactors=FALSE)
	rem <- c()
	for(i in 1:nrow(tab))
	{
		ind <- ind[!(ind[,1] == tab[["ind"]][i] | ind[,2] == tab[["ind"]][i]), ]
		rem <- c(rem, tab[["ind"]][i])
		if(nrow(ind) == 0) break
	}
	rem <- as.numeric(rem)
	ld <- ld[-rem, -rem]
	stopifnot(all(is.finite(ld)))
	return(ld)
}
