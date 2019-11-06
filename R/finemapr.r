#' Generate data for analysis in various finemapping methods
#'
#' Uses the finemapr package https://github.com/variani/finemapr
#'
#' @param region Region of the genome to extract eg 1:109317192-110317192"
#' @param id Array of GWAS studies to query. See \code{gwasinfo} for available studies
#' @param bfile If this is provided then will use the API. Default = NULL
#' @param plink_bin If null and bfile is not null then will detect packaged plink binary for specific OS. Otherwise specify path to plink binary. Default = NULL
#'
#' @export
#' @return Each id will be a list of z score data, ld matrix, and sample size
ieugwasr_to_finemapr <- function(region, id, bfile=NULL, plink_bin=NULL)
{
	id <- unique(id)
	rsid <- ieugwasr::variants_to_rsid(region)
	ld <- ieugwasr::ld_matrix(rsid, bfile, plink_bin, with_alleles=FALSE) %>% greedy_remove()
	rsid_avail <- rownames(ld)
	as <- ieugwasr::associations(rsid_avail, id, proxies=0)
	out <- list()
	for(i in 1:length(unique(id)))
	{
		dat <- list()
		x <- as[as[["name"]] %in% rsid_avail & as[["id"]] == id[i],]
		dat[["z"]] <- dplyr::tibble(snp = x[["name"]], zscore = x[["beta"]] / x[["se"]])
		index <- match(x[["name"]], rsid_avail)
		dat[["ld"]] <- ld[index, index]
		stopifnot(all(x[["name"]] == rownames(dat[["ld"]])))

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

