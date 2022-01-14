#' Generate regional plot for ieugwasr
#'
#' Uses James Staley's gassocplot package https://github.com/jrs95/gassocplot
#'
#' @param chrpos A window range to plot e.g. 16:3349655-3849655
#' @param id Vector of one or more IEU GWAS db study IDs
#'
#' @export
#' @return assoc_plot or stack_assoc_plot if multiple markers given
ieugwasr_to_gassocplot <- function(chrpos, id)
{
	stopifnot(length(chrpos) == 1)
	r1 <- ieugwasr::associations(chrpos, id, proxies=0)
	r1 <- r1[!duplicated(paste(r1[["rsid"]], r1[["id"]])),]
	r1[["z"]] <- r1[["beta"]] / r1[["se"]]
	r1 <- tidyr::spread(subset(r1, 
	                           select=c("rsid", "id", "z", "chr", "position")), 
	                    key="id", value="z")
	message("Found ", nrow(r1), " variants")
	message("Extracting LD matrix for ", nrow(r1), " variants")
	ld <- suppressWarnings(suppressMessages(
		ieugwasr::ld_matrix(r1[["rsid"]], with_alleles=FALSE)
	))
	message("Found ", nrow(ld), " variants in LD reference panel")
	r1 <- r1[match(rownames(ld), r1[["rsid"]]), ]
	stopifnot(all(r1[["rsid"]] == rownames(ld)))
	if(length(id) == 1)
	{
		list(
			data = dplyr::tibble(marker=r1[["rsid"]], chr=r1[["chr"]], 
			                     pos=r1[["position"]], z=r1[[id]]),
			corr = ld
		) %>% return()
	} else {
		list(
			markers = dplyr::tibble(marker=r1[["rsid"]], chr=r1[["chr"]], pos=r1[["position"]]),
			z = subset(r1, select=id),
			corr = ld,
			traits = id
		) %>% return()
	}
}


#' Convert coloc dataset to gassocplot dataset
#'
#' @param coloclist Output from *_to_coloc
#' @param bfile If number of SNPs > 500 then need to provide your own LD reference panel. Provide plink dataset here. 
#' @param plink_bin If number of SNPs > 500 then need to provide your own LD reference panel. Provide plink executable here
#'
#' @export
#' @return List to feed into gassocplot
coloc_to_gassocplot <- function(coloclist, bfile=NULL, plink_bin=NULL)
{
	markers <- dplyr::tibble(
		marker = coloclist$dataset1$snp,
		chr = coloclist$dataset1$chr,
		pos = coloclist$dataset1$pos,
	)
	z <- dplyr::tibble(
		id1 = coloclist$dataset1$z,
		id2 = coloclist$dataset2$z
	)
	message("Extracting LD matrix for ", nrow(markers), " variants")
	ld <- ieugwasr::ld_matrix(markers[["marker"]], with_alleles=FALSE, bfile=bfile, plink_bin=plink_bin)
	message("Found ", nrow(ld), " variants in LD reference panel")
	index <- match(rownames(ld), markers[["marker"]])
	markers <- markers[index, ]
	z <- z[index, ]
	stopifnot(all(markers$marker == rownames(ld)))
	traits <- c(coloclist[["dataset1"]][["id"]], coloclist[["dataset2"]][["id"]])
	names(z) <- traits

	list(markers = markers, z = z, corr = ld, traits = traits) %>% return()
}
