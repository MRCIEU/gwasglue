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
	message("Getting regional variants for ", chrpos)
	variants <- variants_chrpos(chrpos)
	message("Extracting ", nrow(variants), " variants from the following datasets:\n", paste(id, collapse="\n"))
	r1 <- associations(variants[["name"]], id, proxies=0)
	r1 <- dplyr::left_join(r1, variants, by="name")
	r1 <- r1[!duplicated(paste(r1[["name"]], r1[["id"]])),]
	r1[["z"]] <- r1[["beta"]] / r1[["se"]]
	r1 <- tidyr::spread(subset(r1, select=c("name", "id", "z", "chr", "pos")), key="id", value="z")
	message("Found ", nrow(r1), " variants")
	message("Extracting LD matrix for ", nrow(r1), " variants")
	ld <- suppressWarnings(suppressMessages(ld_matrix(r1$name, with_alleles=FALSE)))
	message("Found ", nrow(ld), " variants in LD reference panel")
	r1 <- r1[match(rownames(ld), r1$name), ]
	stopifnot(all(r1$name == rownames(ld)))
	if(length(id) == 1)
	{
		list(
			data = dplyr::tibble(marker=r1$name, chr=r1$chr, pos=r1$pos, z=r1[[id]]),
			corr = ld
		) %>% return()
	} else {
		list(
			markers = dplyr::tibble(marker=r1$name, chr=r1$chr, pos=r1$pos),
			z = subset(r1, select=id),
			corr = ld,
			traits = id
		) %>% return()
	}
}

