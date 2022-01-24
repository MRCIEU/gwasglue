#' Convert output from query to TwoSampleMR format
#'
#' @param x Output from ieugwasr query e.g. associations, tophits, phewas
#' @param type `"exposure"` (default) or `"outcome"`
#'
#' @export
#' @return data frame
ieugwasr_to_TwoSampleMR <- function(x, type="exposure")
{
	stopifnot(type %in% c("exposure", "outcome"))
	stopifnot(is.data.frame(x))
	names(x) <- paste0(names(x), ".", type)
	nom <- names(x)
	names(x)[nom == paste0("rsid.", type)] <- "SNP"
	names(x)[nom == paste0("ea.", type)] <- paste0("effect_allele.", type)
	names(x)[nom == paste0("nea.", type)] <- paste0("other_allele.", type)
	names(x)[nom == paste0("eaf.", type)] <- paste0("eaf.", type)
	names(x)[nom == paste0("p.", type)] <- paste0("pval.", type)
	names(x)[nom == paste0("n.", type)] <- paste0("samplesize.", type)
	names(x)[nom == paste0("trait.", type)]	<- type

	x[[paste0("mr_keep.", type)]] <- 
	  !is.na(x[[paste0("beta.", type)]]) & 
	  !is.na(x[[paste0("se.", type)]]) & 
	  !is.na(x[[paste0("effect_allele.", type)]]) & 
	  !is.na(x[["SNP"]])
	return(x)
}


#' Create exposure or outcome data format for TwoSampleMR from vcf
#'
#' @param vcf VCF object
#' @param type `"exposure"` (default) or `"outcome"`
#'
#' @export
#' @return data frame
gwasvcf_to_TwoSampleMR <- function(vcf, type="exposure")
{
	a <- vcf %>% gwasvcf::vcf_to_granges()
	if (!requireNamespace("S4Vectors", quietly = TRUE)) {
	  stop("Package \"S4Vectors\" must be installed to use this function.",
	       call. = FALSE)
	}
	S4Vectors::mcols(a)[["SNP"]] <- names(a)
	a <- dplyr::as_tibble(a)
	if(!"ES" %in% names(a)) a[["ES"]] <- NA
	if(!"SE" %in% names(a)) a[["SE"]] <- NA
	if(!"LP" %in% names(a)) a[["LP"]] <- NA
	if(!"SS" %in% names(a)) a[["SS"]] <- NA
	if(!"NC" %in% names(a)) a[["NC"]] <- NA
	if(!"id" %in% names(a)) a[["id"]] <- NA
	a[["LP"]] <- 10^-a[["LP"]]
	a[["NCONT"]] <- a[["SS"]] - a[["NC"]]
	TwoSampleMR::format_data(
		a, type=type,
		snp_col="SNP",
		effect_allele_col="ALT",
		other_allele_col="REF",
		eaf_col="AF",
		chr_col="seqnames",
		pos_col="start",
		beta_col="ES",
		se_col="SE",
		pval_col="LP",
		samplesize_col="SS",
		ncase_col="NC",
		ncontrol_col="NCONT",
		phenotype_col="id"
	)
}



#' Create a harmonised dataset from lists of vcf files
#'
#' This mimics the [`TwoSampleMR::make_dat`] function, which automatically looks 
#' up exposure and outcome datasets and harmonises them, 
#' except this function uses GWAS-VCF datasets instead.
#' The supporting reference datasets can be accessed by UoB users on BC4 
#' using [`set_bc4_files()`]
#'
#' @param id1 Exposure datasets. Either an array of vcf files, 
#' or array of IDs if vcfdir is set
#' @param id2 Outcome datasets. Either an array of vcf files, 
#' or array of IDs if vcfdir is set
#' @param proxies Lookup proxies? default=`TRUE` but requires either bfile or 
#' proxydb to be set
#' @param nthreads Parellelise default=`1`
#' @param vcfdir Location of vcf files if id1 and id2 are just IDs. \cr
#' Defaults to `options()$gwasglue.vcfdir`
#' @param proxydb Location of LD proxy database. 
#' Default=`options()$gwasglue.proxydb`
#' @param rsidx Location of rsidx index database. 
#' Default=`options()$gwasglue.rsidx`
#' @param bfile Location of LD reference panel. 
#' Default=`options()$gwasglue.bfile`
#' @param action action argument passed to [`TwoSampleMR::harmonise_data`]. 
#' The level of strictness in dealing with SNPs.
#' @param plink_bin Path to plink. 
#' Default = [`genetics.binaRies::get_plink_binary`]
#'
#' @export
#' @return harmonised dataset
make_TwoSampleMR_dat <- function(id1, id2, proxies=TRUE, nthreads=1, 
                                 vcfdir=options()$gwasglue.vcfdir, 
                                 proxydb=options()$gwasglue.proxydb, 
                                 rsidx=options()$gwasglue.rsidx, 
                                 bfile=options()$gwasglue.bfile, 
                                 action=1, 
                                 plink_bin=genetics.binaRies::get_plink_binary()
                                 )
{

	id1 <- organise_ids(id1, vcfdir)
	id2 <- organise_ids(id2, vcfdir)

	exposure_dat <- parallel::mclapply(1:nrow(id1), function(i)
	{
		message("extracting tophits for ", id1$id[i])
		d <- dirname(id1$filename[i])
		if(file.exists(file.path(d, "clump.txt")))
		{
			tophits <- scan(file.path(d, "clump.txt"), character())
			if(length(tophits) > 0)
			{
				out <- gwasvcf::query_gwas(id1$filename[i], 
				                           rsidx=rsidx, rsid=tophits) %>%
					gwasvcf_to_TwoSampleMR("exposure") %>%
					dplyr::mutate(exposure=id1$trait[i], id.exposure=id1$id[i])
			} else {
				out <- NULL
			}
		} else {
			message("Clumping ", id1$filename[i])
			stopifnot(file.exists(bfile))
			o <- gwasvcf::query_gwas(id1$filename[i], pval=pval)
			p <- gwasvcf_to_TwoSampleMR(o, "exposure")
			names(p)[names(p) == "pval.exposure"] <- "pval"
			names(p)[names(p) == "SNP"] <- "rsid"
			tophits <- ieugwasr::ld_clump(p, bfile=bfile, plink_bin=plink_bin)$rsid
			if(length(tophits) > 0)
			{
				out <- gwasvcf::query_gwas(o, rsid=tophits) %>%
					gwasvcf_to_TwoSampleMR("exposure") %>%
					dplyr::mutate(exposure=id1$trait[i], id.exposure=id1$id[i])
			} else {
				out <- NULL
			}
		}
		return(out)
	}, mc.cores=nthreads) %>% 
		dplyr::bind_rows()

	outcome_dat <- parallel::mclapply(1:nrow(id2), function(i)
	{
		gwasvcf::query_gwas(id2$filename[i], rsidx=rsidx, rsid=exposure_dat$SNP) %>%
			gwasvcf_to_TwoSampleMR("outcome") %>%
			dplyr::mutate(id.outcome=id2$id[i], outcome=id2$trait[i])
	}, mc.cores=nthreads) %>% 
		dplyr::bind_rows()

	# harmonise and perform MR
	dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action=action)
	return(dat)
}


#' Figure out specific files and IDs depending on what files exist and whether 
#' vcfdir is set
#'
#' @param id List of IDs within the vcfdir structure, 
#' or a list of GWAS VCF files, or a mixture
#' @param vcfdir Location of GWAS VCF files, 
#' or `NULL` if id is a list of actual files
#'
#' @return File paths to all datasets
organise_ids <- function(id, vcfdir)
{
	dat <- dplyr::tibble(id=id, trait=id, filename=id)
	# if id is a file
	index <- file.exists(id)
	if(!all(index))
	{
		if(is.null(vcfdir))
		{
			stop("vcfdir is not set, and the following files do not exist\n", 
			     paste(id[!index], collapse="\n"))
		} else {
			message("Constructing file names")
			dat$filename[!index] <- file.path(vcfdir, id[!index], 
			                                  paste0(id[!index], ".vcf.gz"))
			index2 <- file.exists(dat$filename)
			if(!all(index2))
			{
				stop("can't find the following files: \n", paste(id[!index], 
				                                                 collapse="\n"))
			}
			index3 <- index2 & !index
			for(i in which(index3))
			{
				jf <- file.path(dirname(dat$filename[i]), paste0(dat$id[i], ".json"))
				dat$trait[i] <- jsonlite::read_json(jf)$trait
			}
		}
	}
	return(dat)
}

#' Determine locations of useful reference datasets on bluecrystal4
#'
#' This is a convenience function for members at the University of Bristol
#' to automatically set file locations for various reference datasets. 
#' It relates only to paths on bc4.
#' 
#' @export
#' @return NULL
set_bc4_files <- function()
{
	l <- list(
		gwasglue.vcfdir="/mnt/storage/private/mrcieu/data/IGD/data/public",
		gwasglue.proxydb="/mnt/storage/private/mrcieu/research/mr-eve/mr-eve/reference/data_maf0.01_rs_ref.sqlite",
		gwasglue.rsidx="/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/1000g_filtered/annotations.vcf.gz.rsidx",
		gwasglue.bfile="/mnt/storage/private/mrcieu/research/mr-eve/mr-eve/reference/data_maf0.01_rs_ref"
	)
	options(l)
}
