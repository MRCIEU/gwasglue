#' Read in GWAS dataset
#'
#' @param filename Filename
#' @param skip Option passed to [`data.table::fread`]
#' @param delimiter The separator between columns in the file
#' @param gzipped Logical indicating whether the file is compressed in gz format
#' @param snp Numeric indicating column number with SNP name
#' @param nea Numeric indicating column number with NEA
#' @param ea Numeric indicating column number with effect allele
#' @param ea_af Numeric indicating column number with effect allele 
#' allele frequency
#' @param effect Numeric indicating column number with effect
#' @param se Numeric indicating column number with standard error of effect
#' @param pval Numeric indicating column number with p-value of effect
#' @param n Numeric indicating column number with sample size
#' @param info Numeric indicating column number with info
#' @param z Numeric indicating column number with z statistic of effect
#'
#' @export
#' @return data frame with log attributes
read_gwas <- function(filename, skip, delimiter, gzipped, snp, nea, ea, ea_af, 
                      effect, se, pval, n, info, z)
{
	if(gzipped)
	{
		# dat <- data.table::fread(paste0("gunzip -c ", filename), header=FALSE, skip=skip, sep=delimiter)
		dat <- data.table::fread(filename, header=FALSE, skip=skip, sep=delimiter)
	} else {
		dat <- data.table::fread(filename, header=FALSE, skip=skip, sep=delimiter)
	}
	nc <- ncol(dat)

	if(snp == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		snp <- ncol(dat)
	}
	if(nea == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		nea <- ncol(dat)
	}
	if(ea == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		ea <- ncol(dat)
	}
	if(ea_af == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		ea_af <- ncol(dat)
	}
	if(effect == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		effect <- ncol(dat)
	}
	if(se == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		se <- ncol(dat)
	}
	if(pval == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		pval <- ncol(dat)
	}
	if(n == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		n <- ncol(dat)
	}
	if(info == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		info <- ncol(dat)
	}
	if(z == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		z <- ncol(dat)
	}
	o <- TwoSampleMR::format_data(
		dat, 
		type="outcome", 
		phenotype_col="outcome",
		snp_col=names(dat)[snp],
		beta_col=names(dat)[effect],
		se_col=names(dat)[se],
		effect_allele_col=names(dat)[ea],
		other_allele_col=names(dat)[nea],
		eaf_col=names(dat)[ea_af],
		pval_col=names(dat)[pval],
		samplesize_col=names(dat)[n],
		info_col=names(dat)[info],
		z_col=names(dat)[z],
	)
	o[["beta.outcome"]][!is.finite(o[["beta.outcome"]])] <- NA
	o[["se.outcome"]][!is.finite(o[["se.outcome"]])] <- NA
	o[["pval.outcome"]][!is.finite(o[["pval.outcome"]])] <- NA
	o[["eaf.outcome"]][!is.finite(o[["eaf.outcome"]])] <- NA
	o[["info.outcome"]][!is.finite(o[["info.outcome"]])] <- NA
	o[["z.outcome"]][!is.finite(o[["z.outcome"]])] <- NA
	if(all(is.na(o[["info.outcome"]])))
	{
		o <- subset(o, select=-c(info.outcome))
	}
	if(all(is.na(o[["z.outcome"]])))
	{
		o <- subset(o, select=-c(z.outcome))
	}
	ind <- is.finite(o[["beta.outcome"]]) & is.finite(o[["pval.outcome"]])
	o <- o[ind,]
	o <- subset(o, !is.na(pval.outcome))
	return(o)
}



#' Read in reference dataset
#'
#' @param reference_file Reference vcf
#' @param rsid List of variants to read
#' @param chrompos List of chrompos to read
#' @param remove_dup_rsids =`TRUE` Remove duplicates from output
#'
#' @export
#' @return data frame
read_reference <- function(reference_file, rsid=NULL, chrompos=NULL, 
                           remove_dup_rsids=TRUE)
{
	if(!is.null(rsid))
	{
		a <- gwasvcf::query_gwas(reference_file, rsid=rsid)
	} else if(!is.null(chrompos)) {
		a <- gwasvcf::query_gwas(reference_file, chrompos=chrompos)
	} else {
		a <- gwasvcf::query_gwas(reference_file, chrompos=chrompos)
	}
	if(remove_dup_rsids)
	{
		a <- a[!duplicated(names(a)), ]
	}
	return(vcf_to_TwoSampleMR(a))
}

#' Harmonise gwas alleles to be same as reference
#'
#' @param gwas data.frame indicating GWAS
#' @param reference data.frame indicating reference dataset
#'
#' @export
#' @return data frame with attributes
harmonise_against_ref <- function(gwas, reference)
{
	# Check strand
	action <- is_forward_strand(gwas[["SNP"]], gwas[["effect_allele.outcome"]], 
	                            gwas[["other_allele.outcome"]], 
	                            reference[["SNP"]], 
	                            reference[["other_allele.exposure"]], 
	                            reference[["effect_allele.exposure"]], 
	                            threshold=0.9)

	# Harmonise
	dat <- suppressMessages(TwoSampleMR::harmonise_data(reference, gwas, action))
	cols <- c("SNP" = "ID", "effect_allele.exposure" = "ALT", 
	          "other_allele.exposure"="REF", "beta.outcome" = "BETA", 
	          "se.outcome"="SE", "pval.outcome"="PVALUE", "eaf.outcome"="AF", 
	          "samplesize.outcome"="N", "z.outcome"="ZVALUE", 
	          "info.outcome"="INFO")
	if(! "z.outcome" %in% names(dat)) dat[["z.outcome"]] <- NA
	if(! "info.outcome" %in% names(dat)) dat[["info.outcome"]] <- NA

	dat <- subset(dat, select=names(cols))
	names(dat) <- cols

	names(ref)[names(ref) == "chr.exposure"] <- "CHROM"
	names(ref)[names(ref) == "pos.exposure"] <- "POS"

	dat <- dat %>%
		dplyr::inner_join(subset(ref, select=c("SNP","other_allele.exposure",
		                                       "effect_allele.exposure","CHROM",
		                                       "POS")), 
		                  by=c("ID"="SNP", "REF"="other_allele.exposure", 
		                       "ALT"="effect_allele.exposure"))
	return(dat)
}


#' Create format for HPC pipeline 
#'
#' Takes raw files and aligns them to reference. Important if files don't have 
#' chr:pos already
#'
#' @param harmonised Output from /code{harmonise_against_ref}
#' @param path Path to write out json file and txt file
#'
#' @export
#' @return NULL
#' @importFrom utils write.table
write_out <- function(harmonised, path)
{
	j <- list(
		chr_col = 10,
		pos_col = 11,
		snp_col = 0,
		ea_col = 2,
		oa_col = 1,
		beta_col = 3,
		se_col = 4,
		ncontrol_col = 7,
		pval_col = 5,
		eaf_col = 6,
		delimiter = " ",
		header = TRUE,
		build = "GRCh37"
	)
	if(!all(is.na(harmonised[["ZVALUE"]]))) j[["imp_z_col"]] <- 8
	if(!all(is.na(harmonised[["INFO"]]))) j[["imp_z_col"]] <- 9
	if(!all(is.na(harmonised[["NCASE"]]))) {
	  j[["ncase_col"]] <- which(names(harmonised) == "NCASE")
  }
	jsonlite::write_json(j, paste0(path, ".json"), auto_unbox=TRUE)
	if(grepl(".gz$", path))
	{
		gz1 <- gzfile(path, "w")
		utils::write.table(harmonised, gz1, row.names = FALSE, col.names = TRUE, 
		                   quote = FALSE)
		close(gz1)
	} else {
		utils::write.table(harmonised, path, row.names = FALSE, col.names = TRUE, 
		                   quote = FALSE)
	}
}


#' Check a GWAS dataset against a reference known to be on the forward strand
#'
#' 
#' Assuming reference data is all on forward strand, check if 
#' the GWAS is also.
#' Use some threshold e.g. if more than 90% of alleles don't 
#' need to be flipped then it's likely that the dataset is on
#' the forward strand
#'
#' This function can be used to evaluate how strict harmonisation should be
#' The trade off if you assume we are not on the forward strand then palindromic SNPs are dropped within a particular frequency range
#' But you could instead have some small probability of error for whether palindromic SNPs are on the forward strand, and avoid dropping too many variants.
#'
#' @param gwas_snp Vector of SNP names for the dataset being checked
#' @param gwas_a1 Vector of alleles
#' @param gwas_a2 Vector of alleles
#' @param ref_snp Vector of SNP names for the reference dataset
#' @param ref_a1 Vector of alleles
#' @param ref_a2 Vector of alleles
#' @param threshold =0.9 If the proportion of allele strands match is above this threshold, then declare the dataset to be on the forward strand
#'
#' @export
#' @return 1 = Forward strand; 2 = Not on forward strand
is_forward_strand <- function(gwas_snp, gwas_a1, gwas_a2, ref_snp, ref_a1, ref_a2, threshold=0.9)
{
	requireNamespace("dplyr", quietly=TRUE)
	if(is.null(gwas_a1) | is.null(gwas_a2))
	{
		message("No info for both GWAS alleles")
		return(2)
	}

	if(1-(sum(is.na(gwas_a1)) / length(gwas_a1)) < threshold)
	{
		message("Too many missing values for gwas A1")
		return(2)
	}
	if(1-(sum(is.na(gwas_a2)) / length(gwas_a2)) < threshold)
	{
		message("Too many missing values for gwas A2")
		return(2)
	}

	g <- dplyr::data_frame(SNP=gwas_snp, A1=toupper(gwas_a1), A2=toupper(gwas_a2)) %>%
    		subset(!is.na(SNP) & !is.na(A1) & !is.na(A2))
	r <- dplyr::data_frame(SNP=ref_snp, A1=toupper(ref_a1), A2=toupper(ref_a2))

	gr <- dplyr::inner_join(g,r,by="SNP")
	index <- (gr[["A1.x"]] == gr[["A1.y"]] & gr[["A2.x"]] == gr[["A2.y"]]) | (gr[["A1.x"]] == gr[["A2.y"]] & gr[["A2.x"]] == gr[["A1.y"]])
	diff <- gr[!index,]
	diff[["A1.x"]][diff[["A1.x"]] == "C"] <- "g"
	diff[["A1.x"]][diff[["A1.x"]] == "G"] <- "c"
	diff[["A1.x"]][diff[["A1.x"]] == "T"] <- "a"
	diff[["A1.x"]][diff[["A1.x"]] == "A"] <- "t"
	diff[["A2.x"]][diff[["A2.x"]] == "C"] <- "g"
	diff[["A2.x"]][diff[["A2.x"]] == "G"] <- "c"
	diff[["A2.x"]][diff[["A2.x"]] == "T"] <- "a"
	diff[["A2.x"]][diff[["A2.x"]] == "A"] <- "t"
	diff[["A1.x"]] <- toupper(diff[["A1.x"]])
	diff[["A2.x"]] <- toupper(diff[["A2.x"]])

	index2 <- (diff[["A1.x"]] == diff[["A1.y"]] & diff[["A2.x"]] == diff[["A2.y"]]) | (diff[["A1.x"]] == diff[["A2.y"]] & diff[["A2.x"]] == diff[["A1.y"]])

	# Number that match initially
	message("SNPs that match: ", sum(index, na.rm=TRUE))
	message("SNPs that match after flipping: ", sum(index2, na.rm=TRUE))
	message("SNPs that never match: ", sum(!index2, na.rm=TRUE))

	prop <- 1 - sum(index2, na.rm=TRUE) / sum(index, na.rm=TRUE)
	message("Proportion on forward strand: ", prop)

	return(ifelse(prop > threshold, 1, 2))
}



check_null <- function(x, n)
{
	if(is.null(x))
	{
		return(rep(NA))
	}
}



#' Generic harmonisation function
#'
#' Assumes ref and alt alleles available for target and reference datasets, and uses chr:pos for alignment
#'
#' 0: stick
#' 1: swap
#' 2: rename indel
#' 3: rename indel and swap
#' 4: flip
#' 5: flip and swap
#' 6: drop (no match)
#' 7: drop (no reference)
#' 
#'
#' @param chr1 Vector
#' @param pos1 Vector
#' @param ref1 Vector
#' @param alt1 Vector
#' @param chr2 Vector
#' @param pos2 Vector
#' @param ref2 Vector
#' @param alt2 Vector
#' @param rsid2 Optional vector
#' @param indel_recode =FALSE. If TRUE then attempts to recode D/I
#' @param strand_flip =FALSE. If TRUE then attempts to flip strand when alignment is not otherwise possible
#'
#' @export
#' @return Dataframe of outcomes
harmonise <- function(chr1, pos1, ref1, alt1, chr2, pos2, ref2, alt2, rsid2 = NULL, indel_recode=FALSE, strand_flip=FALSE)
{
	chrpos1 <- paste(chr1, pos1)
	chrpos2 <- paste(chr2, pos2)
	target <- dplyr::tibble(chr=chr1, pos=pos1, ref=toupper(ref1), alt=toupper(alt1))
	reference <- dplyr::tibble(chr=chr2, pos=pos2, ref=toupper(ref2), alt=toupper(alt2))
	if(!is.null(rsid2))
	{
		reference$rsid <- rsid2
	}
	reference <- subset(reference, !duplicated(paste(chr, pos, ref, alt)))

	target$keep <- chrpos1 %in% chrpos2
	target$i <- 1:nrow(target)
	dat <- dplyr::left_join(target, reference, by=c("chr", "pos"))
	dat$decision <- 6
	dat$decision[dat$decision == 6 & is.na(dat$ref.y)] <- 7


	a0 <- dat$ref.x == dat$ref.y & dat$alt.x == dat$alt.y
	a1 <- dat$ref.x == dat$alt.y & dat$alt.x == dat$ref.y
	dat$decision[a0] <- 0
	dat$decision[a1] <- 1

	if(indel_recode)
	{
		dat$nref.y <- nchar(dat$ref.y)
		dat$nalt.y <- nchar(dat$alt.y)
		dat$inref <- dat$nref.y > dat$nalt.y
		a2 <- dat$ref.x == "I" & dat$alt.x == "D" & dat$inref | dat$ref.x == "D" & dat$alt.x == "I" & ! dat$inref
		a3 <- dat$ref.x == "I" & dat$alt.x == "D" & ! dat$inref | dat$ref.x == "D" & dat$alt.x == "I" & dat$inref		
		dat$decision[a2] <- 2
		dat$decision[a3] <- 3
	}

	if(strand_flip)
	{
		ref.x.flip <- flip_strand(dat$ref.x)
		alt.x.flip <- flip_strand(dat$alt.x)
		a5 <- dat$decision > 5 & ref.x.flip == dat$ref.y & alt.x.flip == dat$alt.y
		a6 <- dat$decision > 5 & ref.x.flip == dat$alt.y & alt.x.flip == dat$ref.y
		dat$decision[a5] <- 5
		dat$decision[a6] <- 6
	}

	dat <- dat %>% dplyr::arrange(i, decision) %>%
		subset(!duplicated(i))

	message("0: stick = ", sum(dat$decision == 0))
	message("1: swap = ", sum(dat$decision == 1))
	message("2: rename indel = ", sum(dat$decision == 2))
	message("3: rename indel and swap = ", sum(dat$decision == 3))
	message("4: flip = ", sum(dat$decision == 4))
	message("5: flip and swap = ", sum(dat$decision == 5))
	message("6: drop (no match) = ", sum(dat$decision == 6))
	message("7: drop (no reference) = ", sum(dat$decision == 7))
	return(dat)
}



flip_strand <- function(x)
{
	x %>% gsub("A", "t", .) %>%
		gsub("G", "c", .) %>%
		gsub("T", "a", .) %>%
		gsub("C", "g", .) %>%
		toupper()
}

