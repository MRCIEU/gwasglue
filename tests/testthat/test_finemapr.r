context("finemapr")
library(ieugwasr)


# v <- variants_rsid("rs7528419")
# r <- paste0(v[["chr"]], ":", v[["pos"]]-500000, "-", v[["pos"]]+500000)
# a <- create_finemapr_data(r, c("IEU-a-7", "IEU-a-2"), bfile, plink_bin)
# finemapr::run_caviar(a[["IEU-a-7"]]$z, a[["IEU-a-7"]]$ld, args = "-c 3")


# function(id = c("IEU-a-2", "IEU-a-7"), region = r)
# {

# 	a2 <- associations(region, id)
# 	snplist <- Reduce(intersect, lapply(id, function(x) a2[a2[["id"]] == x, ][["name"]]))

# 	a2 <- a2[a2[["name"]] %in% snplist, ]

# 	lapply(id, function(x)
# 	{
# 		a3 <- a2[a2[["id"]] == x, ]
# 		dplyr::tibble(SNP=a3[["name"]], A1=a3[["ea"]], A2=a3[["nea"]], freq=a3[["eaf"]], b=a3[["beta"]], se=a3[["se"]], p=a3[["p"]], n=a3[["n"]]) %>% write.table(., x, row.names=FALSE, col.names=TRUE, quote=FALSE)
# 	})

# }

