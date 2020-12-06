# This script was written by SHC in 2020.
# Its main purpose is to run phylogenetic regressions on ovariole number and egg size / number
library(phytools)
source("analysis/pgls_functions.R")
# Set factor for resampling the body size dataframe
set.seed(12345)
downsample_factor <- 0.5
fam_count_threshold <- 1

# run PGLS over one tree from the posterior distribution
run_all_allometry_pgls <- function(tree){
	### build family body size dataset
	source("analysis/combine_body_egg_datasets.R")

	### Curculionoidea ON, egg, indep body pgls
	allometry_COL1 <- run_all_taxa_resid_pgls(COL1_body_egg_ON_genus %>% 
							select(genus,logvol,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2= logvol, indep = logbody),
								tree)

	### Orthoptera ON, egg, indep body pgls
	allometry_ORT <- run_all_taxa_resid_pgls(ORT_body_egg_ON_genus %>% 
							select(genus,logvol,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2= logvol, indep = logbody),
								tree)

	### Apoidea ON, egg, indep body pgls
	allometry_HYM <- run_all_taxa_resid_pgls(HYM_body_egg_ON_genus %>% 
							select(genus,logvol,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2= logvol, indep = logbody),
								tree)

	### Apoidea ON, egg, indep body pgls
	allometry_DIP <- run_all_taxa_resid_pgls(DIP_body_egg_ON_name %>% 
							select(name,logvol,logon,logbody,group) %>% 
							rename(rank = name, trait1 = logon, trait2= logvol, indep = logbody),
								dros_tree)

	return(list(allometry_COL1,
		allometry_ORT,
		allometry_HYM,
		allometry_DIP))
}

allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

allometry_COL1 <- lapply(allometry_distribution_raw,function(x){x[[1]]})
allometry_ORT <- lapply(allometry_distribution_raw,function(x){x[[2]]})
allometry_HYM <- lapply(allometry_distribution_raw,function(x){x[[3]]})
allometry_DIP <- lapply(allometry_distribution_raw,function(x){x[[4]]})

save.image("analysis/run_posterior_pgls_subclade.RData")