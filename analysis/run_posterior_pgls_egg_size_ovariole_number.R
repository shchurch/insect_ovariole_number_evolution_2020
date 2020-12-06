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

	### PGLS on vs egg volume by insect name
	allometry_egg_ON_name <- run_all_taxa_pgls(egg_ON_name %>% 
							select(genus,logvol,logon,group) %>% 
							rename(rank = genus, trait1 = logon, trait2 = logvol),
								tree)

	### PGLS on vs egg volume by insect genus
	allometry_egg_ON_genus <- run_all_taxa_pgls(egg_ON_genus %>% 
							select(genus,logvol,logon,group) %>% 
							rename(rank = genus, trait1 = logon, trait2 = logvol),
								tree)

	### PGLS egg vs ON, independent body size
	allometry_body_egg_ON_family <- run_all_taxa_resid_pgls(body_egg_ON_full %>% 
							select(rank,logvol,logon,logbody,group) %>% 
							rename(trait1 = logon, trait2= logvol, indep = logbody),
								fam_tree)

	### Gilbert egg vs ON, ind body weight name pgls
	allometry_gilbert_name <- run_all_taxa_resid_pgls(egg_ON_gilbert_name %>% 
							select(genus,logvol,logon,logweight,group) %>% 
							mutate(rank = genus, trait1 = logon, trait2= logvol, indep = logweight),
								tree)

	### Gilbert egg vs ON, ind body weight genus pgls
	allometry_gilbert_genus <- run_all_taxa_resid_pgls(egg_ON_gilbert_genus %>% 
							select(genus,logvol,logon,logweight,group) %>% 
							mutate(rank = genus, trait1 = logon, trait2= logvol, indep = logweight),
								tree)

	return(list(allometry_egg_ON_name,
		allometry_egg_ON_genus,
		allometry_body_egg_ON_family,
		allometry_gilbert_name,
		allometry_gilbert_genus))
}

allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

allometry_egg_ON_name <- lapply(allometry_distribution_raw,function(x){x[[1]]})
allometry_egg_ON_genus <- lapply(allometry_distribution_raw,function(x){x[[2]]})
allometry_body_egg_ON_family <- lapply(allometry_distribution_raw,function(x){x[[3]]})
allometry_gilbert_name <- lapply(allometry_distribution_raw,function(x){x[[4]]})
allometry_gilbert_genus <- lapply(allometry_distribution_raw,function(x){x[[5]]})

save.image("analysis/run_posterior_pgls_egg_size_ovariole_number.RData")