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

	### PGLS body length vs ON
	allometry_body_ON_family <- run_all_taxa_pgls(body_egg_ON_full %>% 
							select(rank,logon,logbody,group) %>% 
							rename(trait1 = logon, trait2 = logbody),
								fam_tree)


	### Gilbert body weight vs ON, species level
	allometry_body_gilbert_name <- run_all_taxa_pgls(egg_ON_gilbert_name %>% 
							select(genus,logon,logweight,group) %>% 
							mutate(rank = genus, trait1 = logon, trait2 = logweight),
								tree)

	### Gilbert body weight, genus level
	allometry_body_gilbert_genus <- run_all_taxa_pgls(egg_ON_gilbert_genus %>% 
							select(genus,logon,logweight,group) %>% 
							mutate(rank = genus, trait1 = logon, trait2 = logweight),
								tree)

	### Curculionoidea body vs ON
	allometry_body_COL1 <- run_all_taxa_pgls(COL1_body_egg_ON_genus %>% 
							select(genus,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2 = logbody),
								tree)

	### Orthoptera body vs ON
	allometry_body_ORT <- run_all_taxa_pgls(ORT_body_egg_ON_genus %>% 
							select(genus,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2 = logbody),
								tree)

	### Apoidea body vs ON
	allometry_body_HYM <- run_all_taxa_pgls(HYM_body_egg_ON_genus %>% 
							select(genus,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2 = logbody),
								tree)

	### Drosophilidae body vs ON
	allometry_body_DIP <- run_all_taxa_pgls (DIP_body_egg_ON_name %>% 
							select(name,logon,logbody,group) %>% 
							rename(rank = name, trait1 = logon, trait2 = logbody),
								dros_tree)

	return(list(allometry_body_ON_family,
		allometry_body_gilbert_name,
		allometry_body_gilbert_genus,
		allometry_body_COL1,
		allometry_body_ORT,
		allometry_body_HYM,
		allometry_body_DIP))
}

allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

allometry_body_ON_family <- lapply(allometry_distribution_raw,function(x){x[[1]]})
allometry_body_gilbert_name <- lapply(allometry_distribution_raw,function(x){x[[2]]})
allometry_body_gilbert_genus <- lapply(allometry_distribution_raw,function(x){x[[3]]})
allometry_body_COL1 <- lapply(allometry_distribution_raw,function(x){x[[4]]})
allometry_body_ORT <- lapply(allometry_distribution_raw,function(x){x[[5]]})
allometry_body_HYM <- lapply(allometry_distribution_raw,function(x){x[[6]]})
allometry_body_DIP <- lapply(allometry_distribution_raw,function(x){x[[7]]})

save.image("analysis/run_posterior_pgls_body_size.RData")