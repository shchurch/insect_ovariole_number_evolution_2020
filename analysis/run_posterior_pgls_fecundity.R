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

	### Gilbert fecundity vs ON name
	allometry_fecundity_name <- run_all_taxa_pgls(ON_gilbert_name %>% 
							select(genus,logfecundity,logon,group) %>% 
							mutate(rank = genus, trait1 = logon, trait2= logfecundity),
								tree)

	### Gilbert fecundity vs ON genus
	allometry_fecundity_genus <- run_all_taxa_pgls(ON_gilbert_genus %>% 
							select(genus,logfecundity,logon,group) %>% 
							mutate(rank = genus, trait1 = logon, trait2= logfecundity),
								tree)

	return(list(allometry_fecundity_name,
		allometry_fecundity_genus))
}

allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

allometry_fecundity_name <- lapply(allometry_distribution_raw,function(x){x[[1]]})
allometry_fecundity_genus <- lapply(allometry_distribution_raw,function(x){x[[2]]})

save.image("analysis/run_posterior_pgls_fecundity.RData")