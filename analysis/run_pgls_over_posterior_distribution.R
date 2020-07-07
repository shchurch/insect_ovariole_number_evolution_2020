# This script was written by SHC in 2019.
# Its main purpose is to run phylogenetic regressions on ovariole number and egg size / number
library(phytools)
source("analysis/pgls_functions.R")
# Set factor for resampling the body size dataframe
downsample_factor <- 0.5
fam_count_threshold <- 1

analysis_names <- c("allometry_egg_ON_name",
	"allometry_egg_ON_genus",
	"allometry_body_egg_ON_family",
	"allometry_body_egg_ON_TRIG_family",
	"allometry_gilbert_name",
	"allometry_gilbert_genus",
	"allometry_fecundity_name",
	"allometry_fecundity_genus",
	"allometry_COL1",
	"allometry_ORT",
	"allometry_HYM",
	"allometry_DIP")

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

	### PGLS egg vs ON, ind body size, removing family Trigonalidae
	allometry_body_egg_ON_TRIG_family <- run_all_taxa_resid_pgls(body_egg_ON_full %>% 
							filter(!(family == "Trigonalidae")) %>% 
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


	return(list(allometry_egg_ON_name,
		allometry_egg_ON_genus,
		allometry_body_egg_ON_family,
		allometry_body_egg_ON_TRIG_family,
		allometry_gilbert_name,
		allometry_gilbert_genus,
		allometry_fecundity_name,
		allometry_fecundity_genus,
		allometry_COL1,
		allometry_ORT,
		allometry_HYM,
		allometry_DIP))
}

allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

# retrieve results from over the posterior distribution
get_allometry_all_taxa_table <- function(pgls) {
	table <- round(data.frame(
	slope_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[2]]})),
	slope_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[2]]})),
	int_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[1]]})),
	int_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[1]]})),
	pval_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]]$tTable[[8]]})),
	pval_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]]$tTable[[8]]})),
	taxa = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]]$dims$N}))),
	4)
	return(table)
}

# write results
allometry_tables <- lapply(seq(1:length(analysis_names)),get_allometry_all_taxa_table)

sink("allometry_tables.txt")
for(i in seq(1:length(analysis_names))) {
	cat(analysis_names[i])
	print(allometry_tables[i])
}
sink()

save.image("analysis/run_pgls_over_posterior_distribution.RData")
