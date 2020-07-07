# This script was written by SHC in 2019.
# Its main purpose is to run phylogenetic regressions on ovariole number and body size
library(phytools)
source("analysis/pgls_functions.R")

# Set factor for resampling the body size dataframe
downsample_factor <- 0.5
fam_count_threshold <- 1

body_analysis_names <- c("allometry_body_egg_ON_family",
	"allometry_gilbert_name",
	"allometry_gilbert_genus",
	"allometry_COL1",
	"allometry_ORT",
	"allometry_HYM",
	"allometry_DIP")

# run PGLS over one tree from the posterior distribution
run_all_allometry_pgls <- function(tree){
	### build family body size dataset, repeat each iteration to account for phenotypic uncertainty
	source("analysis/combine_body_egg_datasets.R")


	### PGLS body length vs ON
	allometry_body_egg_ON_family <- run_all_taxa_pgls(body_egg_ON_full %>% 
							select(rank,logvol,logon,logbody,group) %>% 
							rename(trait1 = logon, trait2 = logbody),
								fam_tree)


	### Gilbert body weight vs ON, species level
	allometry_gilbert_name <- run_all_taxa_pgls(egg_ON_gilbert_name %>% 
							select(genus,logvol,logon,logweight,group) %>% 
							mutate(rank = genus, trait1 = logon, trait2 = logweight),
								tree)

	### Gilbert body weight, genus level
	allometry_gilbert_genus <- run_all_taxa_pgls(egg_ON_gilbert_genus %>% 
							select(genus,logvol,logon,logweight,group) %>% 
							mutate(rank = genus, trait1 = logon, trait2 = logweight),
								tree)

	### Curculionoidea body vs ON
	allometry_COL1 <- run_all_taxa_pgls(COL1_body_egg_ON_genus %>% 
							select(genus,logvol,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2 = logbody),
								tree)

	### Orthoptera body vs ON
	allometry_ORT <- run_all_taxa_pgls(ORT_body_egg_ON_genus %>% 
							select(genus,logvol,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2 = logbody),
								tree)

	### Apoidea body vs ON
	allometry_HYM <- run_all_taxa_pgls(HYM_body_egg_ON_genus %>% 
							select(genus,logvol,logon,logbody,group) %>% 
							rename(rank = genus, trait1 = logon, trait2 = logbody),
								tree)

	### Drosophilidae body vs ON
	allometry_DIP <- run_all_taxa_pgls (DIP_body_egg_ON_name %>% 
							select(name,logvol,logon,logbody,group) %>% 
							rename(rank = name, trait1 = logon, trait2 = logbody),
								dros_tree)


	return(list(allometry_body_egg_ON_family,
		allometry_gilbert_name,
		allometry_gilbert_genus,
		allometry_COL1,
		allometry_ORT,
		allometry_HYM,
		allometry_DIP))
}

body_allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

# retrieve results from over the posterior distribution
get_allometry_all_taxa_table <- function(pgls) {
	table <- round(data.frame(
	slope_min = min(sapply(body_allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[2]]})),
	slope_max = max(sapply(body_allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[2]]})),
	int_min = min(sapply(body_allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[1]]})),
	int_max = max(sapply(body_allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[1]]})),
	pval_min = min(sapply(body_allometry_distribution_raw,function(x) {x[[pgls]]$tTable[[8]]})),
	pval_max = max(sapply(body_allometry_distribution_raw,function(x) {x[[pgls]]$tTable[[8]]})),
	taxa = min(sapply(body_allometry_distribution_raw,function(x) {x[[pgls]]$dims$N}))),
	4)
	return(table)
}

# write results
body_allometry_tables <- lapply(seq(1:length(body_analysis_names)),get_allometry_all_taxa_table)

sink("body_allometry_tables.txt")
for(i in seq(1:length(body_analysis_names))) {
	cat(body_analysis_names[i])
	print(body_allometry_tables[i])
}
sink()

save.image("analysis/run_pgls_body_ovariole_number.RData")
