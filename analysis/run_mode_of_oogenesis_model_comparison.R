# This script was written by SHC in 2019.
# Its main purpose is to compare models of evolution, taking into account the evolution of oogenesis mode
library(OUwie)
library(corHMM)
library(parallel)
library(geiger)
source("analysis/get_ovariole_number_data.R")
g_trees <- genus_trees[1:100] # select 100 trees from the posterior distribution
downsample_factor <- 1 # do not downsample the data
fam_count_threshold <- 1 # use all families
source("analysis/combine_body_egg_datasets.R")

# read in the table from Buning, with records for ovary type by genus
ovary_table <- read.delim("analysis/datafiles/mode_of_oogensis_Buning.tsv",header=T,stringsAsFactors=F)

# set number of cores
cores <- 50

args = commandArgs(trailingOnly=TRUE)


# choose the trait of interest, argument 1
if(args[1] == "egg") { #egg volume
	analysis_name <- "egg"
	trait_dataset <- egg_database[sample(nrow(egg_database)),] %>% mutate(trait = logvol) %>%  #set the trait
		filter(!(is.na(trait))) #remove missing data from the trait dataset
} else if(args[1] == "ar") { #egg aspect ratio
	analysis_name <- "ar"
	trait_dataset <- egg_database[sample(nrow(egg_database)),] %>% mutate(trait = logar) %>% filter(!(is.na(trait)))
} else if(args[1] == "asym") { #egg asymmetry
	analysis_name <- "asym"
	trait_dataset <- egg_database[sample(nrow(egg_database)),] %>% mutate(trait = sqasym) %>% filter(!(is.na(trait))) 
} else if(args[1] == "curv") { #egg curvature
	analysis_name <- "curv"
	trait_dataset <- egg_database[sample(nrow(egg_database)),] %>% mutate(trait = sqcurv) %>% filter(!(is.na(trait))) 
} else {
	analysis_name <- "on" #ovariole number
	trait_dataset <- ON[sample(nrow(ON)),] %>% mutate(trait = logon) %>% filter(!(is.na(trait))) 
}

# The Buning tables are listed by genus, but not all genera are in the final dataset
# However, ovary type is consistent across families within Bunings set
# So when a genus is not present, for phylogenetic purposes we can choose an alternate genus from the same family

tree_dataset <- egg_database[sample(nrow(egg_database)),] %>% #start with the largest dataset of genera, families, and orders
	filter(genus %in% tree$tip.label) #use only tips in the tree

# Find one-to-one genus matches
ovary_table_genera_match <- ovary_table %>% filter(genus %in% tree_dataset$genus)
# Find families without one-to-one genus match
ovary_table_family_match <- ovary_table %>% filter(!(family %in% ovary_table_genera_match$family)) %>% 
	filter(family %in% tree_dataset$family)
# Find orders without genus or family match
ovary_table_order_match <- ovary_table %>% 	filter(!(order %in% ovary_table_genera_match$order)) %>% 
	filter(!(order %in% ovary_table_family_match$order)) %>% 
	filter(order %in% tree_dataset$order)

# Build a new ovary table with these substitutions
# First for families
new_ovary_table_family <- ovary_table_family_match %>% 
	group_by(family) %>% 
	slice(1L) %>% 
	select(-genus) %>% 
	left_join(tree_dataset %>% 
		group_by(family) %>% 
		slice(1L) %>% 
		select(family,genus),by="family") 
# Then for orders
new_ovary_table_order <- ovary_table_order_match %>% 
	group_by(order) %>% 
	slice(1L) %>% 
	select(-genus) %>% 
	left_join(tree_dataset %>% 
		group_by(order) %>% 
		slice(1L) %>% 
		select(order,genus),by="order") 
# Then combine
new_ovary_table <- bind_rows(ovary_table_genera_match,new_ovary_table_family,new_ovary_table_order)

# Remove the spurious Strepsiptera ovary type
new_ovary_table <- new_ovary_table %>% filter(ovary != "Mstar" | is.na(ovary))

### This function gets a random set of rows from a dataframe
### when passed an order name and a number of rows to retrieve 
get_rows_by_order <- function(index,order_name,ds) {
	rows <- ds %>% group_by(genus) %>% #choose one per genus
		slice(1L) %>% 
		ungroup() %>% 
		filter(order == order_name) %>% #get only rows in that insect order
		.[sample(nrow(.)),] %>% #shuffle them
		slice(1:index) #choose the first n rows
	return(rows)
}

### This function is the main model comparison
### run each iteration, over a posterior distribution of trees
run_model_comparison <- function(tree) {
	# Format the data frame for the ancestral state reconstruction on the full tree
	new_tree_dataset <- left_join(tree_dataset,new_ovary_table %>% # join ovary data
												select(ovary,genus),by="genus") %>%
		group_by(genus) %>% #choose one observation per genus
		slice(1L) %>% 
		ungroup() %>% 
		mutate(species = genus, #format for corHMM
			discrete = as.numeric(as.factor(ovary))) %>% #has missing ovary type data
		select(species,discrete) %>%
	 	as.data.frame()

	# Prune a tree to only tips in the dataset (still has NAs though, max observations are used here)
	pp_pruned <- drop.tip(tree,setdiff(tree$tip.label,new_tree_dataset$species))

	# Reconstruct ancestral ovary type (with missing data), Equal Rates model
	pp <- rayDISC(pp_pruned,na.omit(as.data.frame(new_tree_dataset[,c(1,2)])),model="ER",node.states="marginal")

	# Extract estimated tip states
	est_tips <- pp$tip.states %>% as.data.frame() %>%  
				mutate(discrete = ifelse(V1 == pmax(V1,V2,V3),1,
						ifelse(V2 == pmax(V1,V2,V3),2,3)),
						species = pp$phy$tip.label) %>% 
				select(species,discrete)

	# Join the annotated ovary dataset to the trait dataset
	new_tree_cont <- left_join(trait_dataset %>% 
									filter(genus %in% pp$phy$tip.label) %>% 
									group_by(genus) %>% 
									slice(1L) %>% ungroup() %>% 
									mutate(species = genus,continuous = trait),
								est_tips,by="species") %>% 
						select(species,discrete,continuous) %>% 
						as.data.frame()

	# Prune the tree with reconstructions to only include genera in this trait dataset without NAs
	pp_pruned_cont <- drop.tip(pp$phy,setdiff(pp$phy$tip.label,new_tree_cont$species))

	# Compare four models
	model <- c("BM1","BMS","OU1","OUM")
	root <- c(TRUE,FALSE,TRUE,TRUE) #OUwie suggests that root.station should be FALSE for BMS
	ouwie_results <- mclapply(seq(1:4),function(x){OUwie(phy=pp_pruned_cont,dat=new_tree_cont,model=model[x],root.station=root[x],simmap.tree=FALSE,diagn=FALSE)},mc.cores=cores)
	return(ouwie_results)
}

# write the results
ouwie_results <- mclapply(g_trees,run_model_comparison,mc.cores=cores)
save.image(paste(analysis_name,"ovary_model_comparison.RData",sep="_"))
