# This script was written by SHC in 2019.
# Its main purpose is to reconstruct the ancestral mode of oogenesis
library(OUwie)
library(corHMM)
library(parallel)
library(geiger)
source("analysis/get_ovariole_number_data")
egg_database <- read.delim("analysis/datafiles/egg_database_final_values_July17.tsv",header=T,stringsAsFactors=F)

### These commands are used to estimate the evolutionary shifts in ovary type 
# read in table of taxonomic groups, labeled by ovary type, as recorded in Buning 1994
ovary_table <- read.delim("analysis/datafiles/mode_of_oogensis_Buning.tsv",header=T,stringsAsFactors=F)

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



