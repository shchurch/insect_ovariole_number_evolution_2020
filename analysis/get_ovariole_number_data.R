# This script was written by SHC in 2019.
# Its main purpose is to read in and format the ovariole number dataset
library(dplyr)
library(ggtree)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(phytools)

# read in data from google doc
raw_ON <- read.delim(file="analysis/datafiles/ovariole_number.tsv",header=TRUE,stringsAsFactors=FALSE)

# read in list of orders and groups = monophyletic large lineages
groups_orders <- read.delim(file="analysis/datafiles/groups_orders.txt",header=T,stringsAsFactors=F)
group_levels = c("Apterygota","Palaeoptera","Polyneoptera","Condylognatha","Psocodea","Amphiesmenoptera","Antliophora","Neuropteroidea","Hymenoptera")

# set colors for figures
mrk <- c("#9c9c9c","#005487","#7570B3","#679c40","#A6761D","#52bcd4","#E6AB02","#1B9E77","#e37870","#9c9c9c","#2a2a2a")
names(mrk) <- c(group_levels,"above","below")

# read in TaxReformer updated taxonomy
matched_taxonomy <- read.delim(file="analysis/datafiles/ovariole_names_matched_June2020.tsv",header=T,stringsAsFactors=F)

# select only relevant taxonomic names
matched <- matched_taxonomy %>% select(order,family,subfamily,suborder,superfamily,superorder,tribe,subtribe,subsection,updated_genus,updated_species,updated_subspecies,ID,updated_fullname,ncbi_id) %>% 
	rename(updated_order = order) %>%
	filter(!(ID %in% (raw_ON %>% filter(problem == "taxonomy") %>% pull(ID)))) # remove the taxonomic information for entries where the information was found to be wrong (e.g. Diplura entries, where the order is incorrect from global names)

# add a name, combining genus and species
ON <- raw_ON %>% 
	mutate(median_on = ((maximum_on - minimum_on)/2) + minimum_on) %>% # if a range, calculate median
	# choose a value per record, with priority given to average, range, then single record
	mutate(on = ifelse(!(is.na(mean_on)), mean_on,
	       ifelse (!(is.na(median_on)), median_on, single_on))) %>% 
	mutate(logon = log10(on),sqrton = sqrt(on)) %>% # transform values
	left_join(.,matched,by="ID") %>% # match to taxonomic records
	mutate(taxonomy = ifelse((updated_fullname == "" | is.na(updated_fullname)),"original","updated")) %>%
	mutate(genus = ifelse(taxonomy == "updated",updated_genus,genus),
		species = ifelse(taxonomy == "updated",updated_species,species),
		subspecies = ifelse(taxonomy == "updated",updated_subspecies,subspecies),
		order = ifelse(taxonomy == "updated",updated_order,order)) %>% 
	mutate(name = ifelse(species != "",paste(genus,species,sep="_"),genus)) %>% # create name, if species recorded
	mutate(group = plyr::mapvalues(order,groups_orders$order,groups_orders$group,warn_missing=F)) %>% # add lineage groups
	# remove entries that have a flagged problem, such as being a worker individually, or raised in starving conditions
	filter(!(problem %in% c("worker","description")))

# code used for identifying order problems originally
# 	order_problems <- ON %>% 
# 				filter(order != updated_order) %>% 
# 				arrange(publication) %>% 
# 				select(publication,name,updated_fullname,order,updated_order,ID)

# read in phylogenies
## posterior distributions
genus_trees <- read.nexus("phylogeny/mrbayes_stitched_misof_posterior_ultrametric.nxs")

## MCC trees
genus_mcc_tree <- read.nexus("phylogeny/misof_mcc_ultrametric.nxs")

## family level tree from the Rainford et al 2014 publication
fam_tree <- read.nexus("phylogeny/rainford_family_ultrametric.nxs")

## Drosophilidae speices tree (original to this study)
dros_tree <- read.tree("phylogeny/Drosophilidae_time_calibrated.tre")

## choose a tree for downstream analyses
tree <- genus_mcc_tree
trees <- genus_trees

### code for cleaning up trees originally
#	library(ape);library(phytools);library(phangorn)
#	genus_trees <- read.nexus("phylogeny/mrbayes_stitched_misof_posterior.tre")
#	### ladderize them for easy visualization
#	genus_trees <- lapply(genus_trees,ladderize)
#	### sub out the spaces in tip labels
#	for(i in 1:length(genus_trees)){
#		genus_trees[[i]]$tip.label <- gsub("_[0-9]+","",genus_trees[[i]]$tip.label)
#	
#		### force family tree to be ultrametric
#		genus_trees[[i]]<-nnls.tree(cophenetic(genus_trees[[i]]),genus_trees[[i]],rooted=TRUE)
#	}
#	### write ultrametric trees
#	write.nexus(file="phylogeny/mrbayes_stitched_misof_posterior_ultrametric.nxs",genus_trees)

