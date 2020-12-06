source("analysis/get_ovariole_number_data.R")
source("analysis/pgls_functions.R") # functions for building a PGLS dataframe
library(geiger)
library(phylolm)
downsample_factor <- 1;	fam_count_threshold <- 1 # don't downsample family level average data

args = commandArgs(trailingOnly=TRUE) # get arguments, there are 2
analysis_name <- args[1]
slope <- as.numeric(args[2])

# if the analysis is done with the family level tree, update the tree objects
if(analysis_name %in% c("allometry_body_egg_ON_family")) {
	tree <- fam_tree; genus_trees <- rep(list(tree),1000); class(genus_trees) <- "multiPhylo"
}

# function for simulating data and running a PGLS
run_simulated_pgls <- function(tree){
	source("analysis/combine_body_egg_datasets.R") # combine ON and other traits, this will randomly shuffle and choose rep. genera each iteration

	if(analysis_name == "allometry_body_egg_ON_family") {	
		residual_dataset <- body_egg_ON_full %>% 
			mutate(trait1 = logon, 
					trait2 = logvol, # egg volume
					indep = logbody) # body length
		# calculate reisduals for each trait to body size
		dataset_to_simulate <- build_resid_pgls(residual_dataset,tree) %>% mutate(X1 = trait1, 
																				X2 = trait2) # use residuals as traits to simulate
	}
	if(analysis_name == "allometry_gilbert_genus") {	
		residual_dataset <- egg_ON_gilbert_genus %>% mutate(rank = genus, 
															trait1 = logon, 
															trait2= logvol, # egg volume
															indep = logweight) # body weight
		# calculate reisduals for each trait to body size
		dataset_to_simulate <- build_resid_pgls(residual_dataset,tree) %>% mutate(X1 = trait1,
																				X2 = trait2) # use residuals as traits to simulate
	}
	
	# trim the phylogeny to only include tips with both traits
	dataset_tips <- dataset_to_simulate %>% select(rank,X1,X2) %>% na.omit() %>% pull(rank)
	sim_tree <- drop.tip(tree,tree$tip.label[which(!(tree$tip.label %in% dataset_tips))])

	# use the observed values of ovariole number / residual ovariole number as the x variable
	x <- dataset_to_simulate %>% filter(rank %in% sim_tree$tip.label) %>% pull(X1)

	# fit parameters using all available data for the trait to be simulated 
	trait_X2 <- dataset_to_simulate %>% 
					group_by(rank) %>% slice(1L) %>% # choose one representative per genus / family
					mutate(trait = X2) %>% 
					filter(rank %in% tree$tip.label) %>% # use only those taxa in the tree
					select(rank,trait) %>% na.omit() # omit missing data
	t_X2 <- trait_X2$trait;	names(t_X2) <- trait_X2$rank; t_X2 <- t_X2[order(match(names(t_X2),tree$tip.label))] # build vector for geiger
	X2_tree <- drop.tip(tree,setdiff(tree$tip.label,names(t_X2))) # trim out tips without data
	X2_fit <- fitContinuous(X2_tree,t_X2,model="BM") # simulate with a brownian motion model
	
	# use the fitted parameters to simulate the y variable
	y_rate <- X2_fit$opt$sigsq # rate of evolution
	y_anc <- X2_fit$opt$z0 # ancestral state
	y <- y_anc + slope*x + rTrait(n = 1,phy = unroot(sim_tree),model = "BM",parameters = list(sigma2 = y_rate))

	# merge the simulated data with the obesrved data frame
	data_of_interest <- data.frame(rank = sim_tree$tip.label, simX1 = x, simX2 = y)
	simulated_data <- merge(dataset_to_simulate,data_of_interest,by = "rank")

	# run a PGLS using the simulated data
	allometry_simulated <- run_all_taxa_pgls(simulated_data %>% 
							select(rank,simX1,simX2) %>% 
							rename(trait1 = simX1, trait2 = simX2),
								sim_tree)

	return(list(allometry_simulated))
}

# run the simulation and PGLS for each of 1000 trees
allometry_distribution_raw <- lapply(genus_trees,run_simulated_pgls)

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
allometry_tables <- lapply(seq(1:length(analysis_name)),get_allometry_all_taxa_table)

sink(paste(analysis_name,slope,"simulated_tables.txt",sep="_"))
for(i in seq(1:length(analysis_name))) {
	cat(analysis_name[i])
	print(allometry_tables[i])
}
sink()

save.image(paste("analysis/run_pgls_simulated",analysis_name,slope,"image.RData",sep="_"))

