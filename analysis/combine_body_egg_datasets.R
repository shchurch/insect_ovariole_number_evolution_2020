# This script was written by SHC in 2019.
# Its main purpose is to combine the ovariole number dataset with other published data

# read in insect egg size from Church et al 2019 Scientific Data
egg_database <- read.delim("analysis/datafiles/egg_database_final_values_July2019.tsv",header=T,stringsAsFactors=F)

# read in body size data for specific groups with ovariole number
raw_bodysize <- read.delim(file="analysis/datafiles/body_size.tsv",header=TRUE,stringsAsFactors=FALSE)
matched_body_taxonomy <- read.delim(file="analysis/datafiles/body_matched_taxonomy_Sep2019.tsv",header=T,stringsAsFactors=F) %>% filter(problem != "no_taxonomy")

matched_body <- matched_body_taxonomy %>% select(order,family,subfamily,suborder,superfamily,tribe,updated_genus,updated_species,updated_subspecies,ID,updated_fullname,ncbi_id) %>% 
		rename(updated_order = order)

body <- raw_bodysize %>% 
		mutate(median_body = ((maximum_body - minimum_body)/2) + minimum_body) %>% 
		mutate(body = ifelse( !(is.na(mean_body)), mean_body,
              ifelse ( !(is.na(median_body)), median_body, single_body))) %>% 
		mutate(logbody = log10(body)) %>% 
		left_join(.,matched_body,by="ID") %>%
		mutate(taxonomy = ifelse((updated_fullname == "" | is.na(updated_fullname)),"original","updated")) %>%
		mutate(genus = ifelse(taxonomy == "updated",updated_genus,genus),
			species = ifelse(taxonomy == "updated",updated_species,species),
			subspecies = ifelse(taxonomy == "updated",updated_subspecies,subspecies),
			order = ifelse(taxonomy == "updated",updated_order,order)) %>% 
		mutate(name = ifelse(species != "",paste(genus,species,sep="_"),genus)) %>%
		mutate(group = plyr::mapvalues(order,groups_orders$order,groups_orders$group,warn_missing=F)) 

# build data frames combining egg size and ovariole number
egg_name <- egg_database[sample(1:nrow(egg_database)),] %>% # shuffle the dataset
			filter(species != "") %>% # choose only egg datapoints with full name
			filter(!(is.na(vol))) %>% # choose only egg datapoints with volume
			group_by(name) %>% # group by full name
			slice(1L) %>% ungroup() # slice random top row

egg_genus <- egg_database[sample(1:nrow(egg_database)),] %>% # shuffle the dataset
		filter(!(is.na(vol))) %>% # choose only egg datapoints with volume
		group_by(genus) %>% # group by genus
		slice(1L) %>% ungroup() # slice random top row

		# repeat for ovariole number
ON_name <- ON[sample(1:nrow(ON)),] %>%
		filter(species != "") %>%
		filter(!(is.na(on))) %>% 
		group_by(name) %>% 
		slice(1L) %>% ungroup()

ON_genus <- ON[sample(1:nrow(ON)),] %>%
		filter(!(is.na(on))) %>% 
		group_by(genus) %>% 
		slice(1L) %>% ungroup()

body_name <- body[sample(1:nrow(body)),] %>% # shuffle the dataset
		filter(species != "") %>% # choose only datapoints with full name
		filter(!(is.na(body))) %>% # choose only data points with body size
		group_by(name) %>% # group by insect name
		slice(1L) %>% ungroup() # the slice off the first row per group, which is random cause we shuffled the data

body_genus <- body[sample(1:nrow(body)),] %>% # shuffle the dataset
		filter(!(is.na(body))) %>% # choose only data points with body size
		group_by(genus) %>% # group by insect genus
		slice(1L) %>% ungroup() # the slice off the first row per group, which is random cause we shuffled the data


# NAMES WITH EGG AND ON DATA
egg_ON_name <- left_join(ON_name,select(egg_name,c(vol,logvol,name)),by="name") %>% filter(!is.na(vol))

# GENERA WITH EGG AND ON DATA
egg_ON_genus <- left_join(ON_genus,select(egg_genus,c(vol,logvol,genus)),by="genus") %>% filter(!is.na(vol))

# NAMES WITH BODY AND ON DATA
body_ON_name <- left_join(ON_name,select(body_name,c(body,logbody,name)),by="name") %>% filter(body != "")
# NAMES WITH BODY, EGG, AND ON DATA
body_egg_ON_name <- left_join(egg_ON_name,select(body_name,c(body,logbody,name)),by="name") %>% filter(body != "")

# GENERA WITH BODY AND ON DATA
body_ON_genus <- left_join(ON_genus,select(body_genus,c(body,logbody,genus)),by="genus") %>% filter(body != "")
# GENERA WITH BODY, EGG, AND ON DATA
body_egg_ON_genus <- left_join(egg_ON_genus,select(body_genus,c(body,logbody,genus)),by="genus") %>% filter(body != "")

# LINEAGE SPECIFIC DATA
COL1_body_egg_ON_genus <- body_egg_ON_genus %>% filter(group == "Neuropteroidea")
ORT_body_egg_ON_genus <- body_egg_ON_genus %>% filter(group == "Polyneoptera")
HYM_body_egg_ON_genus <- body_egg_ON_genus %>% filter(group == "Hymenoptera")
DIP_body_egg_ON_name <- body_egg_ON_name %>% filter(group == "Antliophora")

### get family level body size data from Rainford et al 2014
family_body_size <- read.delim("analysis/datafiles/family_body_sizes.tsv",stringsAsFactors = F)

# downsample the database to model the effects of the random averaging
### downsample_factor <- 0.50
### The downsample factor is going to be set by the parent R command, when this file is run from source()

# filter the data to entries with egg volume, group by family, and select a random sample
egg_family <- egg_database[sample(nrow(egg_database)),] %>% filter(!(is.na(vol))) %>% group_by(family) %>% sample_frac(downsample_factor)
ON_family <- ON[sample(nrow(ON)),] %>% filter(!(is.na(on))) %>% group_by(family) %>% sample_frac(downsample_factor) 

# reduce to only families with more than 1 observation
### fam_count_threshold <- 1
### The fam_count_threshold is going to be set by the parent R command, when this file is run from source()

egg_fam_count <- egg_family %>% group_by(family) %>% summarize(n = n())
egg_fam_count_n <- egg_fam_count %>% filter(n > fam_count_threshold)
egg_family <- egg_family %>% filter(family %in% egg_fam_count_n$family)

ON_fam_count <- ON_family %>% group_by(family) %>% summarize(n = n())
ON_fam_count_n <- ON_fam_count %>% filter(n > fam_count_threshold)
ON_family <- ON_family %>% filter(family %in% ON_fam_count_n$family)

### log transform the dataframe
family_body_size <- family_body_size %>% mutate(
			logxt = log10(xt),
			logmt = log10(mt),
			body = (xt + mt) / 2,
			logbody = (logxt + logmt) / 2,
			bodyvol = (10^(logbody))^3,
			logbodyvol = log10(bodyvol))

### filter out the order level rows
fam <- family_body_size %>% filter(!(family == ""))

### summarize data by family
egg_family_summary <- egg_family %>% filter(family %in% fam$family) %>%
		group_by(family) %>% 
		summarize_at(vars(vol,ar,logX1,logX2,logar,logvol,X1,X2,asym,curv,sqcurv,sqasym),funs(mean(.,na.rm=T)))

ON_family_summary <- ON_family %>% filter(family %in% fam$family) %>%
		group_by(family) %>% 
		summarize_at(vars(on,logon,sqrton),funs(mean(.,na.rm=T)))

### create full dataset of ovariole number, egg size, and body size by family
body_egg_family <- merge(fam,egg_family_summary,by="family") %>% 
			mutate(egg_ind = vol / bodyvol, 
				logegg_ind = log10(egg_ind), 
				logcub = log10((vol)^(1/3)),
				rank = family)

body_ON_family <- merge(fam,ON_family_summary,by="family") %>% 
			mutate(ind = on / bodyvol,
				rank = family)

body_egg_ON_family <- merge(body_egg_family,ON_family_summary,by="family")


### get order level rows
ordr <- family_body_size %>% filter(family == "") 

### summarize egg data by order for selected orders
egg_order_summary <- egg_family %>% filter(order %in% ordr$order) %>% 
		group_by(order) %>% 
		summarize_at(vars(vol,ar,logX1,logX2,logar,logvol,X1,X2,asym,curv,sqcurv,sqasym),funs(mean(.,na.rm=T)))

ON_order_summary <- ON_family %>% filter(order %in% ordr$order) %>% 
		group_by(order) %>% 
		summarize_at(vars(on,logon,sqrton),funs(mean(.,na.rm=T)))

### combine order level and family level datasets
body_egg_order <- merge(ordr,egg_order_summary,by="order") %>% 
				mutate(egg_ind = vol / bodyvol, 
					logegg_ind = log10(egg_ind), 
					logcub = log10((vol)^(1/3)),
					rank = order)

body_ON_order <- merge(ordr,ON_order_summary,by="order") %>% 
				mutate(ind = on / bodyvol,
					rank = order)

body_egg_full <- rbind(body_egg_family,body_egg_order) # family (and a few orders) summaries of egg size
body_ON_full <- rbind(body_ON_family,body_ON_order) # family (and a few orders) summaries of ON

body_egg_ON_order <- merge(body_egg_full,ON_order_summary,by="order") # combining body-egg with ON for a few orders
body_egg_ON_full <- rbind(body_egg_ON_family,body_egg_ON_order) #family (and a few orders) summaries of all 3

### add in groups
body_egg_ON_full$group <- plyr::mapvalues(body_egg_ON_full$order,groups_orders$order,groups_orders$group,warn_missing=F)


### Read in Gilbert PhD thesis body weight and fecundity data
raw_gilbert <- read.delim("analysis/datafiles/gilbert_thesis_weight_fecundity.tsv",header=T,stringsAsFactors=F)
matched_gilbert <- read.delim(file="analysis/datafiles/gilbert_thesis_weight_matched.tsv",header=T,stringsAsFactors=F) %>% filter(!(problem %in% c("no_taxonomy","no_name")))

# add in taxonomic information from TaxReformer
matched_gilbert <- matched_gilbert %>% 
		select(order,family,subfamily,suborder,superfamily,tribe,updated_genus,updated_species,ID,updated_fullname,ncbi_id) %>% 
		rename(updated_order = order)

# format and transform Gilbert data
gilbert <- raw_gilbert[sample(1:nrow(raw_gilbert)),] %>% # shuffle the data
		filter(family != "Trichogrammatidae") %>% # excluding because lifetime fecundity = 0 or 0.1, typo in dataset?
		mutate(name = paste(genus,species,sep="_"),
			logweight=log10(weight),
			logfecundity = log10(fecundity)) %>% 
		left_join(.,matched_gilbert,by="ID") %>%
		mutate(taxonomy = ifelse((updated_fullname == "" | is.na(updated_fullname)),"original","updated")) %>%
		mutate(genus = ifelse(taxonomy == "updated",updated_genus,genus),
			species = ifelse(taxonomy == "updated",updated_species,species),
			order = ifelse(taxonomy == "updated",updated_order,order)) %>% 
		mutate(name = ifelse(species != "",paste(genus,species,sep="_"),genus)) %>%
		mutate(group = plyr::mapvalues(order,groups_orders$order,groups_orders$group,warn_missing=F)) 


# combine the datasets, following the same approach as above
gilbert_genus <- gilbert %>% group_by(genus) %>% slice(1L) %>% ungroup() # select one representative per genus

gilbert_name <- gilbert %>% group_by(name) %>% slice(1L) %>% ungroup() # select one representative per name

# GENERA WITH BODY, ON, AND EGG DATA
ON_gilbert_genus <- left_join(ON_genus,gilbert_genus %>% select(genus,weight,fecundity,logweight,logfecundity),by="genus")
egg_gilbert_genus <- left_join(egg_genus,gilbert_genus %>% select(genus,weight,fecundity,logweight,logfecundity),by="genus")
egg_ON_gilbert_genus <- left_join(egg_ON_genus,gilbert_genus %>% select(genus,weight,fecundity,logweight,logfecundity),by="genus")

# NAMES WITH BODY, ON, AND EGG DATA
ON_gilbert_name <- left_join(ON_name,gilbert_name %>% select(name,weight,fecundity,logweight,logfecundity),by="name")
egg_gilbert_name <- left_join(egg_name,gilbert_name %>% select(name,weight,fecundity,logweight,logfecundity),by="name")
egg_ON_gilbert_name <- left_join(egg_ON_name,gilbert_name %>% select(name,weight,fecundity,logweight,logfecundity),by="name")


