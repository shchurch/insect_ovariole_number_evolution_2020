# insect_ovariole_number_evolution_2020
Data and code to reproduce the analyses and figures in the manuscript: **"Repeated loss of variation in insect ovary morphology highlights the role of development in life-history evolution"** published in Proceedings of the Royal Society B, https://doi.org/10.1098/rspb.2021.0150 and BioRxiv, https://doi.org/10.1101/2020.07.07.191940.

### manuscript_supplement.Rmd

Code and text to reproduce the supplementary methods, all figures, and all main text figure panels.

### mansucript.Rmd

This file reads the results of manuscript_supplement.RData and outputs the manuscript in pdf

### manuscript_figures/

All figures used in the manuscript. This is where the output of manuscript_supplement.Rmd prints panels

### phylogeny/

Accession numbers, intermediate files, and final tree files

### analysis/

Code and data called in manuscript_supplement.Rmd, used to analyze ovariole number evolution

### INSTRUCTIONS

```
### set up conda environment

# follow steps in set_up_conda_environment.sh

### within R, execute the following files

R 
source("analysis/get_ovariole_number_data.R")
source("analysis/run_posterior_pgls_body_size.R")
source("analysis/run_posterior_pgls_egg_size_ovariole_number.R")
source("analysis/run_posterior_pgls_fecundity.R")
source("analysis/run_posterior_pgls_subclade.R")

### outside of R, execute the following files

Rscript analysis/run_pgls_simulated.R allometry_body_egg_ON_family 0
Rscript analysis/run_pgls_simulated.R allometry_body_egg_ON_family -1
Rscript analysis/run_pgls_simulated.R allometry_gilbert_genus 0
Rscript analysis/run_pgls_simulated.R allometry_gilbert_genus -1

Rscript analysis/reconstruct_ancestral_ovary_morphology_type.R 

# the following expect 50 cores
Rscript analysis/run_mode_of_oogenesis_model_comparison.R egg
Rscript analysis/run_mode_of_oogenesis_model_comparison.R on
Rscript analysis/run_mode_of_oogenesis_model_comparison.R asym
Rscript analysis/run_mode_of_oogenesis_model_comparison.R curv
Rscript analysis/run_mode_of_oogenesis_model_comparison.R ar

cd analysis/BAMM_results_Jan2021

bamm -c bamm.config
bamm -c egg_bamm.config

### update the following files in the folder

cd ..

cp run_posterior_pgls_body_size.RData pgls_results_Jan2021/run_posterior_pgls_body_size.RData
cp run_posterior_pgls_egg_size_ovariole_number.RData pgls_results_Jan2021/run_posterior_pgls_egg_size_ovariole_number.RData
cp run_posterior_pgls_fecundity.RData pgls_results_Jan2021/run_posterior_pgls_fecundity.RData
cp run_posterior_pgls_subclade.RData pgls_results_Jan2021/run_posterior_pgls_subclade.RData

cp run_pgls_simulated_allometry_body_egg_ON_family_-1_image.RData pgls_results_Jan2021/run_pgls_simulated_allometry_body_egg_ON_family_-1_image.RData
cp run_pgls_simulated_allometry_body_egg_ON_family_0_image.RData pgls_results_Jan2021/run_pgls_simulated_allometry_body_egg_ON_family_0_image.RData
cp run_pgls_simulated_allometry_gilbert_genus_-1_image.RData pgls_results_Jan2021/run_pgls_simulated_allometry_gilbert_genus_-1_image.RData
cp run_pgls_simulated_allometry_gilbert_genus_0_image.RData pgls_results_Jan2021/run_pgls_simulated_allometry_gilbert_genus_0_image.RData

cp on_ovary_model_comparison.RData mode_of_oogenesis_results_June2020/on_ovary_model_comparison.RData
cp asym_ovary_model_comparison.RData mode_of_oogenesis_results_June2020/asym_ovary_model_comparison.RData
cp curv_ovary_model_comparison.RData mode_of_oogenesis_results_June2020/curv_ovary_model_comparison.RData
cp egg_ovary_model_comparison.RData mode_of_oogenesis_results_June2020/egg_ovary_model_comparison.RData
cp ar_ovary_model_comparison.RData mode_of_oogenesis_results_June2020/ar_ovary_model_comparison.RData
cp reconstruct_ancestral_ovary_morphology_type.RData mode_of_oogenesis_results_June2020/reconstruct_ancestral_ovary_morphology_type.RData

### build the manuscript files

# knit the file
#    manuscript_supplement.Rmd 

# check that manuscript_supplement.RData was correctly saved

# then knit the file
#    manuscript.Rmd 

```
