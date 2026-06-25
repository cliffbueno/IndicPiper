# Test runIndicPiper()

setwd("/scratch/alpine/clbd1748/IndicPiper/")

source("IndicPiper.R")

required_packages <- c(
  "dplyr",
  "tibble",
  "data.table",
  "tidyr",
  "permute",
  "indicspecies",
  "ggplot2",
  "rlang",
  "R.utils",
  "reshape2"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {
  
  if (!(pkg %in% installed)) {
    
    message("Installing package: ", pkg)
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
  
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

message("All packages loaded successfully.")

runIndicPiper(meta = "myMetadataTable.csv.gz",
              genus = "myGenusTable.csv.gz",
              n_multipatt_perm = 100, # number of iterations within the multipatt function
              n_runs = 100, # number of times you want to subsample habitats and run multipatt
              n_per_habitat = 250, # sample size per habitat per run. samples randomly selected
              run_cut = 100, # multipatt mean % runs cutoff
              p_cut = 0.01, # multipatt mean Pfdr cutoff
              IndVal_cut = 0.5,  # multipatt mean IndVal cutoff
              seed = 1, # seed for reproducibility
              output = "genus_habitat_indicators_v2.csv") # output file name
