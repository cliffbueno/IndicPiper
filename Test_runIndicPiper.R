# Test runIndicPiper()

setwd("~/IndicPiper/")

source("IndicPiper.R")

runIndicPiper(meta = "myMetadataTable.csv.gz",
              genus = "myGenusTable.csv.gz",
              n_multipatt_perm = 100, # number of iterations within the multipatt function
              n_runs = 5, # number of times you want to subsample habitats and run multipatt
              n_per_habitat = 215, # sample size per habitat per run. samples randomly selected
              run_cut = 100, # multipatt mean % runs cutoff
              p_cut = 0.01, # multipatt mean Pfdr cutoff
              IndVal_cut = 0.5,  # multipatt mean IndVal cutoff
              seed = 1, # seed for reproducibility
              output = "genus_habitat_indicators_custom.csv") # output file name
