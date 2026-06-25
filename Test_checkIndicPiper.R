# Test checkIndicPiper

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

# Make Figure 2 for manuscript.
custom_order = c("phyllosphere", "soil", "freshwater sediment", "freshwater water",
                 "marine sediment", "marine water", "hydrothermal vent", "wastewater",
                 "mammalian gut", "human oral", "human skin", "human vaginal", "glacier", "non-indicator")
checkIndicPiper(meta = "meta_test.csv",
                genus = "genus_test.csv",
                ind = "genus_habitat_indicators_v2.csv",
                ouput = "Figure2.pdf",
                ncol = 2,
                keep_unassigned = FALSE,
                keep_nonindicator = TRUE,
                custom_order = custom_order)
