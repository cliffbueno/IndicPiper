# Test prepIndicPiper()

setwd("/scratch/alpine/clbd1748/IndicPiper/")

source("IndicPiper.R")

required_packages <- c(
  "dplyr",
  "tibble",
  "data.table",
  "tidyr",
  "permute",
  "indicspecies",
  "FSA",
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

habitat_list <- c("mammalian gut", "soil", "wastewater", "human oral", "human skin",
                  "seawater", "freshwater sediment", "freshwater water",
                  "marine sediment", "human vaginal", "hydrothermal vent",
                  "phyllosphere", "glacier", "peat", "groundwater")

prepIndicPiper(habitat_list = habitat_list,
               combine_soil_rhizo = TRUE,
               combine_freshwater = TRUE,
               combine_glacier_ice = TRUE,
               combine_mammalian_gut = TRUE,
               combine_saliva_oral = TRUE)
