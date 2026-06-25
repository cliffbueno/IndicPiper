# Test countHabitats()

setwd("/scratch/alpine/clbd1748/IndicPiper")

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
  "R.utils"
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

countHabitats(meta = "Sandpiper_Metadata_Filt_n451568.txt")
