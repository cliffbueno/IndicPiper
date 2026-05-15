# Test prepIndicPiper()

setwd("~/IndicPiper/")

source("IndicPiper.R")

habitat_list <- c("mammalian gut", "soil", "wastewater", "human oral", "human skin",
                  "seawater", "freshwater sediment", "freshwater water",
                  "marine sediment", "human vaginal", "hydrothermal vent",
                  "phyllosphere", "glacier")

prepIndicPiper(habitat_list = habitat_list,
               combine_soil_rhizo = TRUE,
               combine_freshwater = TRUE,
               combine_glacier_ice = TRUE,
               combine_mammalian_gut = TRUE,
               combine_saliva_oral = TRUE)