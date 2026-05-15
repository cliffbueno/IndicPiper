# Test checkIndicPiper

setwd("~/IndicPiper/")

source("IndicPiper.R")

checkIndicPiper(meta = "meta_test.csv",
                genus = "genus_test.csv",
                ind = "genus_habitat_indicators_custom.csv")