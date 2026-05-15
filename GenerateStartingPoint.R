# This script shows how the starting point for IndicPiper was generated
# Take the Sandpiper database and filter it down to 483k samples
# The full Sandpiper database metadata and singleM profiles were downloaded from Zenodo
# Used Sandpiper v1.1.10 which screened 707470 metagenomes
# The filtering steps were as follows:
# 1. single cell sequencing
# 2. very vague samples (e.g., "metagenome")
# 3. < 150M bac+arc bases
# 4. “WARNING” flags
# 5. vague freshwater and marine
# 6. misclassified freshwater water/sediment
# 7. misclassified marine water/sediment
# 8. habitats with at least 50 samples
# 9. more vague samples (e.g., "gut" (doesn't say which gut!), "canine" (doesn't say which part!))

setwd("~/IndicPiper/")
library(dplyr)
library(data.table)
library(tibble)
`%notin%` <- Negate(`%in%`)

freshwaterfilt <- readRDS("freshwaterfilt.rds") # misclassified freshwater water/sediment
marinefilt <- readRDS("marinefilt.rds") # misclassified marine water/sediment
meta_per_acc <- read.csv("sandpiper1.1.10.per_acc_summary.csv.gz") %>% # full database, n  = 712661
  filter(grepl("metagenome", organism)) %>% # keep metagenomes, n = 647132
  filter(organism != "metagenome") %>% # remove vague samples, n = 576660
  mutate(organism = gsub(" metagenome", "", organism)) %>% # now remove the metagenome label
  filter(bacterial_archaeal_bases > 150000000) %>% # minimum bac/arc bases, n = 509425
  filter(!grepl("WARNING", warning)) %>% # remove warnings, n = 506252
  filter(organism != "freshwater") %>% # remove vague freshwater, n = 499755
  filter(organism != "marine") %>% # # remove vague marine, n = 486726
  filter(sample %notin% freshwaterfilt$sample) %>% # remove misclass freshwater, n = 484276
  filter(sample %notin% marinefilt$sample) %>% # remove misclass marine 483309
  mutate(Habitat = organism) # change column name

habitat_counts <- meta_per_acc %>%
  group_by(Habitat) %>%
  summarize(n = n()) %>%
  ungroup()
n50 <- habitat_counts %>%
  filter(n >= 50)

meta_per_acc <- meta_per_acc %>%
  filter(Habitat %in% n50$Habitat) # keep habitats with n ≥ 50, n = 481450

# Still lots of vague Habitats that need filtering though
# For example, those labeled "gut" this is not detailed enough!
# "human" also not detailed enough, where on human?
# "feces" also not detailed enough, from which animal!?
habitat_counts <- meta_per_acc %>%
  group_by(Habitat) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  arrange(desc(n))

meta_per_acc <- meta_per_acc %>%
  filter(Habitat %notin% c("gut", "human", "feces", "sediment", "food",
                           "synthetic", "skin", "clinical", "terrestrial",
                           "vaginal", "oral", "food production", "mixed culture",
                           "factory", "surface", "museum specimen", "metagenomes",
                           "biofilter", "symbiont", "outdoor", "whole organism",
                           "wildlife", "clothing")) # 366538

# We could filter this further but we want users to be able to decide the habitats they want
# We set a minimum cutoff of 50 here so people wouldn't be tempted to attempt IndicPiper with very low sample sizes!
# Save the file. Make it publicly available on FigShare.
# write.table(meta_per_acc, "Sandpiper_Metadata_Filt_n366538.txt", sep = "\t", row.names = F, quote = FALSE)

# Then, make the accompanying singleM coverage table for these 366538
sm <- data.table::fread("sandpiper1.1.0.gtdb.csv.gz", showProgress = TRUE) # Downloaded from Zenodo
# Note: ~167 million rows, large file!
sm_subset <- sm %>%
  filter(sample %in% meta_per_acc$sample)
length(unique(sm_subset$sample)) # 366538
write.table(sm_subset, "sandpiper1.1.0.gtdb_366538.tsv", sep = "\t", row.names = F, quote = FALSE)

# Then, run singleM summarise to convert to relative abundance at genus level
# singlem summarise --input-taxonomic-profile sandpiper1.1.0.gtdb_366538.tsv \
# --output-species-by-site-relative-abundance Sandpiper_GenusRelAbund_Filt_n366538.csv \
# --output-species-by-site-level genus

# This file actually only had 365870 after the singlem summarise step.
meta <- read.delim("Sandpiper_Metadata_Filt_n366538.txt")
genus <- data.table::fread("Sandpiper_GenusRelAbund_Filt_n366538.csv.gz", showProgress = TRUE) %>%
  column_to_rownames(var = "taxonomy")
meta <- meta %>%
  filter(sample %in% colnames(genus))
sum(meta$sample != colnames(genus))
which(rownames(genus) == "unassigned") # unassigned is row 1

# Need to filter out samples with > 60% unassigned anyways, so do that here then save
# Remove samples with < 60% genera unassigned
meta$unassigned <- as.numeric(genus[1,])
meta <- meta %>%
  filter(unassigned < 60)
genus <- genus %>%
  dplyr::select(all_of(meta$sample))
sum(meta$sample != colnames(genus))
genus <- genus %>%
  rownames_to_column(var = "taxonomy")
write.table(meta, "Sandpiper_Metadata_Filt_n358209.txt", sep = "\t", row.names = F, quote = FALSE)
fwrite(genus, file = "Sandpiper_Genus_Filt_n358209.csv.gz", sep = ",", quote = FALSE,
       compress = "gzip", nThread = 4)

# Make these files publicly available on FigShare as the starting point.
# Now users of IndicPiper can download the 2 starting files from FigShare
# Sandpiper_Metadata_Filt_n358209.txt
# Sandpiper_Genus_Filt_n358209.csv.gz