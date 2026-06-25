<<<<<<< HEAD
# This script shows how the starting files for IndicPiper were generated
# Take the Sandpiper database and filter it down to 451568 samples
=======
# This script shows how the starting point for IndicPiper was generated
# Take the Sandpiper database and filter it down to 358209 samples
>>>>>>> 12ba2462e210d25c8fd7e892feac6b1670554427
# The full Sandpiper database metadata and singleM profiles were downloaded from Zenodo
# Used Sandpiper v2.0.1 which screened 913165 metagenomes
# The filtering steps were as follows:
# 1. single cell sequencing
# 2. very vague samples (e.g., "metagenome")
# 3. < 150M bac+arc bases
# 4. “WARNING” flags (in v2.0.0 they added a low_complexity flag too)
# 5. vague freshwater and marine
# 6. misclassified freshwater water/sediment
# 7. misclassified marine water/sediment
# 8. habitats with at least 50 samples
# 9. more vague samples (e.g., "gut" (doesn't say which gut!))

setwd("/scratch/alpine/clbd1748/IndicPiper/")
library(dplyr)
library(data.table)
library(tibble)
`%notin%` <- Negate(`%in%`)
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]

# Make with CheckMisclass.R
freshwaterfilt <- readRDS("freshwaterfilt.rds") # misclassified freshwater water/sediment
marinefilt <- readRDS("marinefilt.rds") # misclassified marine water/sediment
<<<<<<< HEAD

# Load and filter
meta_per_acc <- read.csv("sandpiper2.0.1.per_acc_summary.csv.gz") %>% # full database from Zenodo, n  = 925763
  filter(grepl("metagenome", organism)) %>% # keep metagenomes, n = 839436
  filter(organism != "metagenome") %>% # remove vague samples, n = 738031
=======
meta_per_acc <- read.csv("sandpiper1.1.10.per_acc_summary.csv.gz") %>% # full database from Zenodo, n  = 712661
  filter(grepl("metagenome", organism)) %>% # keep metagenomes, n = 647132
  filter(organism != "metagenome") %>% # remove vague samples, n = 576660
>>>>>>> 12ba2462e210d25c8fd7e892feac6b1670554427
  mutate(organism = gsub(" metagenome", "", organism)) %>% # now remove the metagenome label
  filter(bacterial_archaeal_bases > 150000000) %>% # minimum bac/arc bases, n = 660177
  filter(!grepl("WARNING", warning)) %>% # remove warnings, n = 655437
  filter(low_complexity != "yes") %>% # remove low complexity, n = 626614
  filter(sample %notin% freshwaterfilt$sample) %>% # remove misclass freshwater, n = 620650
  filter(sample %notin% marinefilt$sample) %>% # remove misclass marine 616584
  mutate(Habitat = organism) # change column name

habitat_counts <- meta_per_acc %>%
  group_by(Habitat) %>%
  summarize(n = n()) %>%
  ungroup()
n50 <- habitat_counts %>%
  filter(n >= 50)

meta_per_acc <- meta_per_acc %>%
  filter(Habitat %in% n50$Habitat) # keep habitats with n ≥ 50, n = 614970

# Still lots of vague Habitats that need filtering though
# For example, those labeled "gut" this is not detailed enough!
# "human" also not detailed enough, where on human?
# "feces" also not detailed enough, feces from which animal!?
# Recount and sort
habitat_counts <- meta_per_acc %>%
  group_by(Habitat) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  arrange(desc(n))

meta_per_acc <- meta_per_acc %>%
  filter(Habitat %notin% c("gut", "human", "feces", "sediment", "synthetic", 
                           "mouse", "skin", "viral", "food", "bovine", "clinical", 
                           "terrestrial", "oral", "pig", "manure", 
                           "food production", "fermentation", "vaginal", "canine",
                           "metagenomes", "surface", "mammal", "mixed culture",
                           "upper respiratory tract", "respiratory tract", "milk",
                           "sheep", "bird", "primate", "biofilter",
                           "museum specimen", "fungus", "nervous system",
                           "rat", "factory", "outdoor", "tailings", "water",
                           "snake", "hydrocarbon", "runoff", "fish", "phage",
                           "probiotic", "crustacean", "bat", "ciliate", "meat",
                           "symbiont", "farm", "wildlife", "horse", "leachate", 
                           "clothing", "eye", "solid waste", "cetacean", "rodent",
                           "frog", "sand", "protist")) # 461502

# Check again
habitat_counts <- meta_per_acc %>%
  group_by(Habitat) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  arrange(desc(n))

# We could filter this further but we want users to be able to decide the habitats they want
# We set a minimum cutoff of 50 here so people wouldn't be tempted to attempt IndicPiper with very low sample sizes!
# Save the file. Make it publicly available on Zenodo.
# write.table(meta_per_acc, "Sandpiper_Metadata_Filt_n461502.txt", sep = "\t", row.names = F, quote = FALSE)

# Then, make the accompanying singleM coverage table for these 461502
sm <- data.table::fread("sandpiper2.0.1.gtdb.csv.gz", showProgress = TRUE) # Downloaded from Zenodo
# Note: ~238 million rows, large file!
sm_subset <- sm %>%
  filter(sample %in% meta_per_acc$sample)
length(unique(sm_subset$sample)) # 461502
write.table(sm_subset, "sandpiper2.0.1.gtdb_461502.tsv", sep = "\t", row.names = F, quote = FALSE)

# Then, run singleM summarise to convert to relative abundance at genus level
# singlem summarise --input-taxonomic-profile sandpiper2.0.1.gtdb_461502.tsv \
# --output-species-by-site-relative-abundance Sandpiper_GenusRelAbund_Filt_n461502.csv \
# --output-species-by-site-level genus

# This file actually only had 461060 after the singlem summarise step.
meta <- read.delim("Sandpiper_Metadata_Filt_n461502.txt")
genus <- data.table::fread("Sandpiper_GenusRelAbund_Filt_n461502.csv.gz", showProgress = TRUE) %>%
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
write.table(meta, "Sandpiper_Metadata_Filt_n451568.txt", sep = "\t", row.names = F, quote = FALSE)
fwrite(genus, file = "Sandpiper_Genus_Filt_n451568.csv.gz", sep = ",", quote = FALSE,
       compress = "gzip", nThread = 4)

<<<<<<< HEAD
# Make these files publicly available on Zenodo as the starting point.
# Now users of IndicPiper can download the 2 starting files from Zenodo
# Sandpiper_Metadata_Filt_n451568.txt
# Sandpiper_Genus_Filt_n451568.csv.gz
=======
# Make these files publicly available on FigShare as the starting point.
# Now users of IndicPiper can download the 2 starting files from FigShare
# Sandpiper_Metadata_Filt_n358209.txt
# Sandpiper_Genus_Filt_n358209.csv.gz
>>>>>>> 12ba2462e210d25c8fd7e892feac6b1670554427
