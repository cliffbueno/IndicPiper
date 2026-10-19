# Filter the Sandpiper database
# Generate IndicPiper "starting files"
# These are the two starting files people can use to run IndicPiper
# This script is the updated version of the old "GenerateStartingPoint.R" script
# This approach is more streamlined and robust
# Implements feedback from ISME Comms peer reviewers on the filtering and subsetting
# Also removes repeat BioSample runs

#### 1. Setup ####
library(dplyr)
library(data.table)
library(tibble)
library(vegan)
library(parallel)
`%notin%` <- Negate(`%in%`)
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
setwd("/scratch/alpine/clbd1748/IndicPiperv2.1/")



#### 2. Prelim Metadata Filter ####
meta_kf <- data.table::fread("kingfisher_metadata_subset.tsv.gz") %>% # from Zenodo, n = 991856
  dplyr::select(-latitude, -longitude, -location, -geographic_location)

meta_per_acc <- read.csv("sandpiper2.0.1.per_acc_summary.csv.gz") %>% # full database from Zenodo, n = 925763
  filter(grepl("metagenome", organism)) %>% # keep metagenomes, n = 839436
  filter(organism != "metagenome") %>% # remove vague samples, n = 738031
  mutate(organism = gsub(" metagenome", "", organism)) %>% # now remove the metagenome label
  filter(bacterial_archaeal_bases > 150000000) %>% # minimum bac/arc bases, n = 660177
  filter(!grepl("WARNING", warning)) %>% # remove warnings, n = 655437
  filter(low_complexity != "yes") %>% # remove low complexity, n = 626614
  mutate(Habitat = organism) %>% # change column name
  left_join(., meta_kf, by = c("sample" = "run")) %>% # join Kingfisher data
  filter(!is.na(biosample)) %>% # remove samples without Kingfisher data, n = 623134
  filter(!duplicated(biosample)) # remove duplicates, n = 572811

habitat_counts <- meta_per_acc %>%
  group_by(Habitat) %>%
  summarize(n = n()) %>%
  ungroup()
n56 <- habitat_counts %>%
  filter(n >= 56) # use 56 instead of 50 since after more filtering will be down to 50.

meta_per_acc <- meta_per_acc %>%
  filter(Habitat %in% n56$Habitat) # keep habitats with n ≥ 56, n = 570794

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
                           "frog", "sand", "protist")) # 437089

# Check again
habitat_counts <- meta_per_acc %>%
  group_by(Habitat) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  arrange(desc(n))

# We could filter this further but we want users to be able to decide the habitats they want
# We set a minimum cutoff of 56 here so people wouldn't be tempted to attempt IndicPiper with very low sample sizes!
# Save this file
data.table::fwrite(meta_per_acc, "Sandpiper_Metadata_Filt_n437089.tsv", sep = "\t")
# We still have 2 more filtering steps to go though
# For the next 2 steps, we need the taxonomic profiles though
# 1. Filter samples with > 60% unassigned
# 2. Filter samples with extreme within-habitat Bray-Curtis dissimilarities


# Then, make the accompanying SingleM coverage table for these 437089
sm <- data.table::fread("sandpiper2.0.1.gtdb.csv.gz", 
                        showProgress = TRUE,
                        nThread = 8) # Downloaded from Zenodo
# Note: ~238 million rows, large file!
sm_subset <- sm %>%
  filter(sample %in% meta_per_acc$sample)
length(unique(sm_subset$sample)) # 437089
write.table(sm_subset, "sandpiper2.0.1.gtdb_437089.tsv", sep = "\t", row.names = F, quote = FALSE)

# Then, run SingleM summarise to convert to relative abundance at genus level
# singlem summarise --input-taxonomic-profile sandpiper2.0.1.gtdb_437089.tsv \
# --output-species-by-site-relative-abundance Sandpiper_GenusRelAbund_Filt_n437089.csv \
# --output-species-by-site-level genus



#### 3. Genus unassigned filter ####
# This file actually only had 436655 after the singlem summarise step.
meta <- data.table::fread("Sandpiper_Metadata_Filt_n437089.tsv", data.table = FALSE)
genus <- data.table::fread("Sandpiper_GenusRelAbund_Filt_n437089.csv", 
                           showProgress = TRUE,
                           nThread = 8) %>%
  column_to_rownames(var = "taxonomy")
meta <- meta %>%
  filter(sample %in% colnames(genus))
sum(meta$sample != colnames(genus))
which(rownames(genus) == "unassigned") # unassigned is row 1
print("Dim meta")
dim(meta)
print("Dim genus")
dim(genus)

# Need to filter out samples with > 60% unassigned (keep samples with < 60% genera unassigned)
meta$unassigned <- as.numeric(genus[1,])
meta <- meta %>%
  filter(unassigned < 60)
genus <- genus %>%
  dplyr::select(all_of(meta$sample))
sum(meta$sample != colnames(genus))
genus <- genus %>%
  rownames_to_column(var = "taxonomy")
print("Dim meta")
dim(meta)
print("Dim genus")
dim(genus)
data.table::fwrite(meta, file = "Sandpiper_Metadata_Filt_UA.tsv", sep = "\t")
data.table::fwrite(genus, file = "Sandpiper_Genus_Filt_UA.csv.gz", sep = ",", 
                   quote = FALSE, compress = "gzip", nThread = 4)
# After this step we are now at 428,286 samples



#### 4. Outliers filter ####
# Added this during manuscript revisions
# For each habitat calculate within-habitat Bray-Curtis dissimilarity
# Remove samples with the highest mean Bray-Curtis dissimilarity (> 3 SD)
# Just want to filter out outliers that are likely misclassified
# Ran this as it's own script, but including here so it's part of the broader pipeline.
library(dplyr)
library(data.table)
library(tibble)
library(vegan)
library(parallel)
library(ggplot2)
options(bitmapType = "cairo")

`%notin%` <- Negate(`%in%`)
setwd("/scratch/alpine/clbd1748/IndicPiperv2.1/")

# Load data
meta <- data.table::fread("Sandpiper_Metadata_Filt_UA.tsv", data.table = FALSE)
genus <- data.table::fread("Sandpiper_Genus_Filt_UA.csv.gz",
                           showProgress = TRUE,
                           nThread = 8) %>%
  column_to_rownames(var = "taxonomy")

# Directories
checkpoint_dir <- "checkpoints"
log_file       <- "bc_filtering.log"
plot_dir       <- "."
dir.create(checkpoint_dir, showWarnings = FALSE)
dir.create("histograms",   showWarnings = FALSE)
dir.create("boxplots",     showWarnings = FALSE)

# Logging function
log_msg <- function(msg, log_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line      <- paste0("[", timestamp, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

# Inner BC filter and plot function
run_bc_filter <- function(samples, sd_threshold, plot_suffix = "",
                          habitat_name, hab_clean, hab_genus,
                          plot_dir, log_file) {
  mat <- hab_genus[, samples, drop = FALSE]
  mat <- mat[rowSums(mat) > 0, ]
  bc_mat <- as.matrix(vegdist(t(mat), method = "bray"))
  diag(bc_mat) <- NA
  mean_bc   <- rowMeans(bc_mat, na.rm = TRUE)
  threshold <- mean(mean_bc) + sd_threshold * sd(mean_bc)
  
  rm(bc_mat, mat); gc()
  
  plot_df   <- data.frame(sampleID = names(mean_bc), mean_bc = mean_bc)
  file_base <- paste0(hab_clean, plot_suffix)
  
  # Histogram
  pdf(file.path(plot_dir, "histograms", paste0(file_base, "_histogram.pdf")),
      width = 7, height = 5)
  print(
    ggplot(plot_df, aes(x = mean_bc)) +
      geom_histogram(bins = 30, fill = "grey70", color = "white") +
      geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 0.8) +
      labs(title = paste0(habitat_name, plot_suffix),
           x = "Mean Bray-Curtis dissimilarity", y = "Count") +
      theme_bw()
  )
  dev.off()
  
  # Boxplot
  pdf(file.path(plot_dir, "boxplots", paste0(file_base, "_boxplot.pdf")),
      width = 7, height = 5)
  print(
    ggplot(plot_df, aes(x = habitat_name, y = mean_bc)) +
      geom_boxplot(outlier.color = "red", outlier.size = 3) +
      geom_jitter(pch = 16, width = 0.25, size = 2, alpha = 0.25) +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed", linewidth = 0.8) +
      labs(title = paste0(habitat_name, plot_suffix),
           x = NULL, y = "Mean Bray-Curtis dissimilarity") +
      theme_bw()
  )
  dev.off()
  
  rm(plot_df); gc()
  
  kept_samples     <- names(mean_bc[mean_bc <= threshold])
  filtered_samples <- names(mean_bc[mean_bc >  threshold])
  
  log_msg(paste0("  Chunk '", hab_clean, plot_suffix, "' complete — ",
                 length(kept_samples), " kept, ",
                 length(filtered_samples), " filtered"), log_file)
  
  list(kept = kept_samples, filtered = filtered_samples)
}

# Main filtering function
filter_bc_outliers <- function(habitat_name, meta, genus,
                               sd_threshold   = 3,
                               chunk_size     = 2000,
                               large_thresh   = 4000,
                               seed           = 42,
                               plot_dir       = ".",
                               checkpoint_dir = "checkpoints",
                               log_file       = "bc_filtering.log") {
  
  hab_clean       <- gsub("[^a-zA-Z0-9]", "_", habitat_name)
  checkpoint_file <- file.path(checkpoint_dir, paste0(hab_clean, ".rds"))
  
  # Resume: skip if already done
  if (file.exists(checkpoint_file)) {
    log_msg(paste0("SKIPPING '", habitat_name, "' — checkpoint found"), log_file)
    return(readRDS(checkpoint_file))
  }
  
  log_msg(paste0("STARTING '", habitat_name, "'"), log_file)
  
  hab_meta <- meta %>% filter(Habitat == habitat_name)
  n        <- nrow(hab_meta)
  
  if (n < 3) {
    result <- list(habitat = habitat_name, kept = hab_meta$sample, filtered = character(0))
    saveRDS(result, checkpoint_file)
    log_msg(paste0("DONE '", habitat_name, "' — fewer than 3 samples, all kept"), log_file)
    return(result)
  }
  
  hab_genus <- genus %>%
    dplyr::select(all_of(hab_meta$sample)) %>%
    .[rowSums(.) > 0, ]
  
  all_samples <- colnames(hab_genus)
  
  if (n <= large_thresh) {
    # Small habitat: full BC on all samples
    res <- run_bc_filter(all_samples, sd_threshold,
                         plot_suffix  = "",
                         habitat_name = habitat_name,
                         hab_clean    = hab_clean,
                         hab_genus    = hab_genus,
                         plot_dir     = plot_dir,
                         log_file     = log_file)
    result <- list(habitat  = habitat_name,
                   kept     = res$kept,
                   filtered = res$filtered)
    
  } else {
    # Large habitat: shuffle then partition into equal non-overlapping chunks
    n_chunks <- ceiling(n / chunk_size)
    set.seed(seed)
    shuffled <- sample(all_samples)
    chunks   <- split(shuffled, cut(seq_along(shuffled), breaks = n_chunks, labels = FALSE))
    
    results <- lapply(seq_along(chunks), function(i) {
      run_bc_filter(chunks[[i]], sd_threshold,
                    plot_suffix  = paste0("_", i),
                    habitat_name = habitat_name,
                    hab_clean    = hab_clean,
                    hab_genus    = hab_genus,
                    plot_dir     = plot_dir,
                    log_file     = log_file)
    })
    
    result <- list(
      habitat  = habitat_name,
      kept     = unlist(lapply(results, `[[`, "kept"),     use.names = FALSE),
      filtered = unlist(lapply(results, `[[`, "filtered"), use.names = FALSE)
    )
    rm(results, chunks, shuffled); gc()
  }
  
  # Save checkpoint
  saveRDS(result, checkpoint_file)
  
  log_msg(paste0("DONE '", habitat_name, "' — ",
                 length(result$kept), " kept, ",
                 length(result$filtered), " filtered (",
                 round(100 * length(result$filtered) / n, 1), "% removed)"), log_file)
  
  rm(hab_genus, hab_meta); gc()
  
  result
}

# Habitat tiers
habitat_counts <- meta %>% count(Habitat, name = "n_samples")
small_habs <- habitat_counts %>% filter(n_samples <= 4000) %>% pull(Habitat)
large_habs <- habitat_counts %>% filter(n_samples >  4000) %>% pull(Habitat)

log_msg("=== BC outlier filtering started ===", log_file)
log_msg(paste0("Small habitats: ", length(small_habs),
               " | Large habitats: ", length(large_habs)), log_file)

# --- Small habitats: parallel ---
log_msg("--- Starting small habitats ---", log_file)
results_small <- mclapply(small_habs, function(hab) {
  tryCatch(
    filter_bc_outliers(hab, meta, genus,
                       sd_threshold   = 3,
                       chunk_size     = 2000,
                       large_thresh   = 4000,
                       plot_dir       = plot_dir,
                       checkpoint_dir = checkpoint_dir,
                       log_file       = log_file),
    error = function(e) {
      log_msg(paste0("ERROR in '", hab, "': ", e$message), log_file)
      NULL
    }
  )
}, mc.cores = 32)
names(results_small) <- small_habs

failed_small <- small_habs[sapply(results_small, is.null)]
if (length(failed_small) > 0) {
  log_msg(paste0("WARNING: ", length(failed_small),
                 " small habitats failed: ",
                 paste(failed_small, collapse = ", ")), log_file)
}
results_small <- Filter(Negate(is.null), results_small)
log_msg(paste0("--- Small habitats complete: ",
               length(results_small), "/", length(small_habs), " succeeded ---"), log_file)

# --- Large habitats: sequential, chunked internally ---
log_msg("--- Starting large habitats ---", log_file)
results_large <- lapply(large_habs, function(hab) {
  tryCatch(
    filter_bc_outliers(hab, meta, genus,
                       sd_threshold   = 3,
                       chunk_size     = 2000,
                       large_thresh   = 4000,
                       plot_dir       = plot_dir,
                       checkpoint_dir = checkpoint_dir,
                       log_file       = log_file),
    error = function(e) {
      log_msg(paste0("ERROR in '", hab, "': ", e$message), log_file)
      NULL
    }
  )
})
names(results_large) <- large_habs

failed_large <- large_habs[sapply(results_large, is.null)]
if (length(failed_large) > 0) {
  log_msg(paste0("WARNING: ", length(failed_large),
                 " large habitats failed: ",
                 paste(failed_large, collapse = ", ")), log_file)
}
results_large <- Filter(Negate(is.null), results_large)
log_msg(paste0("--- Large habitats complete: ",
               length(results_large), "/", length(large_habs), " succeeded ---"), log_file)

# --- Combine results ---
all_results <- c(results_small, results_large)

summary_df <- do.call(rbind, lapply(all_results, function(res) {
  n_total <- length(res$kept) + length(res$filtered)
  data.frame(
    Habitat      = res$habitat,
    n_kept       = length(res$kept),
    n_filtered   = length(res$filtered),
    pct_filtered = round(100 * length(res$filtered) / n_total, 2)
  )
}))

data.table::fwrite(summary_df, "bc_filtering_summary.tsv", sep = "\t")
log_msg("Summary table written", log_file)

all_kept     <- unlist(lapply(all_results, `[[`, "kept"),     use.names = FALSE)
all_filtered <- unlist(lapply(all_results, `[[`, "filtered"), use.names = FALSE)

log_msg(paste0("=== Filtering complete — ", length(all_kept), " kept, ",
               length(all_filtered), " filtered across ",
               length(all_results), " habitats ==="), log_file)

# --- Write outputs ---
log_msg("--- Writing filtered outputs ---", log_file)

meta_filtered <- meta %>% filter(sample %in% all_kept)
data.table::fwrite(meta_filtered, "Sandpiper_Metadata_Filt_n426084.txt", sep = "\t")
log_msg("Metadata written", log_file)

genus_filtered <- genus %>%
  dplyr::select(all_of(meta$sample)) %>%
  rownames_to_column(var = "taxonomy")
data.table::fwrite(genus_filtered, "Sandpiper_Genus_Filt_n426084.csv.gz", nThread = 8)
log_msg("Genus table written", log_file)

log_msg("=== ALL DONE ===", log_file)



#### 5. Save ####
# Make these files publicly available on Zenodo as the starting point.
# Users of IndicPiper should download the 2 starting files from Zenodo
# Sandpiper_Metadata_Filt_n426084.txt
# Sandpiper_Genus_Filt_n426084.csv.gz