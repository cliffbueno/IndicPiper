# IndicPiper Main Functions and Support Functions

# First, download the 2 starting files from Zenodo
# Sandpiper_Metadata_Filt_n451568.txt
# Sandpiper_Genus_Filt_n451568.csv.gz

# Put them in your working directory



#### countHabitats ####
countHabitats <- function(meta = "Sandpiper_Metadata_Filt_n451568.txt") {
  logtime("Starting countHabitats()")
  
  # Load libraries
  suppressMessages(library(dplyr))
  suppressMessages(library(tibble))
  suppressMessages(library(data.table))
  logtime("Libraries loaded")
  
  # Load metadata
  meta <- read.delim(meta, stringsAsFactors = FALSE)
  logtime("Metadata table loaded")
  
  habitat_counts <- meta %>%
    group_by(Habitat) %>%
    summarize(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  print(habitat_counts, n = Inf)
  flush.console()
  logtime("Habitats counted")
}

#### prepIndicPiper ####
# Prepare data for IndicPiper. Then you can decide n_per_habitat for runIndicPiper()
prepIndicPiper <- function(habitat_list,
                           combine_soil_rhizo = TRUE,
                           combine_freshwater = TRUE,
                           combine_glacier_ice = TRUE,
                           combine_mammalian_gut = TRUE,
                           combine_saliva_oral = TRUE)
{
  
  logtime("Starting prepIndicPiper")
  
  # Load libraries
  suppressMessages(library(dplyr))
  suppressMessages(library(tibble))
  suppressMessages(library(data.table))
  logtime("Libraries loaded")
  
  # Load metadata
  meta <- read.delim("Sandpiper_Metadata_Filt_n451568.txt", stringsAsFactors = FALSE)
  logtime("Metadata table loaded")
  
  # Change marine to seawater
  meta <- meta %>%
    mutate(Habitat = ifelse(Habitat == "marine", "seawater", Habitat))
  
  # Print starting habitat counts
  habitat_counts <- meta %>%
    group_by(Habitat) %>%
    summarize(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  print(habitat_counts, n = Inf)
  flush.console()
  logtime("Habitats counted")
  
  # Load genus table
  genus <- fread("Sandpiper_Genus_Filt_n451568.csv.gz", showProgress = TRUE, nThread = 2)
  logtime("Raw genus table loaded")
  # Move taxonomy column to rownames
  rownames(genus) <- genus[[1]]
  genus[[1]] <- NULL
  logtime("Taxonomy column converted to rownames")
  
  # Transpose genus table
  logtime("Starting transpose")
  genus_table_t <- transpose(genus)
  rownames(genus_table_t) <- colnames(genus)
  colnames(genus_table_t) <- rownames(genus)
  sum(meta$sample != rownames(genus_table_t))
  rm(genus)
  logtime("Transpose complete")
  
  # Combine habitats
  if (combine_soil_rhizo) {
    meta$Habitat[meta$Habitat %in% "rhizosphere"] <- "soil"
    logtime("Soil and rhizosphere combined")
  }
  
  if (combine_freshwater) {
    meta$Habitat[
      meta$Habitat %in% c("freshwater", "aquatic", "pond", "lake water", "riverine")
    ] <- "freshwater water"
    logtime("Freshwater habitats combined")
  }
  
  if (combine_glacier_ice) {
    meta$Habitat[meta$Habitat %in% "ice"] <- "glacier"
    logtime("Glacier and ice combined")
  }
  
  if (combine_mammalian_gut) {
    meta$Habitat[
      meta$Habitat %in% c("human gut", "pig gut", "bovine gut", "sheep gut",
        "canine gut", "goat gut"
      )
    ] <- "mammalian gut"
    logtime("Mammalian guts combined")
  }
  
  if (combine_saliva_oral) {
    meta$Habitat[meta$Habitat %in% "human saliva"] <- "human oral"
    logtime("Human saliva and oral combined")
  }
  
  # Filter metadata
  meta <- meta %>%
    filter(Habitat %in% habitat_list)
  logtime("Metadata table filtered to habitats of interest")
  cat("Remaining metadata rows:", nrow(meta), "\n")
  flush.console()
  
  # Filter genus table to samples
  genus_table_t <- genus_table_t %>% 
    filter(rownames(.) %in% meta$sample)
  if (sum(meta$sample != rownames(genus_table_t)) > 0) { 
    stop("Metadata table and genus table sample order do not match.") } 
  logtime("Genus table filtered to habitats of interest")
  
  # Remove zero-prevalence taxa
  logtime("Starting prevalence filtering")
  cat("Genus table dimensions before filtering:\n")
  print(dim(genus_table_t))
  flush.console()
  
  genus_prevalence <- data.table(Genus = colnames(genus_table_t), 
                                 Prevalence = colSums(genus_table_t > 0))
  tokeep <- genus_prevalence %>% 
    filter(Prevalence > 0) %>% # Very lenient, just get rid of zeroes. 
    filter(Genus != "unassigned") 
  genus_table_t <- genus_table_t %>% 
    dplyr::select(all_of(tokeep$Genus))
  
  logtime("Prevalence filtering complete")
  cat("Genus table dimensions after filtering:\n")
  print(dim(genus_table_t))
  flush.console()
  
  # Save outputs
  meta <- meta %>%
    mutate(Habitat = gsub("seawater", "marine water", Habitat))
  logtime("Changed 'seawater' to 'marine water'")
  logtime("Starting file writes")
  fwrite(meta, "myMetadataTable.csv.gz", compress = "gzip", nThread = 2, row.names = FALSE)
  genus_table_t$sampleID <- rownames(genus_table_t)
  fwrite(genus_table_t, "myGenusTable.csv.gz", compress = "gzip", nThread = 2, row.names = FALSE)
  logtime("Metadata table and genus table saved")
  
  # Final habitat counts
  habitat_counts2 <- meta %>%
    group_by(Habitat) %>%
    summarize(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  print(habitat_counts2, n = Inf)
  flush.console()
  logtime("Habitats recounted")
  cat("prepIndicPiper complete\n")
  flush.console()
}



#### runIndicPiper #####
# Run IndicPiper using the data you prepped with prepIndicPiper
runIndicPiper <- function(meta = "myMetadataTable.csv.gz",
                          genus = "myGenusTable.csv.gz",
                          n_multipatt_perm = 100, # number of iterations within the multipatt function
                          n_runs = 100, # number of times you want to subsample habitats and run multipatt
                          n_per_habitat = 250, # sample size per habitat per run. samples randomly selected
                          run_cut = 100, # multipatt mean % runs cutoff
                          p_cut = 0.01, # multipatt mean Pfdr cutoff
                          IndVal_cut = 0.5,  # multipatt mean IndVal cutoff
                          seed = 1, # seed for reproducibility
                          output = "genus_habitat_indicators_custom.csv") # output file name
{
  logtime("Starting runIndicPiper()")
  
  # Load libraries
  suppressMessages(library(tidyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(tibble))
  suppressMessages(library(data.table))
  suppressMessages(library(indicspecies))
  logtime("Libraries loaded")
  
  # Import data produced by prepIndicPiper.R
  meta <- fread(meta, showProgress = TRUE, nThread = 2)
  n_habitats <- length(levels(as.factor(meta$Habitat)))
  genus_table_t_start <- fread(genus, showProgress = TRUE, nThread = 16) %>%
    column_to_rownames(var = "sampleID")
  if (sum(meta$sample != rownames(genus_table_t_start)) > 0) { 
    stop("Metadata table and genus table sample order do not match.") } 
  logtime("Data imported")
  
  # Make subsets, given the n_runs and n_per_habitat you want
  subs <- make_subsamples(meta, n_runs = n_runs, n_per_habitat = n_per_habitat, seed = seed)
  saveRDS(subs, "subs.rds")
  logtime("Subsets for multipatt made and saved as subs.rds")
  
  # Run multipatt analysis given the n_multipatt_perm you want
  mp <- vector("list", n_runs)
  for (i in 1:n_runs) {
    
    # Prep the genus table for the run (subset to the samples and filter zeroes)
    genus_table_t <- genus_table_t_start %>%
      filter(rownames(.) %in% subs[[i]]$sample) %>%
      arrange(match(rownames(.), subs[[i]]$sample))
    genus_prevalence <- data.table(Genus = colnames(genus_table_t),
                                   Prevalence = colSums(genus_table_t > 0)) %>%
      arrange(Prevalence)
    tokeep <- genus_prevalence %>%
      filter(Prevalence > 0) %>% # Very lenient, just get rid of zeroes. 
      filter(Genus != "unassigned")
    genus_table_t <- genus_table_t %>%
      dplyr::select(all_of(tokeep$Genus))
    sum(subs[[i]]$sample != rownames(genus_table_t))
    
    # Now, run multipatt
    set.seed(seed)
    mp[[i]] <- multipatt(x = genus_table_t, 
                         cluster = subs[[i]]$Habitat, 
                         func = "IndVal.g", 
                         duleg = TRUE,
                         control = how(nperm = n_multipatt_perm))
    logtime(paste0("Finished iteration ", i))
  }
  
  logtime(paste0("Finshed all ", n_runs, " iterations"))
  
  # Save one iteration (in this case the last one) of meta and genus
  write.csv(subs[[n_runs]], "meta_test.csv")
  write.csv(genus_table_t, "genus_test.csv")
  logtime("Saved one iteration for checkIndicPiper()")
  
  # Merge the filter results of the n_runs of multipatt
  s <- vector("list", n_runs)
  for (i in 1:n_runs) {
    s[[i]] <- mp[[i]]$sign %>%
      rownames_to_column(var = "Genus") %>%
      dplyr::select(-index) %>%
      rename(IndVal = stat) %>%
      pivot_longer(cols = c(2:(n_habitats+1)), values_to = "value", names_to = "Habitat") %>%
      filter(value == 1) %>%
      dplyr::select(-value) %>%
      mutate(Habitat = gsub("s\\.", "", Habitat)) %>%
      mutate(Genus = gsub("Root; |d__|p__|c__|o__|f__|g__", "", Genus)) %>%
      rename(Taxonomy = Genus) %>%
      mutate(Run = i)
  }
  
  comb <- s[[1]]
  for (i in 2:n_runs) {
    comb <- rbind(comb, s[[i]])
  }
  taxonomy_counts <- as.data.frame(table(comb$Taxonomy)) %>%
    mutate(Perc = Freq / n_runs * 100)
  tax100 <- taxonomy_counts %>%
    filter(Perc == 100)
  result <- comb %>%
    filter(Taxonomy %in% tax100$Var1) %>% # Genera present in all runs
    group_by(Taxonomy, Habitat) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Taxonomy) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>% # Most frequent habitat
    rename(nRunsInd = n) %>%
    mutate(PercRunsInd = nRunsInd / n_runs * 100) %>%
    left_join(., comb, by = c("Taxonomy", "Habitat")) %>%
    group_by(Taxonomy, Habitat, PercRunsInd) %>%
    summarise(IndVal_mean = mean(IndVal),
              IndVal_se = se(IndVal),
              IndVal_min = min(IndVal),
              IndVal_max = max(IndVal),
              PVal_mean = mean(p.value),
              PVal_se = se(p.value),
              PVal_min = min(p.value),
              PVal_max = max(p.value),
              .groups = "drop") %>%
    filter(PercRunsInd >= run_cut & IndVal_mean >= IndVal_cut & PVal_mean < p_cut)
  logtime("Finished merging and filtering multipatt output")
  
  # Check numbrer of indicators per habitat
  table(result$Habitat)
  
  # Save the taxa list - these are your habitat indicators!
  write.csv(result, output, row.names = FALSE)
  logtime("Saved genus habitat indicator list. runIndicPiper() complete")
}



#### checkIndicPiper ####
checkIndicPiper <- function(meta = "meta_test.csv",
                            genus = "genus_test.csv",
                            ind = "genus_habitat_indicators_custom.csv",
                            ouput = "IndicPiper_SummedAbund.pdf",
                            ncol = 2,
                            keep_unassigned = FALSE,
                            keep_nonindicator = TRUE,
                            custom_order = NULL) {
  logtime("Starting checkIndicPiper()")
  
  # Load libraries
  suppressMessages(library(data.table))
  suppressMessages(library(tidyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(tibble))
  suppressMessages(library(ggplot2))
  suppressMessages(library(rlang))
  logtime("Libraries loaded")
  
  # Import data
  meta_test <- read.csv(meta, row.names = 1)
  genus_test <- fread(genus, showProgress = TRUE, nThread = 2) %>%
    column_to_rownames(var = "V1")
  ind_merge <- read.csv(ind) %>%
    dplyr::select(Taxonomy, Habitat)
  logtime("Test data imported")
  
  cs1 <- data.frame(sample = rownames(genus_test),
                    sum = rowSums(genus_test)) %>%
    mutate(unassigned = 100 - sum)
  # Use this to add Unassigned
  
  # Make mctoolsr object for easier analyses
  taxonomy_loaded <- data.frame(Taxonomy = colnames(genus_test)) %>%
    add_row(Taxonomy  = "Root; d__Unassigned; p__Unassigned; c__Unassigned; o__Unassigned; f__Unassigned; g__Unassigned") %>%
    mutate(Taxonomy = gsub("Root; |d__|p__|c__|o__|f__|g__", "", Taxonomy)) %>%
    separate(Taxonomy, sep = "; ", remove = F,
             into = c("taxonomy1", "taxonomy2", "taxonomy3", "taxonomy4", "taxonomy5", "taxonomy6")) %>%
    left_join(., ind_merge, by = "Taxonomy") %>%
    rename(taxonomy7 = Habitat) %>%
    mutate(taxonomy7 = ifelse(is.na(taxonomy7), "non-indicator", taxonomy7)) %>%
    mutate(taxonomy7 = ifelse(taxonomy6 == "Unassigned", "unassigned", taxonomy7)) %>%
    dplyr::select(-Taxonomy) %>%
    mutate(ID = taxonomy6) %>%
    column_to_rownames(var = "ID")
  map_loaded <- meta_test %>%
    mutate(sampleID = sample) %>%
    column_to_rownames(var = "sample")
  data_loaded <- genus_test %>%
    mutate("Root; d__Unassigned; p__Unassigned; c__Unassigned; o__Unassigned; f__Unassigned; g__Unassigned"= cs1$unassigned) %>%
    t() %>%
    as.data.frame() %>%
    set_names(meta_test$sample) %>%
    rownames_to_column(var = "Taxonomy") %>%
    mutate(Taxonomy = gsub("Root; |d__|p__|c__|o__|f__|g__", "", Taxonomy)) %>%
    separate(Taxonomy, sep = "; ", remove = T,
             into = c("taxonomy1", "taxonomy2", "taxonomy3", "taxonomy4", "taxonomy5", "taxonomy6")) %>%
    dplyr::select(-taxonomy1, -taxonomy2, -taxonomy3, -taxonomy4, -taxonomy5) %>%
    column_to_rownames(var = "taxonomy6")
  if (sum(rownames(data_loaded) != rownames(taxonomy_loaded)) > 0) { 
    stop("Taxonomy order in data_loaded and taxonomy_loaded doesn't match.") } 
  if (sum(colnames(data_loaded) != rownames(map_loaded)) > 0) { 
    stop("Sample order in data_loaded and map_loaded doesn't match.") } 
  input <- list()
  input$map_loaded <- map_loaded
  input$data_loaded <- data_loaded
  input$taxonomy_loaded <- taxonomy_loaded
  logtime("Finished making mctoolsr-style input")
  
  # Summarize taxonomy
  tax_sum_ind <- summarize_taxonomy(input = input,
                                    level = 7,
                                    relative = FALSE,
                                    report_higher_tax = FALSE)
  logtime("Finished summarizing taxonomy")
  
  # Summarize data for plotting
  plot_data <- plot_taxa_bars(tax_table = tax_sum_ind,
                              metadata_map = input$map_loaded,
                              type_header = "sampleID",
                              num_taxa = 100,
                              data_only = TRUE) %>%
    left_join(., input$map_loaded, by = c("group_by" = "sampleID")) %>%
    mutate(match = ifelse(taxon == Habitat, "target", "off-target")) %>%
    mutate(taxon = gsub("unassigned", "genus unassigned", taxon))
  
  if (!keep_unassigned) {
    plot_data <- plot_data %>%
      dplyr::filter(taxon != "genus unassigned")
    logtime("Genus unassigned removed from plot_data")
  }
  
  if (!keep_nonindicator) {
    plot_data <- plot_data %>%
      dplyr::filter(taxon != "non-indicator")
    logtime("Non-indicators removed from plot_data")
  }
  
  if (!is.null(custom_order)) {
    plot_data$taxon <- factor(plot_data$taxon, levels = custom_order)
    plot_data$Habitat <- factor(plot_data$Habitat, levels = custom_order)
    logtime("Re-ordered the habitats for plotting")
  }
  
  logtime("Finished preparing plotting data")
  
  # Print median and mean for target and non-target
  info <- plot_data %>%
    group_by(taxon, match) %>%
    summarise(Mean = mean(mean_value),
              Median = median(mean_value)) %>%
    ungroup()
  print(info, n = 100)
  
  # Plot
  g <- ggplot(plot_data, aes(Habitat, mean_value)) +
    geom_boxplot(outliers = F, colour = "blue") +
    geom_jitter(size = 1, alpha = 0.25, pch = 16, width = 0.25, height = 0,
                aes(colour = match)) +
    scale_colour_manual(values = c("black", "red")) +
    labs(x = NULL,
         y = "% abundance") +
    facet_wrap(~ taxon, ncol = ncol) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  pdf(output, width = 8, height = 8)
  print(g)
  dev.off()
  logtime("Finished diagnostic plot. checkIndicPiper() complete")
}



#### Support Functions ####
# Support functions
se <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  stats::sd(x) / sqrt(sum(!is.na(x)))
}

logtime <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))

make_subsamples <- function(df, n_runs = 100, n_per_habitat = 215, seed = 1){
  set.seed(seed)
  split_ids <- split(df$sample, df$Habitat)
  res <- vector("list", n_runs)
  for(i in seq_len(n_runs)) res[[i]] <- data.frame(sample=character(),Habitat=character(),stringsAsFactors=FALSE)
  for(h in names(split_ids)){
    ids <- sample(split_ids[[h]])
    N <- length(ids); needed <- n_runs*n_per_habitat
    if(N < needed) ids <- rep(ids, length.out=needed)
    for(i in seq_len(n_runs)){
      idx <- ((i-1)*n_per_habitat+1):(i*n_per_habitat)
      res[[i]] <- rbind(res[[i]], data.frame(sample=ids[idx],Habitat=h,stringsAsFactors=FALSE))
    }
  }
  res
}

# 2 functions from the mctoolsr package (Leff et al. 2017)
summarize_taxonomy = function(input, level, relative = TRUE,
                              report_higher_tax = TRUE) {
  if (report_higher_tax)
    taxa_strings = apply(input$taxonomy_loaded[1:level], 1,
                         paste0, collapse = '; ')
  else
    taxa_strings = input$taxonomy_loaded[, level]
  no_taxa = length(unique(taxa_strings))
  tax_sum = apply(input$data_loaded, 2, function(x)
    by(x, taxa_strings, sum))
  if (no_taxa == 1)
    tax_sum = data.frame(t(tax_sum),
                         row.names = unique(taxa_strings))
  else
    tax_sum = as.data.frame(tax_sum)
  if (relative) {
    output = convert_to_relative_abundances(tax_sum)
    # warn if NAs produced from trying to divide by 0
    na_samples = names(output)[colSums(is.na(output)) > 0]
    if (length(na_samples) > 0) {
      warning(
        paste(
          'The following samples produced NAs:',
          paste(na_samples, collapse = ', '),
          '\nThis might be because they had no observation data.'
        )
      )
    }
    output
  } else
    tax_sum
}

plot_taxa_bars <- function(tax_table, metadata_map, type_header, num_taxa,
                           data_only = FALSE) {
  
  tax_table$taxon <- row.names(tax_table)
  
  tax_table_melted <- reshape2::melt(
    tax_table,
    variable.name = "Sample_ID",
    id.vars = "taxon"
  )
  
  group_by_levels <- metadata_map[
    match(tax_table_melted$Sample_ID, row.names(metadata_map)),
    type_header
  ]
  
  tax_table_melted$group_by <- group_by_levels
  
  mean_tax_vals <- tax_table_melted %>%
    dplyr::group_by(group_by, taxon) %>%
    dplyr::summarise(
      mean_value = mean(value),
      .groups = "drop"
    )
  
  # get top taxa and convert others to "Other"
  mean_tax_vals_sorted <- mean_tax_vals[
    order(mean_tax_vals$mean_value, decreasing = TRUE),
  ]
  
  top_taxa <- unique(mean_tax_vals_sorted$taxon)[1:num_taxa]
  
  mean_tax_vals_sorted$taxon[
    !mean_tax_vals_sorted$taxon %in% top_taxa
  ] <- "Other"
  
  to_plot <- mean_tax_vals_sorted %>%
    dplyr::group_by(group_by, taxon) %>%
    dplyr::summarise(
      mean_value = sum(mean_value),
      .groups = "drop"
    )
  
  # make it so "Other" appears at bottom of key
  o <- c(
    unique(to_plot$taxon)[unique(to_plot$taxon) != "Other"],
    "Other"
  )
  
  to_plot$taxon <- factor(
    x = to_plot$taxon,
    levels = o
  )
  
  to_plot <- to_plot[
    order(
      to_plot$group_by,
      as.numeric(to_plot$taxon),
      decreasing = TRUE
    ),
  ]
  
  if (data_only) {
    
    return(to_plot)
    
  } else {
    
    return(
      ggplot2::ggplot(
        to_plot,
        ggplot2::aes(
          x = group_by,
          y = mean_value,
          fill = taxon
        )
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::ylab("") +
        ggplot2::xlab("") +
        ggplot2::theme(
          legend.title = ggplot2::element_blank()
        )
    )
  }
}

