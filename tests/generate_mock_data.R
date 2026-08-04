# Generate small mock input data for the IndicPiper test suite.
#
# The real IndicPiper inputs live on Zenodo and are enormous (hundreds of
# thousands of samples, needing ~275 Gb RAM to process). For automated testing
# we instead synthesise a tiny dataset with the same *structure* so the four
# IndicPiper functions can be exercised end-to-end in seconds.
#
# Two files are produced, with the exact names prepIndicPiper() expects:
#   Sandpiper_Metadata_Filt_n451568.txt   - tab-delimited, cols: sample, Habitat
#   Sandpiper_Genus_Filt_n451568.csv.gz   - CSV, first col "taxonomy" (6-level
#                                            GTDB string), then one column per
#                                            sample of % relative abundances.
#
# The mock data is designed so a handful of "indicator" genera are strongly and
# reproducibly associated with each habitat, guaranteeing runIndicPiper()
# produces a non-empty indicator table.

suppressMessages({
  library(data.table)
})

set.seed(42)

# Where to write the mock inputs. Defaults to tests/tmp relative to repo root.
out_dir <- Sys.getenv("INDICPIPER_TEST_DIR", unset = "tests/tmp")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta_file  <- file.path(out_dir, "Sandpiper_Metadata_Filt_n451568.txt")
genus_file <- file.path(out_dir, "Sandpiper_Genus_Filt_n451568.csv.gz")

# ---- Metadata ------------------------------------------------------------
# Four final habitats. We deliberately include the "rhizosphere" and
# "human saliva" source labels so prepIndicPiper()'s combine logic is exercised
# (rhizosphere -> soil, human saliva -> human oral).
habitat_spec <- list(
  soil            = c(soil = 25, rhizosphere = 15),
  `marine sediment` = c(`marine sediment` = 40),
  `human oral`    = c(`human oral` = 25, `human saliva` = 15),
  glacier         = c(glacier = 40)
)

meta <- rbindlist(lapply(names(habitat_spec), function(final_hab) {
  labels <- habitat_spec[[final_hab]]
  rbindlist(lapply(names(labels), function(src) {
    data.table(
      final_habitat = final_hab,
      Habitat = src,
      n = labels[[src]]
    )
  }))
}))
meta <- meta[rep(seq_len(nrow(meta)), meta$n)]
meta[, n := NULL]
meta[, sample := sprintf("SRR%06d", seq_len(.N))]

# Final habitat used only to drive the abundance signal below; the written
# metadata keeps the *source* Habitat label (combine logic reproduces the rest).
final_habitat <- meta$final_habitat
samples <- meta$sample
n_samples <- length(samples)
final_habitats <- unique(final_habitat)

# ---- Genus abundance table ----------------------------------------------
# 5 indicator genera per habitat + 10 background genera present everywhere.
make_tax <- function(domain, phylum, class, order, family, genus) {
  sprintf("Root; d__%s; p__%s; c__%s; o__%s; f__%s; g__%s",
          domain, phylum, class, order, family, genus)
}

indicator_taxa <- list()
for (h in final_habitats) {
  hslug <- gsub("[^A-Za-z]", "", h)
  for (i in seq_len(5)) {
    g <- sprintf("%s_ind%d", hslug, i)
    indicator_taxa[[make_tax("Bacteria", paste0("P_", hslug),
                             paste0("C_", hslug), paste0("O_", hslug),
                             paste0("F_", hslug), g)]] <- h
  }
}
background_taxa <- vapply(seq_len(10), function(i) {
  make_tax("Bacteria", "P_Background", "C_Background", "O_Background",
           "F_Background", sprintf("Background_%02d", i))
}, character(1))

all_taxa <- c(names(indicator_taxa), background_taxa)
n_taxa <- length(all_taxa)

# Build the abundance matrix (taxa x samples), values are percentages.
mat <- matrix(0, nrow = n_taxa, ncol = n_samples,
              dimnames = list(all_taxa, samples))

# Indicator genera: high in their habitat, near-zero elsewhere.
for (tax in names(indicator_taxa)) {
  target <- indicator_taxa[[tax]]
  is_target <- final_habitat == target
  mat[tax, is_target]  <- runif(sum(is_target), min = 10, max = 25)
  mat[tax, !is_target] <- runif(sum(!is_target), min = 0, max = 0.4)
}
# Background genera: low, present everywhere.
for (tax in background_taxa) {
  mat[tax, ] <- runif(n_samples, min = 0.5, max = 3)
}

genus_dt <- data.table(taxonomy = all_taxa)
genus_dt <- cbind(genus_dt, as.data.table(mat))

# ---- Write ---------------------------------------------------------------
write.table(meta[, .(sample, Habitat)], meta_file,
            sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(genus_dt, genus_file, compress = "gzip", quote = FALSE)

cat(sprintf("Wrote %d samples across %d habitats and %d genera to:\n  %s\n  %s\n",
            n_samples, length(final_habitats), n_taxa, meta_file, genus_file))
