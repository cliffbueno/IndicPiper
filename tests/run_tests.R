# IndicPiper test suite.
#
# Runs all four IndicPiper functions end-to-end against the small mock dataset
# produced by tests/generate_mock_data.R, checking that each stage completes and
# produces output with the expected structure. Exits non-zero on any failure so
# it can be used directly in CI.
#
# Usage:  Rscript tests/run_tests.R      (run from the repo root)

# ---- Locate repo root and mock data --------------------------------------
# Resolve the directory of this script so the test can be run from anywhere.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/run_tests.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."))

indicpiper_R <- file.path(repo_root, "IndicPiper.R")
test_dir <- Sys.getenv("INDICPIPER_TEST_DIR",
                       unset = file.path(repo_root, "tests", "tmp"))
test_dir <- normalizePath(test_dir, mustWork = FALSE)

meta_file  <- file.path(test_dir, "Sandpiper_Metadata_Filt_n451568.txt")
genus_file <- file.path(test_dir, "Sandpiper_Genus_Filt_n451568.csv.gz")
if (!file.exists(meta_file) || !file.exists(genus_file)) {
  stop("Mock data not found. Run tests/generate_mock_data.R first ",
       "(or `pixi run mock-data`).")
}

# ---- Tiny test harness ---------------------------------------------------
n_pass <- 0L
n_fail <- 0L
check <- function(desc, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    message("    error: ", conditionMessage(e)); FALSE
  })
  if (ok) {
    n_pass <<- n_pass + 1L
    cat(sprintf("  PASS: %s\n", desc))
  } else {
    n_fail <<- n_fail + 1L
    cat(sprintf("  FAIL: %s\n", desc))
  }
}

# The IndicPiper functions read/write files in the working directory, so run
# everything inside the mock-data directory.
setwd(test_dir)
source(indicpiper_R)

# ---- 1. countHabitats() --------------------------------------------------
cat("== countHabitats() ==\n")
check("countHabitats runs without error", {
  countHabitats(meta = basename(meta_file)); TRUE
})

# ---- 2. prepIndicPiper() -------------------------------------------------
cat("== prepIndicPiper() ==\n")
habitat_list <- c("soil", "marine sediment", "human oral", "glacier")
check("prepIndicPiper runs and writes prepared tables", {
  prepIndicPiper(habitat_list = habitat_list,
                 combine_soil_rhizo = TRUE,
                 combine_freshwater = TRUE,
                 combine_glacier_ice = TRUE,
                 combine_mammalian_gut = TRUE,
                 combine_saliva_oral = TRUE)
  file.exists("myMetadataTable.csv.gz") && file.exists("myGenusTable.csv.gz")
})
check("combine logic folded rhizosphere/saliva into final habitats", {
  m <- data.table::fread("myMetadataTable.csv.gz")
  # source labels rhizosphere and human saliva must be gone; final set is the 4
  # target habitats (with seawater renamed to marine water where relevant).
  all(sort(unique(m$Habitat)) == sort(c("soil", "marine sediment",
                                        "human oral", "glacier")))
})

# ---- 3. runIndicPiper() --------------------------------------------------
cat("== runIndicPiper() ==\n")
# Small parameters keep the run fast; the mock signal is strong enough that
# indicators are still recovered in every run. n_multipatt_perm must be large
# enough that the smallest achievable p-value (1/(nperm+1)) is below p_cut.
check("runIndicPiper runs and writes an indicator table", {
  runIndicPiper(meta = "myMetadataTable.csv.gz",
                genus = "myGenusTable.csv.gz",
                n_multipatt_perm = 199,
                n_runs = 3,
                n_per_habitat = 20,
                run_cut = 100,
                p_cut = 0.05,
                IndVal_cut = 0.5,
                seed = 1,
                output = "genus_habitat_indicators_test.csv")
  file.exists("genus_habitat_indicators_test.csv")
})
check("indicator table has expected columns and at least one indicator", {
  ind <- read.csv("genus_habitat_indicators_test.csv")
  expected_cols <- c("Taxonomy", "Habitat", "PercRunsInd", "IndVal_mean",
                     "PVal_mean")
  nrow(ind) > 0 && all(expected_cols %in% colnames(ind))
})
check("recovered indicators map to the correct habitats", {
  ind <- read.csv("genus_habitat_indicators_test.csv")
  # Each mock indicator genus is named "<HabitatSlug>_ind<N>"; the last
  # taxonomy field should encode the habitat it was planted in.
  slug <- function(h) gsub("[^A-Za-z]", "", h)
  genus <- sub(".*; ", "", ind$Taxonomy)
  expected_slug <- vapply(ind$Habitat, slug, character(1))
  # every recovered genus name starts with the slug of its assigned habitat
  all(mapply(function(g, s) startsWith(g, s), genus, expected_slug))
})
check("runIndicPiper wrote per-iteration files for checkIndicPiper", {
  file.exists("meta_test.csv") && file.exists("genus_test.csv")
})

# ---- 4. checkIndicPiper() ------------------------------------------------
cat("== checkIndicPiper() ==\n")
check("checkIndicPiper runs and writes the diagnostic PDF", {
  checkIndicPiper(meta = "meta_test.csv",
                  genus = "genus_test.csv",
                  ind = "genus_habitat_indicators_test.csv",
                  output = "IndicPiper_diagnostic_test.pdf",
                  ncol = 2,
                  keep_unassigned = FALSE,
                  keep_nonindicator = TRUE)
  file.exists("IndicPiper_diagnostic_test.pdf")
})

# ---- Summary -------------------------------------------------------------
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1, save = "no")
cat("All IndicPiper tests passed.\n")
