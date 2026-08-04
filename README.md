# IndicPiper
IndicPiper: microbial indicator taxa analysis using Sandpiper and multipatt

This repo contains the IndicPiper database made by Cliff Bueno de Mesquita based on 13 habitats, using 100 runs, 250 random samples per habitat per run, and cutoffs of indicator in 100% of runs, mean p-value < 0.01, and mean IndVal > 0.5. For many cases, you can just use this database for your projects. Whether you have metagenomes or 16S sequencing, just use GTDB taxonomy and perform exact name matching at the genus level. IndicPiper v2 uses GTDB r232 taxonomy. For GTDB r226 taxonomy, use IndicPiper v1 (see releases). Then, for example, you can aggregate relative abundances by indicator habitat, and perform statistics and plotting. 

The pre-built database is distributed as part of this repository rather than via Zenodo. It is committed as `genus_habitat_indicators_v2.csv` (GTDB r232 taxonomy) and `genus_habitat_indicators_v1.csv` (GTDB r226 taxonomy, see the `v1.0.0` release/tag). Cloning the repository or downloading a release ZIP as described in [Installation](#installation) gives you these files directly at the top level of the repo.

For example, `Example_PlotSandpiperSample.R` takes a GTDB genus-level taxonomic profile for sample SRR34514425 — a human gut metagenome from (Fernandes et al. 2026) — exact-matches its genera against `genus_habitat_indicators_v2.csv`, aggregates relative abundance by indicator habitat, and plots the result. The profile is committed to this repo as `SRR34514425_condensed.tsv` so the example doesn't depend on the Sandpiper API; it was obtained with:

```bash
wget -O SRR34514425_condensed.tsv \
  "https://sandpiper.qut.edu.au/api/condensed_csv_with_extras/SRR34514425?taxonomy_type=gtdb"
```

The key steps are:

```r
# Load the pre-built IndicPiper database (habitat indicator genera)
ind <- read.csv("genus_habitat_indicators_v2.csv")

# Load the GTDB condensed taxonomic profile for this sample
profile <- read.delim("SRR34514425_condensed.tsv")

# Keep genus-level rows and strip the rank prefixes so the taxonomy strings
# match the format used in genus_habitat_indicators_v2.csv
genus <- profile %>%
  filter(level == "genus") %>%
  mutate(Taxonomy = gsub("^Root; |d__|p__|c__|o__|f__|g__", "", taxonomy))

# Exact-match genera to their indicator habitat; anything not in the
# database is "Non-indicator"
genus_annotated <- genus %>%
  left_join(ind %>% select(Taxonomy, Habitat), by = "Taxonomy") %>%
  mutate(Habitat = ifelse(is.na(Habitat), "Non-indicator", Habitat))

# Aggregate relative abundance by indicator habitat
habitat_abund <- genus_annotated %>%
  group_by(Habitat) %>%
  summarise(relative_abundance = sum(relative_abundance), .groups = "drop")
```

Run the full script (including the plotting code) with:

```bash
pixi run Rscript Example_PlotSandpiperSample.R
```

which produces:

![Relative abundance of SRR34514425 by IndicPiper indicator habitat](img/SRR34514425_habitat_abundance.png)

There are also functions (in IndicPiper.R) to generate your own database based on habitats of interest or different parameters. The starting metadata and taxaonomic profile files are available on Zenodo (https://zenodo.org/records/20855888), and these were generated with GenerateStartingPoint.R. You can then supply the functions with your habitats of interest and cutoffs you want to use. We recommend not going any less stringent than the cutoffs we used, but you could potentially try more stringent cutoffs to get only the strongest associations. We also recommend focusing on habitats that have good sample sizes (ideally in the hundreds of samples). We have removed habitats with < 50 samples.  

There is also a function to generate a diagnostic plot (`checkIndicPiper()`) so you can see to what relative abundance the indicator taxa sum to in the target habitat as well as how much they spill over into other habitats.

## Installation

IndicPiper is distributed only through GitHub (https://github.com/cliffbueno/IndicPiper). There is no CRAN/Bioconda package to install. Get the code by either cloning the repository:

```bash
git clone https://github.com/cliffbueno/IndicPiper.git
cd IndicPiper
```

or downloading a release (or the current source) as a ZIP from the [GitHub page](https://github.com/cliffbueno/IndicPiper) and unzipping it. The functions live in `IndicPiper.R`, which you load with `source("IndicPiper.R")`.

To run IndicPiper, you only need R and a few R packages. IndicPiper was developed with R 4.5.2.

Install the required packages with:

```r
install.packages(c(
  "tidyr", "dplyr", "tibble", "permute", "indicspecies",
  "ggplot2", "data.table", "rlang",
  "R.utils", "reshape2"
))
```

Alternatively, if you use [pixi](https://pixi.sh), a `pixi.toml` is included that pins R and all of the required packages. Running any pixi command in the repository (for example `pixi shell`, or `pixi run test`) will install R and the packages into an isolated environment for you — no `install.packages()` needed:

```bash
pixi shell     # drop into a shell with R and all packages available, then run R / Rscript
```

You will also need the two starting input files, which can be downloaded from Zenodo:

https://zenodo.org/records/20855888

or from the command line using:

```bash
wget "https://zenodo.org/records/20855888/files/Sandpiper_Genus_Filt_n451568.csv.gz?download=1" \
  -O Sandpiper_Genus_Filt_n451568.csv.gz

wget "https://zenodo.org/records/20855888/files/Sandpiper_Metadata_Filt_n451568.txt?download=1" \
  -O Sandpiper_Metadata_Filt_n451568.txt
```

## Usage
To use the provided database, generate GTDB taxonomic abundance profiles from metagenomes or 16S rRNA gene sequencing and then exact match by genus name to assign genera as "non-indicator" or as indicators of the habitats according to the IndicPiper output. Then you can just aggregate by indicator taxa and plot relative abundances as you would for any other aggregated taxonomic level like phylum.

## Generating a custom database

You can also generate your own database (for example, if you need a habitat not in the 13 provided).

To do a custom run of IndicPiper, download the two input files from Zenodo https://zenodo.org/records/20855888.
Then load all of the functions with source(IndicPiper.R). You could then use the functions interactively on RStudio Server, or make an .R script in which you set the working directory, run source(IndicPiper.R), and then run one of the functions, providing your arguments. Such a file could then be run in the terminal with Rscript YourFile.R. For examples of each function, I have provided 4 "Test" scripts, which demonstrate this. I've also provided some example slurm scripts for how I ran this on the University of Colorado supercomputer; those would have to be updated according to your computer. Note that due to the large file sizes, 275 Gb RAM is needed for prepIndicPiper(). Less RAM is needed for runIndicPiper(). For example, for the 13 habitats I made the database for, 192 Gb RAM were needed.
IndicPiper has 4 main functions:\
  `countHabitats()`\
  `prepIndicPiper()`\
  `runIndicPiper()`\
  `checkIndicPiper()`

1. Use `countHabitats()` to import the provided metadata table and count number of samples by habitat. This will help you decide which habitats you can use and which ones you want to test.

   Arguments:

   - `meta`: the metadata file path.  
     Default = `"Sandpiper_Metadata_Filt_n451568.txt"`

2. Use the `prepIndicPiper()` function to filter the metadata table and genus relative abundance table to your habitats of interest and optionally combine certain habitats. The combine arguments won't affect you if you aren't using those habitats. Don't combine if you need a specific habitat that is otherwise combined by default.

   Arguments:

   - `habitat_list`: a vector of habitat names.

   - `combine_soil_rhizo`: TRUE/FALSE.  
     Default = `TRUE`

   - `combine_freshwater`: TRUE/FALSE.  
     Default = `TRUE`

   - `combine_glacier_ice`: TRUE/FALSE.  
     Default = `TRUE`

   - `combine_mammalian_gut`: TRUE/FALSE.  
     Default = `TRUE`

   - `combine_saliva_oral`: TRUE/FALSE.  
     Default = `TRUE`

3. Use the `runIndicPiper()` function to perform the random subsetting, multipatt analysis, and merging/filtering of multipatt output tables. This builds an output table with genera and their indicator habitats, as well as information about the strength of the association.

   Arguments:

   - `meta`: the metadata table output from `prepIndicPiper()`.  
     Default = `"myMetadataTable.csv.gz"`

   - `genus`: the genus relative abundance table output from `prepIndicPiper()`.  
     Default = `"myGenusTable.csv.gz"`

   - `n_multipatt_perm`: the number of iterations within the multipatt function.  
     Default = `100`

   - `n_runs`: the number of times you want to subsample habitats and run multipatt.  
     Default = `100`

   - `n_per_habitat`: the sample size per habitat per run. Samples are randomly selected. Must be less than the lowest number of samples in any one habitat from your `habitat_list`. Check the habitat counts at the end of `prepIndicPiper()`.  
     Default = `250`

   - `run_cut`: cutoff for filtering multipatt output. The percent of runs that a given genus is an indicator of the same habitat.  
     Default = `100`

   - `p_cut`: cutoff for filtering multipatt output. The multipatt mean p-value.  
     Default = `0.01`

   - `IndVal_cut`: cutoff for filtering multipatt output. The multipatt mean IndVal.  
     Default = `0.5`

   - `seed`: seed for reproducibility.  
     Default = `1`

4. Use the `checkIndicPiper()` function to plot the summed abundances of the indicator taxa.

   Arguments:

   - `meta`: a metadata table from one of your runs (e.g., made with `runIndicPiper()`).  
     Default = `"meta_test.csv"`

   - `genus`: a genus relative abundance table from one of your runs (e.g., made with `runIndicPiper()`).  
     Default = `"genus_test.csv"`

   - `ind`: an IndicPiper database (e.g., made with `runIndicPiper()`).  
     Default = `"genus_habitat_indicators_custom.csv"`

   - `output`: name of the output plot.  
     Default = `"IndicPiper_SummedAbund.pdf"`

   - `ncol`: number of columns in the `facet_wrap()`.  
     Default = `2`

   - `keep_unassigned`: whether to keep the **Unassigned** facet.  
     Default = `FALSE`

   - `keep_nonindicator`: whether to keep the **Non-indicator** facet.  
     Default = `TRUE`

   - `custom_order`: a custom order for the facets. Supply a character vector of habitat names in the desired order.  
     Default = `NULL`

## Testing
A test suite exercises all four IndicPiper functions end-to-end. Because the real Zenodo inputs are very large (hundreds of thousands of samples, needing ~275 Gb RAM), the tests instead run against a small synthetic dataset generated on the fly, so the whole suite finishes in seconds.

The R dependencies are managed with [pixi](https://pixi.sh). To run the tests locally:

```bash
pixi run test
```

This generates the mock data (`tests/generate_mock_data.R`) and then runs the test suite (`tests/run_tests.R`). The tests also run automatically on every push and pull request via GitHub Actions (see `.github/workflows/test.yml`).

## Resources
We recommend running IndicPiper on a server or supercomputer due to the size of the databases and the heavy computation needed to run all of the iterations of multipatt on the large input tables. IndicPiper v2 was developed on a supercomputer with 275 Gb RAM. `countHabitats` took 1 minute. `prepIndicPiper` took 44 minutes with 20 cores. `runIndicPiper` took 7 hours with 16 cores for 100 runs. `checkIndicPiper` took 22 seconds.

## References
IndicPiper will be described in a forthcoming publication. If you need to cite IndicPiper before the paper is out, please cite this GitHub repository. Please also mention that IndicPiper relies on the Sandpiper database and cite the Sandpiper/SingleM paper too (Woodcroft et al. 2025). We also suggest citing the paper associated with the `multipatt()` R function used by IndicPiper (De Cáceres and Legendre, 2009). Lastly, if you use `checkIndicPiper()`, it uses functions from Jonathan Leff's mctoolsr package (Leff, 2022).

De Cáceres, M. and Legendre, P. (2009), Associations between species and groups of sites: indices and statistical inference. *Ecology*, 90: 3566–3574. https://doi.org/10.1890/08-1823.1

Leff, J. 2022. mctoolsr: Microbial Community Data Analysis Tools. R package version 0.1.1.9. <https://github.com/leffj/mctoolsr>

Woodcroft, B.J., Aroney, S.T.N., Zhao, R. *et al.* Comprehensive taxonomic identification of microbial species in metagenomic data using SingleM and Sandpiper. *Nat Biotechnol* (2025). https://doi.org/10.1038/s41587-025-02738-1

Fernandes, R., Jabbarizadeh, B., Rajeh, A. *et al.* Fecal microbiota transplantation plus immunotherapy in metastatic renal cell carcinoma: the phase 1 PERFORM trial. *Nat Med* 32, 1325–1336 (2026). https://doi.org/10.1038/s41591-025-04183-8
