# IndicPiper
IndicPiper: microbial indicator taxa analysis using Sandpiper and multipatt

This repo contains the IndicPiper database made by Cliff Bueno de Mesquita based on 13 habitats, using 100 runs, 215 random samples per habitat per run, and cutoffs of indicator in 100% of runs, mean p-value < 0.01, and mean IndVal > 0.5. For many cases, you can just use this database for your projects. Whether you have metagenomes or 16S sequencing, just use GTDB taxonomy and perform exact name matching at the genus level. Then, for example, you can aggregate relative abundances by indicator habitat, and perform statistics and plotting. 

There are also functions (in IndicPiper.R) to generate your own database based on habitats of interest or different parameters. The input metadata and taxaonomic profile are available on FigShare, and these were generated with GenerateStartingPoint.R. You can then supply the function with your habitats of interest and cutoffs you want to use. We recommend not going any less stringent than the cutoffs we used. We also recommend focusing on habitats that have good sample sizes (ideally in the hundreds of samples). We have removed habitats with < 50 samples.  

There is also a function to generate a diagnostic plot (`checkIndicPiper()`) so you can see to what relative abundance the indicator taxa sum to in the target habitat as well as how much they spill over into other habitats.

## Installation
To run IndicPiper, you just need R and few libraries. IndicPiper was developed with R 4.5.2. Required libraries are tidyr, dplyr, tibble, permute, indicspecies, ggplot2, data.table, FSA, rlang and R.utils. These can be installed with:\
`install.packages(c("tidyr", "dplyr", "tibble", "permute", "indicspecies", "ggplot2", "data.table", "FSA", "rlang", "R.utils"))`.\
You also need to the two starting input files, which can be downloaded here: https://doi.org/10.6084/m9.figshare.32305302.

## Usage
To use the provided database, generate GTDB taxonomic abundance profiles from metagenomes or 16S rRNA gene sequencing and then exact match by genus name to assign genera as "non-indicator" or as indicators of the habitats according to the IndicPiper output. Then you can just aggregate by indicator taxa and plot relative abundances as you would other taxa or guilds.
 
You can also generate your own database (for example, if you need a habitat not in the 13 provided).

To do a custom run of IndicPiper, download the two input files from FigShare (https://doi.org/10.6084/m9.figshare.32305302).
Then load all of the functions with source(IndicPiper.R). You could then use the functions interactively on RStudio Server, or make an .R script in which you set the working directory, run source(IndicPiper.R), and then run one of the functions, providing your arguments. Such a file could then be run in the terminal with Rscript YourFile.R. For examples of each function, I have provided 4 "Test" scripts, which demonstrate this.
IndicPiper has 4 main functions:\
  `countHabitats()`\
  `prepIndicPiper()`\
  `runIndicPiper()`\
  `checkIndicPiper()`

1. Use `countHabitats()` to import the provided metadata table and count number of samples by habitat. This will help you decide which habitats you can use and which ones you want to test.

   Arguments:

   - `meta`: the metadata file path.  
     Default = `"Sandpiper_Metadata_Filt_n358209.txt"`

2. Use the `prepIndicPiper()` function to filter the metadata table and genus relative abundance table to your habitats of interest and optionally combine certain habitats.

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
     Default = `215`

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

## Resources
We recommend running IndicPiper on a server or supercomputer due to the size of the databases and the heavy computation needed to run all of the iterations of multipatt on the large input tables. IndicPiper was developed on a server with 250 Gb RAM and 32 cores. `countHabitats` took 1 minute. `prepIndicPiper` took 24 minutes. `runIndicPiper` took 7 hours for 100 runs or 15 minutes for 5 runs (testing). `checkIndicPiper` took 30 seconds.

## References
IndicPiper will be described in a forthcoming publication. If you need to cite IndicPiper before the paper is out, please cite this GitHub repository. Please also mention that IndicPiper relies on the Sandpiper database and cite the Sandpiper/SingleM paper too (Woodcroft et al. 2025). We also suggest citing the paper associated with the `multipatt()` R function used by IndicPiper (De Cáceres and Legendre, 2009). Lastly, if you use `checkIndicPiper()`, it uses functions from Jonathan Leff's mctoolsr package (Leff, 2022).

De Cáceres, M. and Legendre, P. (2009), Associations between species and groups of sites: indices and statistical inference. *Ecology*, 90: 3566–3574. https://doi.org/10.1890/08-1823.1

Leff, J. 2022. mctoolsr: Microbial Community Data Analysis Tools. R package version 0.1.1.9. <https://github.com/leffj/mctoolsr>

Woodcroft, B.J., Aroney, S.T.N., Zhao, R. *et al.* Comprehensive taxonomic identification of microbial species in metagenomic data using SingleM and Sandpiper. *Nat Biotechnol* (2025). https://doi.org/10.1038/s41587-025-02738-1
