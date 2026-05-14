# IndicPiper
IndicPiper: microbial indicator taxa analysis using Sandpiper and multipatt

This repo contains the IndicPiper database made by Cliff Bueno de Mesquita based on 13 habitats, using 100 runs, 215 random samples per habitat per run, and cutoffs of indicator in 100% of runs, mean p-value < 0.01, and mean IndVal > 0.5. For many cases, you can just use this database for your projects. Whether you have metagenomes or 16S sequencing, just use GTDB taxonomy and perform exact name matching at the genus level. Then, for example, you can aggregate relative abundances by indicator habitat. 

There is also a function to generate your own database based on habitats of interest or different parameters. The input metadata and taxaonomic profile are available on FigShare. You can then supply the function with your habitats of interest and cutoffs you want to use. We recommend not going any less stringent than the cutoffs we used. We also recommend focusing on habitats that have good sample sizes (ideally in the hundreds of samples). 

There is also a function to generate a diagnostic plot so you can see to what relative abundance the indicator taxa sum to in the target habitat as well as how much they spillover into other habitats.

## Installation
To run IndicPiper, you just need R and few libraries. IndicPiper was developed with R 4.5.2.
Required libraries are dplyr, tibble, permute, indicspecies, ggplot2, and data.table. These can be installed with:
`install.packages(c("dplyr", "tibble", "permute", "indicspecies", "ggplot2", "data.table")`

## Usage
To use the provided database, generate GTDB taxonomic abundance profiles from metagenomes or 16S rRNA gene sequencing and then exact match by genus name to assign genera as "non-indicator" or as indicators of the habitats according to the IndicPiper output. Then you can just aggregate by indicator taxa and plot relative abundances as you would other taxa or guilds.
 
You can also generate your own database (for example, if you need a habitat not in the 13 provided).

To do a custom run of IndicPiper, download the two input files from FigShare.
Then load all of the functions with source(IndicPiper.R).
IndicPiper has 4 main functions:
  `countHabitats()`
  `prepIndicPiper()`
  `runIndicPiper()`
  `checkIndicPiper()`

1. Use `countHabitats()` to import the provided metadata table and count number of samples by habitat. This will help you decide which habitats you can use and which ones you want to test. The arguments are: 
  meta: the metadata file path, default = Sandpiper_Metadata_Filt_n358209.txt

2. Use the `prepIndicPiper()` function to filter the metadata table and genus relative abundance table to your habitats of interest and also combine certain habitats if you want. The arguments are:
   habitat_list: a vector of habitat names.
   combine_soil_rhizo: TRUE/FALSE, default = TRUE
   combine_freshwater: TRUE/FALSE, default = TRUE
   combine_glacier_ice: TRUE/FALSE, default = TRUE
   combine_mammalian_gut: TRUE/FALSE, default = TRUE
   combine_saliva_oral: TRUE/FALSE, default = TRUE
   
3. Use the `runIndicPiper()` function to perform the random subsetting, multipatt analysis, and merging and filtering of multipatt output tables. This builds and output table with genera and their indicator habitats, as well as information about the strength of the association The arguments are:
   meta: the metadata table, output from `prepIndicPiper()`, default = "myMetadataTable.csv.gz",
   genus: the genus relative abundance table, output from `prepIndicPiper()`, default = "myGenusTable.csv.gz",
   n_multipatt_perm: the number of iterations within the multipatt function, default = 100
   n_runs: the number of times you want to subsample habitats and run multipatt, default = 100
   n_per_habitat: the sample size per habitat per run. samples are randomly selected. must be less than the lowest number of samples in any one of your habitat_list. Check the habitat counts at the end of `prepIndicPiper()`. default = 215
   run_cut: cutoff for filtering multipatt output. the percent of runs that a given genera is an indicator or the same habitat. default = 100
   p_cut: cutoff for filtering multipatt output. the multipatt mean p-value. default = 0.01
   IndVal_cut: cutoff for filtering multipatt output. the multipatt mean IndVal. default = 0.5
   seed: seed for reproducibility. default = 1.

5. Use the `checkIndicPiper()` function to plot the summed abundances of the indicator taxa. The arguments are:

Once you have your database

## Resources
We recommend running IndicPiper on a server or supercomputer due to the size of the databases and the heavy computation needed to run all of the iterations of multipatt on the large input tables. IndicPiper was developed on a server with 250 Gb RAM and 32 cores. `countHabitats` took 1 minute. `prepIndicPiper` took 24 minutes. `runIndicPiper` took 7 hours. `checkIndicPiper` took 3 minutes.

## Reference
IndicPiper will be described in a forthcoming publication.
If you need to cite IndicPiper before the paper is out, please cite this GitHub repository.
Please also mention that IndicPiper relies the Sandpiper database and cite the Sandpiper paper too. 
