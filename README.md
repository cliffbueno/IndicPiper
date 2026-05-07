# IndicPiper
IndicPiper: microbial indicator taxa analysis using Sandpiper and multipatt

This repo contains the IndicPiper database made by Cliff Bueno de Mesquita based on 13 habitats, using 100 runs, 215 random samples per habitat per run, and cutoffs of indicator in 100% of runs, mean p-value < 0.01, and mean IndVal > 0.5. For many cases, you can just use this database for your projects. Whether you have metagenomes or 16S sequencing, just use GTDB taxonomy and perform exact name matching at the genus level. Then, for example, you can aggregate relative abundances by indicator habitat. 

There is also a function to generate your own database based on habitats of interest or different parameters. The input metadata and taxaonomic profile are available on FigShare. You can then supply the function with your habitats of interest and cutoffs you want to use. We recommend not going any less stringent than the cutoffs we used. We also recommend focusing on habitats that have good sample sizes (ideally in the hundreds of samples). 

There is also a function to generate a diagnostic plot so you can see to what relative abundance the indicator taxa sum to in the target habitat as well as how much they spillover into other habitats.

## Installation
To run IndicPiper, you just need R and few libraries. IndicPiper was developed with R 4.5.2.
Required libraries are dplyr, permute, indicspecies, ggplot2. These can be installed with:
install.packages(c("dplyr"), "permute", "indicspecies", "ggplot2")

## Usage
To use the provided database, generate GTDB taxonomic abundance profiles from metagenomes or 16S sequencing and then match by genus name. Then you can just aggregate by indicator taxa and plot relative abundances as you would other taxa or guilds.
Or, you can generate your own database (for example, if you need a habitat not in the 13 provided).

To run IndicPiper, download the two input files from FigShare.
Then load all of the functions with source(IndicPiper.R)
Then use the prepIndicPiper() function. The arguments are:

Then use the runIndicPiper() function. The arguments are:

Then use the checkIndicPiper() function to plot the summed abundances of the indicator taxa. The arguments are:

Once you have your database
