# IndicPiper
IndicPiper: microbial indicator taxa analysis using Sandpiper and multipatt

This repo contains the IndicPiper database made by Cliff Bueno de Mesquita based on 13 habitats, using 100 runs, 215 random samples per habitat per run, and cutoffs of indicator in 100% of runs, mean p-value < 0.01, and mean IndVal > 0.5. For many cases, you can just use this database for your projects. Whether you have metagenomes or 16S sequencing, just use GTDB taxonomy and perform exact name matching at the genus level. Then, for example, you can aggregate relative abundances by indicator habitat. 

There is also a function to generate your own database based on habitats of interest or different parameters. The input metadata and taxaonomic profile are available on FigShare. You can then supply the function with your habitats of interest and cutoffs you want to use. We recommend not going any less stringent than the cutoffs we used. We also recommend focusing on habitats that have good sample sizes (ideally in the hundreds of samples). 
