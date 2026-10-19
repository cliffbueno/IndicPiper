# 3 applications of IndicPiper
# For ISME Comms manuscript
# By Cliff Bueno de Mesquita, April-June 2026, revised October 2026
# Note v1 uses GTDB r226, v2 uses GTDB r232
# v2.1.0 was updated during ISME Comms revisions



#### 1. Setup ####
# Libraries
library(plyr)
library(tibble)
library(tidyr)
library(dplyr)
library(rlang)
library(indicspecies)
library(data.table)
library(readr)
library(Matrix)
library(ggplot2)
library(vegan)
library(scales)
library(plotly)
library(mctoolsr)
library(car)

# Working directory
setwd("~/Desktop/Fierer/AEGIS/IndicPiper/")

# Functions
`%notin%` <- Negate(`%in%`)
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]

# Indicator database
InDB <- read.csv("genus_habitat_indicators_v1.csv")
table(InDB$Habitat)
InDB_formerge <- InDB %>%
  dplyr::select(Genus, Habitat)

InDBv2 <- read.csv("genus_habitat_indicators_v2.1.0.csv") %>%
  separate(Taxonomy, remove = F, sep = "; ",
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
table(InDBv2$Habitat)
InDBv2_formerge <- InDBv2 %>%
  dplyr::select(Genus, Habitat)
table(InDBv2_formerge$Habitat)

InDB_pal <- data.frame(Habitat = c("freshwater sediment",
                                   "freshwater water",
                                   "glacier",
                                   "human oral",
                                   "human skin",
                                   "human vaginal",
                                   "hydrothermal vent",
                                   "mammalian gut",
                                   "marine sediment",
                                   "marine water",
                                   "phyllosphere",
                                   "soil",
                                   "wastewater"),
                       Color = c(brewer_pal(palette = "Paired")(12)[1], 
                                 brewer_pal(palette = "Paired")(12)[2], 
                                 "cyan",
                                 brewer_pal(palette = "Paired")(12)[8],
                                 brewer_pal(palette = "Paired")(12)[5],
                                 brewer_pal(palette = "Paired")(12)[6],
                                 brewer_pal(palette = "Paired")(12)[11],
                                 brewer_pal(palette = "Paired")(12)[7],
                                 brewer_pal(palette = "Paired")(12)[9],
                                 brewer_pal(palette = "Paired")(12)[10],
                                 brewer_pal(palette = "Paired")(12)[3],
                                 brewer_pal(palette = "Paired")(12)[4],
                                 brewer_pal(palette = "Paired")(12)[12]))
InDB_cols <- setNames(InDB_pal$Color, InDB_pal$Habitat)



#### 2. Figure 3 ####
# Make new figure 3 showing counts by phylum and habitat
ip2 <- read.csv("genus_habitat_indicators_v2.1.0.csv") %>%
  separate(Taxonomy, remove = F, sep = "; ",
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
phylum_counts <- ip2 %>%
  group_by(Habitat, Phylum) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Habitat, Phylum, fill = list(Count = 0)) %>%
  mutate(Habitat = factor(Habitat,
                          levels = c("phyllosphere", "soil", "freshwater sediment",
                                     "freshwater water", "marine sediment", "marine water",
                                     "hydrothermal vent", "wastewater", "mammalian gut",
                                     "human oral", "human skin", "human vaginal",
                                     "glacier")))
length(unique(phylum_counts$Phylum))
summary(phylum_counts$Count) # Range is 0 - 71
hist(phylum_counts$Count)
fig3 <- ggplot(phylum_counts, aes(Habitat, Phylum, fill = Count)) +
  geom_tile(colour = "grey50", linewidth = 0.1) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradientn(colours = c("white", viridisLite::viridis(100)),
                       values = scales::rescale(sqrt(c(0, seq(1, 66, length.out = 100)))),
                       limits = c(0, 66),
                       breaks = c(1, 10, 20, 30, 40, 50, 60, 66)) +
  labs(fill = "Indicator\ngenera") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8, margin = margin(r = 1), vjust = 0.4),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())
pdf("Manuscript/Figure3.pdf", width = 6, height = 8)
fig3
dev.off()
png("Manuscript/Figure3.png", width = 6, height = 8, units = "in", res = 300)
fig3
dev.off()



#### 3. Lake Tjörnin ####
# Use v1.0.0
taxonomy_loaded <- read.delim("taxonomy.tsv") %>%
  dplyr::select(domain, phylum, class, order, family, genus, group) %>%
  mutate(ID = genus) %>%
  column_to_rownames(var = "ID") %>%
  left_join(., InDB_formerge, by = c("genus" = "Genus")) %>%
  rename(taxonomy8 = Habitat) %>%
  mutate(taxonomy8 = ifelse(is.na(taxonomy8), "non-indicator", taxonomy8)) %>%
  mutate(OTU_ID = genus) %>%
  column_to_rownames(var = "OTU_ID") %>%
  set_names(c("taxonomy1", "taxonomy2", "taxonomy3", "taxonomy4", "taxonomy5",
              "taxonomy6", "taxonomy7", "taxonomy8"))
data_loaded <- read.delim("otu_table.tsv") %>%
  column_to_rownames(var = "X.OTU.ID")
map_loaded <- data.frame(sampleID = colnames(data_loaded)) %>%
  mutate(ID = sampleID) %>%
  column_to_rownames(var = "ID") %>%
  mutate(Year = gsub("X|\\.0", "", sampleID)) %>%
  mutate(Year = as.numeric(Year))
sum(rownames(taxonomy_loaded) != rownames(data_loaded))
sum(rownames(map_loaded) != colnames(data_loaded))
input <- list()
input$map_loaded <- map_loaded
input$data_loaded <- data_loaded
input$taxonomy_loaded <- taxonomy_loaded
sort(rowSums(input$data_loaded))
sort(colSums(input$data_loaded))
# Range from 1,404,437 to 88,497,629 (1.4 million to 88.4 million)
# Don't need to rarefy though
# We can just look at relative abundance without rarefying
# But, let's just do the relative abundance out of prokaryotes, not total
prok <- taxonomy_loaded %>%
  filter(taxonomy7 == "prokaryotes")
input_prok <- filter_taxa_from_input(input = input,
                                     taxa_IDs_to_keep = prok$taxonomy6)

# Basic info
nrow(input_prok$taxonomy_loaded) # 1285 prok genera in this dataset
table(input_prok$taxonomy_loaded$taxonomy1) # 1155 bac, 117 arc
ind_count <- as.data.frame(table(input_prok$taxonomy_loaded$taxonomy8)) 

# Summarize
tax_sum_ind <- summarize_taxonomy(input = input_prok,
                                  level = 8,
                                  relative = TRUE,
                                  report_higher_tax = FALSE)
colSums(tax_sum_ind) # Add up to 1 (not 100). These are proportions not percents.

meta_sort <- input_prok$map_loaded %>%
  mutate(Year = as.numeric(Year)) %>%
  arrange(Year)

# Default plot
ggplotly(plot_taxa_bars(tax_table = tax_sum_ind,
                        metadata_map = input_prok$map_loaded,
                        type_header = "Year",
                        num_taxa = 15,
                        data_only = FALSE) +
           theme_classic() +
           theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)))

plot_data <- plot_taxa_bars(tax_table = tax_sum_ind,
                            metadata_map = input_prok$map_loaded,
                            type_header = "Year",
                            num_taxa = 15,
                            data_only = TRUE)

# Check max abund of habitats
ind_max <- plot_data %>%
  group_by(taxon) %>%
  summarise(max = max(mean_value)*100,
            mean = mean(mean_value)*100) %>%
  ungroup() %>%
  arrange(desc(mean))

plot_data_abund <- plot_data %>%
  filter(taxon %in% c("freshwater sediment", "wastewater", "marine sediment",
                      "freshwater water", "soil", "marine water", "phyllosphere",
                      "glacier")) %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = c("freshwater sediment", "freshwater water", 
                                   "marine sediment", "marine water",
                                   "glacier", "phyllosphere",
                                   "soil", "wastewater")))
facet_text <- c(
  "freshwater sediment" = "freshwater\nsediment",
  "freshwater water" = "freshwater\nwater",
  "marine sediment" = "marine\nsediment",
  "marine water" = "marine\nwater",
  "glacier" = "glacier",
  "phyllosphere" = "phyllosphere",
  "soil" = "soil",
  "wastewater" = "wastewater"
)
epi <- data.frame(y = c(600, 1050, 1500,1950),
                  x = rep(0.45, 4),
                  taxon = rep("marine water", 4),
                  label = c("Episode I", "Episode II", "Episode III", "Episode IV"))
fig4 <- ggplot(plot_data_abund, aes(group_by, mean_value*100, fill = taxon)) +
  geom_area() +
  geom_line(linewidth = 0.1) +
  geom_vline(xintercept = 1703, linetype = "dashed") +
  geom_vline(xintercept = 1220, linetype = "dashed") +
  geom_vline(xintercept = 810, linetype = "dashed") +
  geom_text(data = epi, aes(x = y, y = x, label = label)) +
  labs(x = "Year (CE)",
       y = "% indicator genera abundance") +
  facet_wrap(~factor(taxon, levels = levels(plot_data_abund$taxon)), 
             nrow = 1, scales = "free_x", labeller = as_labeller(facet_text)) +
  scale_fill_manual(values = InDB_cols) +
  scale_y_continuous(breaks = pretty_breaks(n = 2)) +
  scale_x_continuous(breaks = seq(from = 200, to = 2000, by = 100),
                     limits = c(200, 2020)) +
  
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none",
        panel.grid.minor.y = element_blank())
fig4
pdf("Manuscript/Figure4.pdf", width = 8, height = 4)
fig4
dev.off()
png("Manuscript/Figure4.png", width = 8, height = 4, res = 300, units = "in")
fig4
dev.off()



#### 4. Dust ####
# Use v2.1.0

# Downloaded the repset and OTU table from Figshare
# Did not rerun DADA2 from raw data
# Just imported those 2 files and reassigned taxonomy with SBDI GTDB
# Import mctoolsr
tax_table_fp <- "dust_seqtab_wTax232_mctoolsr.txt"
map_fp <- "dust_metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 2854 samples loaded
nrow(input$taxonomy_loaded) # 115241
sum(colSums(input$data_loaded)) # 286,358,199 total seqs

# Already filtered out Chloroplast and Mitochondria with SILVA
# Did this immediately after assignTaxonomy, before saving the seqtab_wTax
# Note - could also use SILVA to filter out Euk.
# But these are few and probably Domain NA in GTDB assignment anyways

# Filter out Domain NA and Eukaryota
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 230 removed

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.)) # 0 because already removed on server prior to assignTaxonomy()

# Rarefaction?
# Check all samples including NTC, EB
# NTC and EB have very few seqs, good.
# rarecurve(t(input_filt$data_loaded), step = 10000, label = F)
# abline(v = c(50000), col = "red")
# hist(colSums(input_filt$data_loaded))
# sort(colSums(input_filt$data_loaded))

# Hmm this is not even close to saturating.
# Since I don't care about richness, let's not rarefy
# Can also check the blanks for human-associated indicators!
# But, still should have a very basic minimum cutoff - remove samples < 1200 reads
min_reads <- data.frame(sampleID = rownames(input_filt$map_loaded),
                        reads = colSums(input_filt$data_loaded)) %>%
  filter(reads > 1200) # 2732 of 2854
input_filt$map_loaded$sampleID <- rownames(input_filt$map_loaded)
input_filt <- filter_data(input_filt,
                          filter_cat = "sampleID",
                          keep_vals = min_reads$sampleID) # 2732 remaining
table(input_filt$map_loaded$Type)
table(input_filt$map_loaded$Location)
# Get only sample numbers with both I and O samples (i.e., n = 2)
paired_samples <- input_filt$map_loaded %>%
  group_by(sampleNum) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(count == 2) %>%
  dplyr::select(sampleNum)
control_samples <- input_filt$map_loaded %>%
  filter(grepl("NTC|Ext", sampleID)) %>%
  dplyr::select(sampleNum)
toKeep <- rbind(paired_samples, control_samples) # 1335 (1276 + 59)
input_filt <- filter_data(input_filt,
                          filter_cat = "sampleNum",
                          keep_vals = toKeep$sampleNum) # 2615

# Now analyze
# First make indic column in taxonomy loaded
# But there's an issue where species is assigned but genus is NA
input_filt$taxonomy_loaded <- input_filt$taxonomy_loaded %>%
  mutate(ID = taxonomy8) %>% # Use this to restore rownames after left_join
  separate(taxonomy7, into = c("Genus", "Species"), remove = F, sep = " ") %>%
  mutate(taxonomy6 = ifelse(taxonomy6 == "NA", Genus, taxonomy6)) %>%
  dplyr::select(-Genus, -Species) %>%
  left_join(., InDBv2_formerge, by = c("taxonomy6" = "Genus")) %>%
  rename(taxonomy9 = Habitat) %>%
  mutate(taxonomy9 = ifelse(taxonomy6 == "NA", "NA", taxonomy9)) %>%
  mutate(taxonomy9 = ifelse(is.na(taxonomy9), "non-indicator", taxonomy9)) %>%
  column_to_rownames(var = "ID")

# Now summarize by habitat and plot indoor vs outdoor (also check controls)
# Basic info
nrow(input_filt$taxonomy_loaded) # 107825 OTUs
table(input_filt$taxonomy_loaded$taxonomy1) # 107718 bac, 107 arc
ind_count <- as.data.frame(table(input_filt$taxonomy_loaded$taxonomy9))

# Summarize
tax_sum_gen <- summarize_taxonomy(input = input_filt,
                                  level = 6,
                                  relative = TRUE,
                                  report_higher_tax = FALSE)
ind_gen <- tax_sum_gen %>%
  filter(rownames(.) %in% InDBv2_formerge$Genus) %>%
  rownames_to_column(var = "Genus") %>%
  left_join(., InDBv2_formerge, by = "Genus")
table(ind_gen$Habitat)

tax_sum_ind <- summarize_taxonomy(input = input_filt,
                                  level = 9,
                                  relative = TRUE,
                                  report_higher_tax = FALSE)
colSums(tax_sum_ind) # Add up to 1 (not 100). These are proportions not percents.

# Default plot
table(input_filt$map_loaded$Type)
table(input_filt$map_loaded$Location)
# Need to make new varaible with NTC, ExtB, Indoor, Outdoor
input_filt$map_loaded$Treatment <- 
  ifelse(grepl("NTC", input_filt$map_loaded$sampleID), "NTC",
         ifelse(grepl("Ext", input_filt$map_loaded$sampleID), "EB",
                ifelse(input_filt$map_loaded$Location == "I", "Indoor", "Outdoor")))
table(input_filt$map_loaded$Treatment)
ggplotly(plot_taxa_bars(tax_table = tax_sum_ind,
                        metadata_map = input_filt$map_loaded,
                        type_header = "Treatment",
                        num_taxa = 15,
                        data_only = FALSE) +
           theme_classic() +
           theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)))

plot_data <- plot_taxa_bars(tax_table = tax_sum_ind,
                            metadata_map = input_filt$map_loaded,
                            type_header = "sampleID",
                            num_taxa = 15,
                            data_only = TRUE) %>%
  left_join(., input_filt$map_loaded, by = c("group_by" = "sampleID"))

# Now focus on samples
plot_data <- plot_data %>%
  filter(Treatment %in% c("Indoor", "Outdoor")) %>%
  droplevels()
levels(plot_data$Treatment) # It was character anyways

# Check mean and max abund of habitats
ind_max <- plot_data %>%
  group_by(taxon) %>%
  summarise(max = max(mean_value)*100,
            mean = mean(mean_value)*100) %>%
  ungroup() %>%
  arrange(desc(mean))
ind_max
# Note: All 13 habitats detected! Some veryyy rare though.
# Let's plot 5 (with minimum mean of 0.5%)

plot_data_abund <- plot_data %>%
  filter(taxon %in% c("phyllosphere", "soil",
                      "human skin", "human oral", "mammalian gut")) %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = c("phyllosphere", "soil",
                                   "human skin", "human oral", "mammalian gut")))
phy <- plot_data_abund %>%
  filter(taxon == "phyllosphere")
median(phy$mean_value, na.rm = T)

# Stats - need paired t-test or LME model with sampleNum as random
phyI <- plot_data_abund %>%
  filter(taxon == "phyllosphere") %>%
  filter(Treatment == "Indoor")
phyO <- plot_data_abund %>%
  filter(taxon == "phyllosphere") %>%
  filter(Treatment == "Outdoor")
t.test(phyI$mean_value, phyO$mean_value, paired = TRUE) # p < 0.001
wilcox.test(phyI$mean_value, phyO$mean_value, paired = TRUE) # p < 0.001

soiI <- plot_data_abund %>%
  filter(taxon == "soil") %>%
  filter(Treatment == "Indoor")
soiO <- plot_data_abund %>%
  filter(taxon == "soil") %>%
  filter(Treatment == "Outdoor")
t.test(soiI$mean_value, soiO$mean_value, paired = TRUE) # p < 0.001
wilcox.test(soiI$mean_value, soiO$mean_value, paired = TRUE) # p < 0.001

skiI <- plot_data_abund %>%
  filter(taxon == "human skin") %>%
  filter(Treatment == "Indoor")
skiO <- plot_data_abund %>%
  filter(taxon == "human skin") %>%
  filter(Treatment == "Outdoor")
t.test(skiI$mean_value, skiO$mean_value, paired = TRUE) # p < 0.001
wilcox.test(skiI$mean_value, skiO$mean_value, paired = TRUE) # p < 0.001

oraI <- plot_data_abund %>%
  filter(taxon == "human oral") %>%
  filter(Treatment == "Indoor")
oraO <- plot_data_abund %>%
  filter(taxon == "human oral") %>%
  filter(Treatment == "Outdoor")
t.test(oraI$mean_value, oraO$mean_value, paired = TRUE) # p < 0.001
wilcox.test(oraI$mean_value, oraO$mean_value, paired = TRUE) # p < 0.001

# wasI <- plot_data_abund %>%
#   filter(taxon == "wastewater") %>%
#   filter(Treatment == "Indoor")
# wasO <- plot_data_abund %>%
#   filter(taxon == "wastewater") %>%
#   filter(Treatment == "Outdoor")
# t.test(wasI$mean_value, wasO$mean_value, paired = TRUE) # p < 0.001
# wilcox.test(wasI$mean_value, wasO$mean_value, paired = TRUE) # p < 0.001

mamI <- plot_data_abund %>%
  filter(taxon == "human oral") %>%
  filter(Treatment == "Indoor")
mamO <- plot_data_abund %>%
  filter(taxon == "human oral") %>%
  filter(Treatment == "Outdoor")
t.test(mamI$mean_value, mamO$mean_value, paired = TRUE) # p < 0.001
wilcox.test(mamI$mean_value, mamO$mean_value, paired = TRUE) # p < 0.001

facet_text <- c(
  "phyllosphere" = "phyllosphere",
  "soil" = "soil",
  "human skin" = "human\nskin", 
  "human oral" = "human\noral",
  #"wastewater" = "wastewater",
  "mammalian gut" = "mammalian\ngut"
)
ind_max
label_df <- data.frame(taxon = c("phyllosphere", "soil",
                                 "human skin", "human oral", 
                                 "mammalian gut"),
                       x = c(1.5, 1.5, 1.5, 1.5, 1.5),
                       y = c(55, 21, 22, 50, 25),
                       label = c("***", "***",
                                 "***", "***", "***"))
fig5 <- ggplot(plot_data_abund, aes(Treatment, mean_value*100)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 1, pch = 16, width = 0.25, alpha = 0.25, 
              aes(colour = taxon)) +
  geom_text(data = label_df,
            aes(x = x, y = y, label = label)) +
  labs(x = "Location",
       y = "% indicator genera abundance") +
  facet_wrap(~factor(taxon, levels = levels(plot_data_abund$taxon)), 
             nrow = 1, labeller = as_labeller(facet_text),
             scales = "free_y") +
  scale_colour_manual(values = InDB_cols) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.grid.minor.y = element_blank())
fig5
pdf("Manuscript/Figure5.pdf", width = 8, height = 4)
fig5
dev.off()
png("Manuscript/Figure5.png", width = 8, height = 4, res = 300, units = "in")
fig5
dev.off()



#### 5. Willamette ####
# Use v2.1.0

# Make metadata
river_SRA <- read.csv("River_SraRunTable.csv") %>%
  filter(keep == "keep") %>%
  rename(sampleID = Run)
#write.table(river_SRA, "river_metadata232.txt", row.names = F, quote = F, sep = "\t")

# Import mctoolsr
tax_table_fp <- "river_seqtab_wTax232_mctoolsr.txt"
map_fp <- "river_metadata232.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 69 samples loaded
check <- river_SRA %>%
  filter(sampleID %notin% rownames(input$map_loaded)) # Missed that one. Oh well. Same coords as others.
nrow(input$taxonomy_loaded) # 10363
sum(colSums(input$data_loaded)) # 881940 total seqs

# Already filtered chloroplasts and mitochondria
# Here check for Domain NA and Euk
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 306 removed
# input_filt <- filter_taxa_from_input(input_filt,
#                                      taxa_to_remove = "Eukaryota",
#                                      at_spec_level = 1) # Not found

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_IDs_to_remove = singdoub$ASV) # 913 removed

# Rarefaction
hist(colSums(input_filt$data_loaded))
sort(colSums(input_filt$data_loaded))
rarecurve(t(input_filt$data_loaded), step = 1000, label = F)
abline(v = c(2757), col = "red")
abline(v = c(3484), col = "red")
abline(v = c(4148), col = "red")
abline(v = c(5348), col = "red")
set.seed(100)
input_filt_rare <- single_rarefy(input_filt, 5348) # n = 63
sort(colSums(input_filt_rare$data_loaded))
nrow(input_filt_rare$taxonomy_loaded) # 8304

# Add rarefied richness and Shannon
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, 
                                              MARGIN = 2)
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", 
                                                MARGIN = 2)

# Some locations have replicates at the same exact sampling location
# Group by coords and just use 1 sample per coordinate
uniq <- input_filt_rare$map_loaded %>%
  group_by(lat_lon, long) %>%
  slice_head(n = 1) %>%
  ungroup()
input_filt_rare <- filter_data(input = input_filt_rare,
                               filter_cat = "Sample.Name",
                               keep_vals = uniq$Sample.Name) # n = 38

# Now analyze
# First make indic column in taxonomy loaded
# But there's an issue where species is assigned but genus is NA
input_filt_rare$taxonomy_loaded <- input_filt_rare$taxonomy_loaded %>%
  mutate(ID = taxonomy8) %>% # Use this to restore rownames after left_join
  separate(taxonomy7, into = c("Genus", "Species"), remove = F, sep = " ") %>%
  mutate(taxonomy6 = ifelse(taxonomy6 == "NA", Genus, taxonomy6)) %>%
  dplyr::select(-Genus, -Species) %>%
  left_join(., InDBv2_formerge, by = c("taxonomy6" = "Genus")) %>%
  rename(taxonomy9 = Habitat) %>%
  mutate(taxonomy9 = ifelse(taxonomy6 == "NA", "NA", taxonomy9)) %>%
  mutate(taxonomy9 = ifelse(is.na(taxonomy9), "non-indicator", taxonomy9)) %>%
  column_to_rownames(var = "ID")

# Now summarize by habitat and plot over the watershed
# Basic info
nrow(input_filt_rare$taxonomy_loaded) # 7144 ASVs
table(input_filt_rare$taxonomy_loaded$taxonomy1) # 7028 bac, 116 arc
ind_count <- as.data.frame(table(input_filt_rare$taxonomy_loaded$taxonomy9))

# Summarize
tax_sum_gen <- summarize_taxonomy(input = input_filt_rare,
                                  level = 6,
                                  relative = TRUE,
                                  report_higher_tax = FALSE)
ind_gen <- tax_sum_gen %>%
  filter(rownames(.) %in% InDBv2_formerge$Genus) %>%
  rownames_to_column(var = "Genus") %>%
  left_join(., InDBv2_formerge, by = "Genus")
table(ind_gen$Habitat)
tax_sum_ind <- summarize_taxonomy(input = input_filt_rare,
                                  level = 9,
                                  relative = TRUE,
                                  report_higher_tax = FALSE)
colSums(tax_sum_ind) # Add up to 1 (not 100). These are proportions not percents.

# Default plot
ggplotly(plot_taxa_bars(tax_table = tax_sum_ind,
                        metadata_map = input_filt_rare$map_loaded,
                        type_header = "Sample.Name",
                        num_taxa = 15,
                        data_only = FALSE) +
           theme_classic() +
           theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)))

input_filt_rare$map_loaded$sampleID <- rownames(input_filt_rare$map_loaded)
plot_data <- plot_taxa_bars(tax_table = tax_sum_ind,
                            metadata_map = input_filt_rare$map_loaded,
                            type_header = "sampleID",
                            num_taxa = 15,
                            data_only = TRUE) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))

# Check mean and max abund of habitats
ind_max <- plot_data %>%
  group_by(taxon) %>%
  summarise(max = max(mean_value)*100,
            mean = mean(mean_value)*100) %>%
  ungroup() %>%
  arrange(desc(mean))
ind_max
# Note: 0 human oral, human skin, human vaginal detected

plot_data_abund <- plot_data %>%
  filter(taxon %in% c("freshwater water", "freshwater sediment", "glacier", "wastewater",
                      "soil")) %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = c("freshwater water", "freshwater sediment",
                                   "glacier", "soil", "wastewater")))

fww <- plot_data_abund %>%
  filter(taxon == "freshwater water")
m <- aov(mean_value ~ Location.1, data = fww)
summary(m)
Anova(m)
TukeyHSD(m) # a b b
fws <- plot_data_abund %>%
  filter(taxon == "freshwater sediment")
m <- aov(mean_value ~ Location.1, data = fws)
summary(m)
Anova(m)
TukeyHSD(m) # a b ab
gla <- plot_data_abund %>%
  filter(taxon == "glacier")
m <- aov(mean_value ~ Location.1, data = gla)
summary(m)
Anova(m)
TukeyHSD(m) # a b ab
was <- plot_data_abund %>%
  filter(taxon == "wastewater")
m <- aov(mean_value ~ Location.1, data = was)
summary(m)
Anova(m)
TukeyHSD(m) # N.S.D
soi <- plot_data_abund %>%
  filter(taxon == "soil")
m <- aov(mean_value ~ Location.1, data = soi)
summary(m)
Anova(m)
TukeyHSD(m) # a b b

label <- data.frame(taxon = c("freshwater water", "freshwater water", "freshwater water",
                              "freshwater sediment", "freshwater sediment", "freshwater sediment",
                              "glacier", "glacier", "glacier",
                              "wastewater", "wastewater", "wastewater",
                              "soil", "soil", "soil"),
                    Location.1 = c("Mountains", "Plains", "Urban", 
                                   "Mountains", "Plains", "Urban", 
                                   "Mountains", "Plains", "Urban", 
                                   "Mountains", "Plains", "Urban", 
                                   "Mountains", "Plains", "Urban"),
                    y = c(31,31,31, 1.75,1.75,1.75, 14,14,14, 0.75,0.75,0.75, 3,3,3),
                    label = c("a", "b", "b",
                              "a", "b", "ab",
                              "", "N.S.D.", "",
                              "", "N.S.D.", "",
                              "a", "b", "b"))
facet_text <- c(
  "freshwater water" = "freshwater\nwater",
  "freshwater sediment" = "freshwater\nsediment",
  "glacier" = "glacier",
  "wastewater" = "wastewater",
  "soil" = "soil"
)
fig6 <- ggplot(plot_data_abund, aes(Location.1, mean_value*100)) +
  geom_boxplot(outliers = F) +
  geom_jitter(pch = 21, size = 2, aes(fill = taxon), width = 0.25, height = 0) +
  geom_text(data = label, aes(x = Location.1, y = y, label = label)) +
  labs(x = NULL,
       y = "% indicator genera abundance") +
  facet_wrap(~factor(taxon, levels = levels(plot_data_abund$taxon)),
             nrow = 1, labeller = as_labeller(facet_text),
             scales = "free_y") +
  scale_fill_manual(values = InDB_cols) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.grid.minor.y = element_blank())
pdf("Manuscript/Figure6.pdf", width = 8, height = 4)
fig6
dev.off()
png("Manuscript/Figure6.png", width = 8, height = 4, res = 300, units = "in")
fig6
dev.off()



