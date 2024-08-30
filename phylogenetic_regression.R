# Load required libraries
library(data.table)
library(dplyr)
library(ape)

setwd("/home/vonui/KEGG_analysis/figshare_data_figshare")

# Read the input data files
genecount_wsize <- fread("./phylogenetic_regression/wsize_genecount_4groups.tsv")
genomesize_dt <- fread("./phylogenetic_regression/genomesize_4groups.tsv")

# Perform analysis for yeast group
# Read yeast tree and information file
yeast_tree <- read.tree("./trees_4groups/yeast_tree.txt")
yeast_info_original <- fread("./phylogenetic_regression/info_yeast_original.tsv")

# Filter yeast genome data based on tree and information file
yeast_genome_dt <- genomesize_dt[group == "yeast"] %>%
  .[, Tip_id := yeast_info_original$Tip_ID[match(ID, yeast_info_original$assembly_fullID_updated)]] %>%
  .[match(yeast_tree$tip.label, Tip_id)]

# Calculate PIC (phylogenetic independent contrast) for yeast genomesize
yeast_pic_genomesize <- pic(phy = yeast_tree, yeast_genome_dt$Genome_size)

# Filter yeast gene count and weighted size data based on tree and information file
yeast_dt <- genecount_wsize[Group == "yeast"] %>%
  .[, Tip_id := yeast_info_original$Tip_ID[match(ID, yeast_info_original$assembly_fullID_updated)]] %>%
  .[match(yeast_tree$tip.label, Tip_id)]

# Calculate PIC for yeast gene count and weighted size
yeast_pic_genecount <- pic(phy = yeast_tree, yeast_dt$Gene_count)
yeast_pic_wsize <- pic(phy = yeast_tree, yeast_dt$Weighted_average_size)

# Perform correlation and linear regression analysis for yeast group
cor.test(yeast_pic_genomesize, yeast_pic_wsize, method = "spearman")
cor.test(yeast_pic_genecount, yeast_pic_wsize, method = "spearman")
yeast_slope <- lm(yeast_pic_wsize ~ yeast_pic_genecount + 0)

# Perform analysis for plant group
# Read plant tree and filter genome data based on tree
plant_tree <- read.tree("./trees_4groups/plant_tree.txt")
plant_tree_genome <- drop.tip(plant_tree, plant_tree$tip.label[!plant_tree$tip.label %in% genomesize_dt[group == "plant"]$ID])

# Filter plant gene count and weighted size data based on tree
plant_dt <- genecount_wsize[Group == "plant"] %>%
  .[match(plant_tree$tip.label, ID)]

# Filter plant genome data based on tree, only use the species that have genome data (31 genomes)
plant_genome_dt <- genomesize_dt[group == "plant"] %>%
  .[match(plant_tree_genome$tip.label, ID)]

# Calculate PIC for plant genomesize and weighted size
plant_pic_genomesize <- pic(phy = plant_tree_genome, plant_genome_dt$Genome_size)
plant_pic_wsize_genome <- pic(phy = plant_tree_genome, plant_genome_dt$Weighted_average_size)

# Calculate PIC for plant genecount and weighted size
plant_pic_genecount <- pic(phy = plant_tree, plant_dt$Gene_count)
plant_pic_wsize <- pic(phy = plant_tree, plant_dt$Weighted_average_size)

# Perform correlation analysis for plant group
cor.test(plant_pic_genomesize, plant_pic_wsize_genome, method = "spearman")
cor.test(plant_pic_genecount, plant_pic_wsize, method = "spearman")
lm(plant_pic_wsize ~ plant_pic_genecount - 1)

# Perform analysis for pezizo group
pezizo_tree <- read.tree("./trees_4groups/pezizo_tree.txt")

# Filter pezizo gene count and weighted size data based on tree
pezizo_dt <- genecount_wsize[Group == "pezizo"] %>%
  .[match(pezizo_tree$tip.label, ID)]

# Filter pezizo genome data based on tree
pezizo_genome_dt <- genomesize_dt[group == "pezizo"] %>%
  .[match(pezizo_tree$tip.label, ID)]

# Calculate PIC for pezizo genomesize, gene count, and weighted size
pezizo_pic_genomesize <- pic(phy = pezizo_tree, pezizo_genome_dt$Genome_size)
pezizo_pic_genecount <- pic(phy = pezizo_tree, pezizo_dt$Gene_count)
pezizo_pic_wsize <- pic(phy = pezizo_tree, pezizo_dt$Weighted_average_size)

# Perform correlation and linear regression analysis for pezizo group
cor.test(pezizo_pic_genomesize, pezizo_pic_wsize, method = "spearman")
cor.test(pezizo_pic_genecount, pezizo_pic_wsize, method = "spearman")
lm(pezizo_pic_wsize ~ pezizo_pic_genecount - 1)

# Perform analysis for animal group
animal_tree <- read.tree("./trees_4groups/animal_tree.txt")

# Filter animal gene count and weighted size data based on tree
animal_dt <- genecount_wsize[Group == "animal"] %>%
  .[match(animal_tree$tip.label, ID)]

# Filter animal genome data based on tree
animal_genome_dt <- genomesize_dt[group == "animal"] %>%
  .[match(animal_tree$tip.label, ID)]

# Calculate PIC for animal genomesize, gene count, and weighted size
animal_pic_genomesize <- pic(phy = animal_tree, animal_genome_dt$Genome_size)
animal_pic_genecount <- pic(phy = animal_tree, animal_dt$Gene_count)
animal_pic_wsize <- pic(phy = animal_tree, animal_dt$Weighted_average_size)

# Perform correlation and linear regression analysis for animal group
cor.test(animal_pic_genomesize, animal_pic_wsize, method = "spearman")
cor.test(animal_pic_genecount, animal_pic_wsize, method = "spearman")
lm(animal_pic_wsize ~ animal_pic_genecount - 1)

#############
plot_pic_dt <- data.table(
  pic_wsize = c(yeast_pic_wsize, plant_pic_wsize, pezizo_pic_wsize, animal_pic_wsize),
  pic_genecount = c(yeast_pic_genecount, plant_pic_genecount, pezizo_pic_genecount, animal_pic_genecount),
  group = c(rep("Yeasts", length(yeast_pic_wsize)), rep("Plants", length(plant_pic_wsize)), rep("Filamentous ascomycetes", length(pezizo_pic_wsize)), rep("Animals", length(animal_pic_wsize)))
)

library(ggplot2)

ggplot(plot_pic_dt, aes(x = pic_genecount, y = pic_wsize, color = group)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed") +
  scale_color_manual(values = c("Filamentous ascomycetes" = "#5B8FA8FF", "Yeasts" = "#800000FF", "Animals" = "#FFB547FF", "Plants" = "#ADB17DFF")) +
  labs(x = "PICs for Total number of genes", y = "PICs for Avergae gene family size", color = "") +
  # scale_x_continuous(labels = c("0", "10", "20", "30", "40", "50")) +
  theme_classic(base_size = 18, base_line_size = 0.5, base_family = "Arial") +
  theme(axis.title = element_text(face = "bold")) +
  theme(legend.position = "none") +
  facet_wrap(~group)
fwrite(plot_pic_dt, "/home/vonui/KEGG_analysis/figshare_data_figshare/phylogenetic_regression/PICs_wsize_genecount_4groups.tsv", sep = "\t")
