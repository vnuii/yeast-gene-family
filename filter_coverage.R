# This script performs coverage filtering on gene count data for different clades (yeast, pezizomycotina, animal, and plant).
# It reads input files containing information about the clades and their corresponding gene count data.
# The script calculates the mean coverage and gene number for each orthogroup within each clade.
# It then filters the orthogroups based on a coverage threshold (0.1 and 0.5) and writes the filtered data to separate output files.
library(data.table)
library(dplyr)

setwd("PATH/TO/THIS/DOWNLOAD/ROOT/DIRACTORY")

# File paths for input data
info_yeast_path <- "./info_yeast.tsv"
yeast_og_count_path <- "./orthofinder_results/yeast_Orthogroups.GeneCount.tsv"
info_pezizo_path <- "./info_pezizomycotina.tsv"
pezizo_og_count_path <- "./orthofinder_results/pezizomycotina_Orthogroups.GeneCount.tsv"
info_animal_path <- "./info_animal.tsv"
animal_og_count_path <- "./orthofinder_results/animal_Orthogroups.GeneCount.tsv"
info_plant_path <- "./info_plant.tsv"
plant_og_count_path <- "./orthofinder_results/plant_Orthogroups.GeneCount.tsv"

## yeast
# Read and process information about yeast clade
info_yeast <- fread(info_yeast_path) %>%
  select("assembly_fullID_updated", "Order") %>%
  setnames(c("id", "clade")) %>%
  setkey("clade")

# Read yeast orthogroup count data
yeast_og_count <- fread(yeast_og_count_path, drop = "Total")

# Create a table to store coverage information for yeast orthogroups
covertable_yeast <- data.table(
  Orthogroup = yeast_og_count$Orthogroup,
  mean_coverage = NA,
  genenumber = NA
)

# Calculate coverage for each clade within yeast
for (clade in unique(info_yeast$clade)) {
  covertable_yeast[, paste0(clade, "_coverage") :=
    rowSums(yeast_og_count[, info_yeast[clade]$id,
      with = FALSE
    ] > 0) /
      length(info_yeast[clade]$id)]
}

# Calculate mean coverage and gene number for yeast orthogroups
covertable_yeast[, mean_coverage :=
  rowMeans(
    covertable_yeast[, 2:ncol(covertable_yeast),
      with = FALSE
    ],
    na.rm = TRUE
  )]

covertable_yeast[, gene_number := rowSums(yeast_og_count[, -1, with = FALSE])]

# Filter yeast orthogroups based on coverage threshold (0.1 and 0.5) and write to output files
og_0.1_yeast <- covertable_yeast[mean_coverage >= 0.1]$Orthogroup
og_0.5_yeast <- covertable_yeast[mean_coverage >= 0.5]$Orthogroup

fwrite(covertable_yeast,
  "./filtering_coverage/output/yeast_coverage.tsv",
  sep = "\t"
)
fwrite(yeast_og_count[Orthogroup %in% og_0.1_yeast],
  "./orthofinder_results/yeast_0.1_Orthogroups.GeneCount.tsv",
  sep = "\t"
)
fwrite(yeast_og_count[Orthogroup %in% og_0.5_yeast],
  "./orthofinder_results/yeast_0.5_Orthogroups.GeneCount.tsv",
  sep = "\t"
)

## pezizomycotina
# Read and process information about pezizomycotina clade
info_pezizo <- fread(info_pezizo_path) %>%
  select("old_taxonID_linked_genome_sequence", "clade/class") %>%
  setnames(c("id", "clade")) %>%
  setkey("clade")

# Read pezizomycotina orthogroup count data
pezizo_og_count <- fread(pezizo_og_count_path, drop = "Total")

# Create a table to store coverage information for pezizomycotina orthogroups
covertable_pezizo <- data.table(
  Orthogroup = pezizo_og_count$Orthogroup,
  mean_coverage = NA,
  genenumber = NA
)

# Calculate coverage for each clade within pezizomycotina
for (clade in unique(info_pezizo$clade)) {
  covertable_pezizo[, paste0(clade, "_coverage") :=
    rowSums(pezizo_og_count[, info_pezizo[clade]$id,
      with = FALSE
    ] > 0) /
      length(info_pezizo[clade]$id)]
}

# Calculate mean coverage and gene number for pezizomycotina orthogroups
covertable_pezizo[, mean_coverage :=
  rowMeans(
    covertable_pezizo[, 2:ncol(covertable_pezizo),
      with = FALSE
    ],
    na.rm = TRUE
  )]

covertable_pezizo[, gene_number := rowSums(pezizo_og_count[, -1, with = FALSE])]

# Filter pezizomycotina orthogroups based on coverage threshold (0.1) and write to output files
og_0.1_pezizo <- covertable_pezizo[mean_coverage >= 0.1]$Orthogroup

fwrite(covertable_pezizo,
  "./filtering_coverage/output/pezizo_coverage.tsv",
  sep = "\t"
)
fwrite(pezizo_og_count[Orthogroup %in% og_0.1_pezizo],
  "./orthofinder_results/pezizo_0.1_Orthogroups.GeneCount.tsv",
  sep = "\t"
)

## animal
# Read and process information about animal clade
info_animal <- fread(info_animal_path) %>%
  select("Species", "Phylum") %>%
  setnames(c("id", "clade")) %>%
  setkey("clade")

# Read animal orthogroup count data
animal_og_count <- fread(animal_og_count_path, drop = "Total")

# Create a table to store coverage information for animal orthogroups
covertable_animal <- data.table(
  Orthogroup = animal_og_count$Orthogroup,
  mean_coverage = NA,
  genenumber = NA
)

# Calculate coverage for each clade within animal
for (clade in unique(info_animal$clade)) {
  covertable_animal[, paste0(clade, "_coverage") :=
    rowSums(animal_og_count[, info_animal[clade]$id,
      with = FALSE
    ] > 0) /
      length(info_animal[clade]$id)]
}

# Calculate mean coverage and gene number for animal orthogroups
covertable_animal[, mean_coverage :=
  rowMeans(
    covertable_animal[, 2:ncol(covertable_animal),
      with = FALSE
    ],
    na.rm = TRUE
  )]

covertable_animal[, gene_number := rowSums(animal_og_count[, -1, with = FALSE])]

# Filter animal orthogroups based on coverage threshold (0.1) and write to output files
og_0.1_animal <- covertable_animal[mean_coverage >= 0.1]$Orthogroup

fwrite(covertable_animal,
  "./filtering_coverage/output/animal_coverage.tsv",
  sep = "\t"
)
fwrite(animal_og_count[Orthogroup %in% og_0.1_animal],
  "./orthofinder_results/animal_0.1_Orthogroups.GeneCount.tsv",
  sep = "\t"
)

## plant
# Read and process information about plant clade
info_plant <- fread(info_plant_path) %>%
  select("1KP Index ID", "Clade") %>%
  setnames(c("id", "clade")) %>%
  setkey("clade")

# Read plant orthogroup count data
plant_og_count <- fread(plant_og_count_path, drop = "Total")

# Create a table to store coverage information for plant orthogroups
covertable_plant <- data.table(
  Orthogroup = plant_og_count$Orthogroup,
  mean_coverage = NA,
  genenumber = NA
)

# Calculate coverage for each clade within plant
for (clade in unique(info_plant$clade)) {
  covertable_plant[, paste0(clade, "_coverage") :=
    rowSums(plant_og_count[, info_plant[clade]$id,
      with = FALSE
    ] > 0) /
      length(info_plant[clade]$id)]
}

# Calculate mean coverage and gene number for plant orthogroups
covertable_plant[, mean_coverage :=
  rowMeans(
    covertable_plant[, 2:ncol(covertable_plant),
      with = FALSE
    ],
    na.rm = TRUE
  )]

covertable_plant[, gene_number := rowSums(plant_og_count[, -1, with = FALSE])]

# Filter plant orthogroups based on coverage threshold (0.1) and write to output files
og_0.1_plant <- covertable_plant[mean_coverage >= 0.1]$Orthogroup

fwrite(covertable_plant,
  "./filtering_coverage/output/plant_coverage.tsv",
  sep = "\t"
)
fwrite(plant_og_count[Orthogroup %in% og_0.1_plant],
  "./orthofinder_results/plant_0.1_Orthogroups.GeneCount.tsv",
  sep = "\t"
)
