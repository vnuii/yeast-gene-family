library(data.table)
library(dplyr)
# also need to install stringr package

# Set the working directory to the specified path
setwd("/share/ME-58T/1kp/a_main_graph/figshare_data")

# Define the file paths for the input data
info_yeast_path <- "./info_yeast.tsv"
yeast_0.1_og_count_path <- "./orthofinder_results/yeast_0.1_Orthogroups.GeneCount.tsv"
yeast_0.5_og_count_path <- "./orthofinder_results/yeast_0.5_Orthogroups.GeneCount.tsv"
info_pezizo_path <- "./info_pezizomycotina.tsv"
pezizo_0.1_og_count_path <- "./orthofinder_results/pezizo_0.1_Orthogroups.GeneCount.tsv"
info_animal_path <- "./info_animal.tsv"
animal_0.1_og_count_path <- "./orthofinder_results/animal_0.1_Orthogroups.GeneCount.tsv"
info_plant_path <- "./info_plant.tsv"
plant_0.1_og_count_path <- "./orthofinder_results/plant_0.1_Orthogroups.GeneCount.tsv"

# Function to calculate the weighted average size
getWeightedAverageSize <- function(og_count_path, info_path, threshold) {
  # Read the orthogroup count data from the specified path
  og_count <- fread(og_count_path)

  # Determine the group based on the info path
  group <- gsub(".+_(\\w+)\\..+", "\\1", info_path)

  # Process the info data based on the group
  if (group == "yeast") {
    info_dt <- fread(info_path) %>%
      select("assembly_fullID_updated", "Order") %>%
      setnames(c("id", "clade"))
  }
  if (group == "pezizomycotina") {
    info_dt <- fread(info_path) %>%
      select("old_taxonID_linked_genome_sequence", "clade/class") %>%
      setnames(c("id", "clade"))
  }
  if (group == "animal") {
    info_dt <- fread(info_path) %>%
      select("Species", "Phylum") %>%
      setnames(c("id", "clade"))
  }
  if (group == "plant") {
    info_dt <- fread(info_path) %>%
      select("1KP Index ID", "Clade") %>%
      setnames(c("id", "clade"))
  }

  # Calculate the weights based on the maximum value in each row of the orthogroup count data
  weights <- 1 / apply(og_count[, -1, with = FALSE], 1, max)

  # Calculate the mean of the maximum values
  mean_max <- apply(og_count[, -1, with = FALSE], 1, max) %>% mean()

  # Calculate the weighted average size for each column of the orthogroup count data
  weighted_average_size <- apply(
    og_count[, -1], 2,
    function(x) {
      a <- x * weights
      mean(a) * mean_max
    }
  )

  # Create a data table with the column names and corresponding weighted average sizes
  dt <- data.table(
    id = colnames(og_count)[-1],
    weighted_average_size = weighted_average_size
  )

  # Add the clade information from the info data to the data table
  dt[, clade := info_dt$clade[match(dt$id, info_dt$id)]]

  # Calculate the mean weighted average size for each clade
  mean_clade <- dt %>%
    group_by(clade) %>%
    summarise(mean_weighted_average_size = mean(weighted_average_size)) %>%
    as.data.table()

  # Add the mean weighted average size to the data table based on the clade
  dt[, mean_clade := mean_clade$mean_weighted_average_size[match(dt$clade, mean_clade$clade)]]

  # Write the data table to a file with the specified name
  fwrite(dt, paste0("/share/ME-58T/1kp/a_main_graph/figshare_data/weighted_average_size/output/", group, "_", threshold, "_weighted_average_size.tsv"), sep = "\t")
}

# Call the getWeightedAverageSize function for different input data
getWeightedAverageSize(yeast_0.1_og_count_path, info_yeast_path, 0.1)
getWeightedAverageSize(yeast_0.5_og_count_path, info_yeast_path, 0.5)
getWeightedAverageSize(pezizo_0.1_og_count_path, info_pezizo_path, 0.1)
getWeightedAverageSize(animal_0.1_og_count_path, info_animal_path, 0.1)
getWeightedAverageSize(plant_0.1_og_count_path, info_plant_path, 0.1)
