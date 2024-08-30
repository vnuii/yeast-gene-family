# This script runs the FEL and SEL analyses using BayesTraits software for the Dipodascales dataset.
# It loads the necessary libraries, reads the carbon growth table, and defines the carbon traits names.
# Then, it runs the FEL and SEL analyses for each carbon trait, both in the free and restrict models.

# Required Libraries:
# - dplyr
# - data.table

# Required Files:
# - BayesTraitsV4: Path to the BayesTraits software executable
# - carbon_growth_table_1154.tsv: Path to the carbon growth table file
# - Dipodascales_FEL_tree_nexus.txt: Path to the FEL tree file
# - Dipodascales_SEL_tree_nexus.txt: Path to the SEL tree file

# FEL Analysis:
# - The FEL analysis is performed for each carbon trait in the carbon_traits_names_18 list.
# - For each trait, a command file and a traits file are generated based on the trait name.
# - The system command is used to run BayesTraits with the FEL analysis command file.

# SEL Analysis:
# - The SEL analysis is performed for each carbon trait in the carbon_traits_names_18 list.
# - For each trait, a command file and a traits file are generated based on the trait name.
# - The system command is used to run BayesTraits with the SEL analysis command file.

# Note: Make sure to update the file paths according to your specific directory structure.

library(dplyr)
library(data.table)

bayestraits <- "/home/vonui/KEGG_analysis/figshare_data_figshare/bayestraits/BayesTraitsV4.0.0-Linux-Threaded/BayesTraitsV4"
carbon_growth_table <- fread("/home/vonui/KEGG_analysis/figshare_data_figshare/bayestraits/carbon18_science_Dipodascales.tsv")
carbon_traits_names_18 <- colnames(carbon_growth_table)[-c(1:6)]

#### 跑Dipodascales的18个traits的FEL和SEL
# # FEL
# setwd("/home/vonui/KEGG_analysis/figshare_data_figshare/bayestraits/BayesTraitsV4.0.0-Linux-Threaded/Dipodascales_carbon_18traits/FEL/output")

# command_file_basepath <- "/home/vonui/KEGG_analysis/figshare_data_figshare/bayestraits/BayesTraitsV4.0.0-Linux-Threaded/Dipodascales_carbon_18traits/FEL/command_files/"
# traits_file_basepath <- "/home/vonui/KEGG_analysis/figshare_data_figshare/bayestraits/BayesTraitsV4.0.0-Linux-Threaded/Dipodascales_carbon_18traits/FEL/data/"

# # free
# for (i in carbon_traits_names_18) {
# 	command_file <- paste0(command_file_basepath, "command_Dipodascales_FEL_free_", i, ".txt")
# 	traits_file <- paste0(traits_file_basepath, "Dipodascales_FEL_carbontrait_", i, ".tsv")
# 	tree_file <- paste0(traits_file_basepath, "Dipodascales_FEL_tree_nexus_", i, ".txt")
# 	system(paste(bayestraits, tree_file, traits_file, "<", command_file))
# }

# # restrict
# for (i in carbon_traits_names_18) {
# 	command_file <- paste0(command_file_basepath, "command_Dipodascales_FEL_restrict_", i, ".txt")
# 	traits_file <- paste0(traits_file_basepath, "Dipodascales_FEL_carbontrait_", i, ".tsv")
# 	tree_file <- paste0(traits_file_basepath, "Dipodascales_FEL_tree_nexus_", i, ".txt")
# 	system(paste(bayestraits, tree_file, traits_file, "<", command_file))
# }


## SEL
setwd("/home/vonui/KEGG_analysis/figshare_data_figshare/bayestraits/BayesTraitsV4.0.0-Linux-Threaded/Dipodascales_carbon_18traits/SEL/output/")

command_file_basepath <- "/home/vonui/KEGG_analysis/figshare_data_figshare/bayestraits/BayesTraitsV4.0.0-Linux-Threaded/Dipodascales_carbon_18traits/SEL/command_files/"
traits_file_basepath <- "/home/vonui/KEGG_analysis/figshare_data_figshare/bayestraits/BayesTraitsV4.0.0-Linux-Threaded/Dipodascales_carbon_18traits/SEL/data/"

# free
for (i in carbon_traits_names_18) {
  command_file <- paste0(command_file_basepath, "command_Dipodascales_SEL_free_", i, ".txt")
  traits_file <- paste0(traits_file_basepath, "Dipodascales_SEL_carbontrait_", i, ".tsv")
  tree_file <- paste0(traits_file_basepath, "Dipodascales_SEL_tree_nexus_", i, ".txt")
  system(paste(bayestraits, tree_file, traits_file, "<", command_file))
}

# restrict
for (i in carbon_traits_names_18) {
  command_file <- paste0(command_file_basepath, "command_Dipodascales_SEL_restrict_", i, ".txt")
  traits_file <- paste0(traits_file_basepath, "Dipodascales_SEL_carbontrait_", i, ".tsv")
  tree_file <- paste0(traits_file_basepath, "Dipodascales_SEL_tree_nexus_", i, ".txt")
  system(paste(bayestraits, tree_file, traits_file, "<", command_file))
}
