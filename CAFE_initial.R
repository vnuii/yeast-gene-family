# This script initializes the files required for CAFE analysis.

# Load necessary libraries
library(data.table)
library(dplyr)
library(tibble)

# Set the directory paths
OG_counts_dir <- "./CAFE/1154_cladeCAFE/origin_og_table_dir"
destin_dir <- "./CAFE/1154_cladeCAFE"
tree_file_path <- "./CAFE/timetree_1154.nh"
gotree_path <- "~/anaconda3/envs/basic/bin/gotree"
id_tipID_dt <- fread("./CAFE/id_tipID.tsv")

# Create necessary directories
dirs_2make <- c("node_tree_file_dir", "parameter_og_table_dir", "sp_name_dir")
lapply(dirs_2make, function(x) {
  dir.create(file.path(destin_dir, x), showWarnings = FALSE)
})

# Process each OG counts file
OG_counts_files <- list.files(OG_counts_dir, full.names = TRUE)
lapply(OG_counts_files, function(x) {
  OG_counts <- fread(x)
  name_ <- basename(x) %>%
    gsub(".tsv", "", .)

  # Add a "Des" column if it doesn't exist
  if (!"Des" %in% colnames(OG_counts)) {
    OG_counts <- add_column(OG_counts, Des = OG_counts$Orthogroup, .before = 1)
  }

  ori_ids <- colnames(OG_counts)[-(1:2)]

  # Check if any original IDs are missing in the id_tipID_dt table
  if (!any(ori_ids %in% id_tipID_dt$timetree_id)) {
    tip_ids <- id_tipID_dt$timetree_id[match(ori_ids, id_tipID_dt$assembly_fullID_updated)]
    colnames(OG_counts)[-(1:2)] <- tip_ids

    # Write the modified OG counts to the destination directory
    fwrite(OG_counts, file.path(destin_dir, "origin_og_table_dir", paste0(name_, ".tsv")), sep = "\t")
  }

  # Write the OG counts to the destination directories
  fwrite(OG_counts, file.path(destin_dir, "origin_og_table_dir", paste0(name_, ".tsv")), sep = "\t")
  fwrite(OG_counts, file.path(destin_dir, "parameter_og_table_dir", paste0(name_, ".tsv")), sep = "\t")

  # Write the tip IDs to a file in the sp_name_dir
  tip_ids <- colnames(OG_counts)[-(1:2)]
  fwrite(list(tip_ids), file.path(destin_dir, "sp_name_dir", paste0(name_, ".txt")), sep = "\n")

  # Generate the gotree command and execute it
  gotree_command <- paste(
    gotree_path, "prune -i", tree_file_path,
    "-o", file.path(destin_dir, "node_tree_file_dir", paste0(name_, ".nh")),
    "-r", paste0("`cat ", file.path(destin_dir, "sp_name_dir", paste0(name_, ".txt")), "`")
  )
  system(gotree_command)
})
