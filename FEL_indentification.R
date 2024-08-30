library(data.table)
library(dplyr)
library(diptest)
library(ggplot2)
library(dbscan)

setwd("PATH/TO/THIS/DOWNLOAD/ROOT/DIRACTORY")

# Read the yeast weighted average size data from a TSV file
yeast_01_wsize_length <- fread("./FEL_identification/yeast_0.1_weighted_average_size.tsv")

# Perform multimodality test for each clade
pval_table <- data.table(
  clade = unique(yeast_01_wsize_length$clade),
  p_value_length = as.numeric(NA),
  p_value_wsize = as.numeric(NA)
)

for (i in pval_table$clade) {
  # Perform dip test on branch length
  pvalue_length <- yeast_01_wsize_length[clade == i, ]$branch_length %>%
    dip.test() %>%
    .$p.value
  # Perform dip test on weighted average size
  pvalue_wsize <- yeast_01_wsize_length[clade == i, ]$weighted_average_size %>%
    dip.test() %>%
    .$p.value
  # Store the p-values in the pval_table
  pval_table[clade == i, p_value_wsize := pvalue_wsize]
  pval_table[clade == i, p_value_length := pvalue_length]
}

# Sort the pval_table by descending p_value_length
pval_table <- arrange(pval_table, desc(p_value_length))

# Select the clades with p_value_length < 0.05
orders_interested <- pval_table[p_value_length < 0.05]$clade

# Write the pval_table to a TSV file
fwrite(pval_table, "/share/ME-58T/1kp/a_main_graph/figshare_data/FEL_identification/output/yeast_multimodality_table.tsv", sep = "\t")

# Perform clustering on the selected clades
plot_dbscan_table <- yeast_01_wsize_length[clade %in% orders_interested]
plot_dbscan_table$dbscan_cluster <- as.numeric(NA)

for (i in orders_interested) {
  # Perform DBSCAN clustering on branch length
  dbscan_result <- as.matrix(plot_dbscan_table[clade == i, ]$branch_length) %>%
    dbscan(eps = 0.1, minPts = 5, borderPoints = TRUE) %>%
    .$cluster
  # Store the clustering results in the plot_dbscan_table
  plot_dbscan_table[clade == i, dbscan_cluster := dbscan_result]
}

# Create a scatter plot with jittered points colored by cluster
ggplot(plot_dbscan_table, aes(x = clade, y = branch_length)) +
  geom_jitter(aes(color = as.factor(dbscan_cluster)), alpha = 0.8) +
  scale_color_manual(values = c("#EE6677FF", "#5C80FAFF", "#CFBE31FF")) +
  theme_classic(base_size = 18, base_line_size = 0.5, base_family = "Arial") +
  facet_wrap(~clade, ncol = 3, scales = "free") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(y = "Branch length", x = "Order", color = "Cluster")

# Save the scatter plot as an image file
ggsave("/share/ME-58T/1kp/a_main_graph/figshare_data/FEL_identification/output/yeast_dbscan_cluster.png")

# Assign cluster 0 to the nearest cluster based on the cluster plot
# Identify branch types in Triopsidales, Dipodascales, and Saccharomycodales
branch_type_table <- plot_dbscan_table %>%
  .[clade %in% c("Trigonopsidales", "Dipodascales", "Saccharomycodales")] %>%
  .[clade == "Saccharomycodales" & dbscan_cluster == 0, dbscan_cluster := 1] %>%
  .[, branch_type := ifelse(dbscan_cluster == 2, "FEL", "SEL")]

# Write the branch_type_table to a TSV file
fwrite(branch_type_table, "/share/ME-58T/1kp/a_main_graph/figshare_data/FEL_identification/output/yeast_branch_type_table.tsv", sep = "\t")
fwrite(plot_dbscan_table, "/share/ME-58T/1kp/a_main_graph/figshare_data/FEL_identification/output/yeast_dbscan_table.tsv", sep = "\t")
