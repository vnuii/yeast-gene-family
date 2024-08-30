library(data.table)
library(dplyr)

# Read the branch type table and gene count table
branch_type_dt <- fread("./FEL_identification/output/yeast_branch_type_table.tsv")
yeast_01_og_count <- fread("./orthofinder_results/yeast_0.1_Orthogroups.GeneCount.tsv")

# Function to calculate trimmed row means
trimmed_row_means <- function(dt, columns, trim_threshold = 0.05) {
  apply(dt[, ..columns], 1, function(x) mean(x, trim = trim_threshold))
}

# Function to calculate fold change and perform statistical tests
fc_fun <- function(high_cluster_number, low_cluster_number, clade_need) {
  high_cluster_spnames <- branch_type_dt[clade == clade_need & branch_type == high_cluster_number, ]$id
  low_cluster_spnames <- branch_type_dt[clade == clade_need & branch_type == low_cluster_number, ]$id
  high_cluster_table <- yeast_01_og_count[, ..high_cluster_spnames]
  high_cluster_table$row_mean <- trimmed_row_means(high_cluster_table, high_cluster_spnames)
  low_cluster_table <- yeast_01_og_count[, ..low_cluster_spnames]
  low_cluster_table$row_mean <- trimmed_row_means(low_cluster_table, low_cluster_spnames)
  fc_table <- data.table(Orthogroup = yeast_01_og_count$Orthogroup, fc_value_lowhigh = low_cluster_table$row_mean / high_cluster_table$row_mean)

  pvalues <- numeric(nrow(fc_table))
  for (i in 1:nrow(fc_table)) {
    test_result <- ks.test(unlist(unname(high_cluster_table[i, -"row_mean"])), unlist(unname(low_cluster_table[i, -"row_mean"])))
    pvalues[i] <- test_result$p.value
  }
  fc_table$pvalues <- pvalues
  fc_table$sig_state <- ifelse(fc_table$pvalues > 0.05, "no_change",
    ifelse(fc_table$fc_value_lowhigh >= 1.5 & is.finite(fc_table$fc_value_lowhigh), "expansion",
      ifelse(fc_table$fc_value_lowhigh <= 1 / 1.5 & fc_table$fc_value_lowhigh != 0, "contraction",
        ifelse(is.infinite(fc_table$fc_value_lowhigh), "birth",
          ifelse(fc_table$fc_value_lowhigh == 0, "death", "no_change")
        )
      )
    )
  )
  return(fc_table)
}

# Calculate fold change tables for different clades
Dipodascales_fc_table <- fc_fun("SEL", "FEL", "Dipodascales")
Saccharomycodales_fc_table <- fc_fun("SEL", "FEL", "Saccharomycodales")
Trigonopsidales_fc_table <- fc_fun("SEL", "FEL", "Trigonopsidales")

# Adjust p-values using the Benjamini-Hochberg method
Dipodascales_fc_table$p.adjust <- p.adjust(Dipodascales_fc_table$pvalues, method = "BH")
Saccharomycodales_fc_table$p.adjust <- p.adjust(Saccharomycodales_fc_table$pvalues, method = "BH")
Trigonopsidales_fc_table$p.adjust <- p.adjust(Trigonopsidales_fc_table$pvalues, method = "BH")

# Assign significance states based on adjusted p-values and fold change values
Dipodascales_fc_table$sig_state_adjust <- ifelse(Dipodascales_fc_table$p.adjust > 0.05, "no_change",
  ifelse(Dipodascales_fc_table$fc_value_lowhigh >= 1.5 & is.finite(Dipodascales_fc_table$fc_value_lowhigh), "expansion",
    ifelse(Dipodascales_fc_table$fc_value_lowhigh <= 1 / 1.5 & Dipodascales_fc_table$fc_value_lowhigh != 0, "contraction",
      ifelse(is.infinite(Dipodascales_fc_table$fc_value_lowhigh), "birth",
        ifelse(Dipodascales_fc_table$fc_value_lowhigh == 0, "death", "no_change")
      )
    )
  )
)

Saccharomycodales_fc_table$sig_state_adjust <- ifelse(Saccharomycodales_fc_table$p.adjust > 0.05, "no_change",
  ifelse(Saccharomycodales_fc_table$fc_value_lowhigh >= 1.5 & is.finite(Saccharomycodales_fc_table$fc_value_lowhigh), "expansion",
    ifelse(Saccharomycodales_fc_table$fc_value_lowhigh <= 1 / 1.5 & Saccharomycodales_fc_table$fc_value_lowhigh != 0, "contraction",
      ifelse(is.infinite(Saccharomycodales_fc_table$fc_value_lowhigh), "birth",
        ifelse(Saccharomycodales_fc_table$fc_value_lowhigh == 0, "death", "no_change")
      )
    )
  )
)

Trigonopsidales_fc_table$sig_state_adjust <- ifelse(Trigonopsidales_fc_table$p.adjust > 0.05, "no_change",
  ifelse(Trigonopsidales_fc_table$fc_value_lowhigh >= 1.5 & is.finite(Trigonopsidales_fc_table$fc_value_lowhigh), "expansion",
    ifelse(Trigonopsidales_fc_table$fc_value_lowhigh <= 1 / 1.5 & Trigonopsidales_fc_table$fc_value_lowhigh != 0, "contraction",
      ifelse(is.infinite(Trigonopsidales_fc_table$fc_value_lowhigh), "birth",
        ifelse(Trigonopsidales_fc_table$fc_value_lowhigh == 0, "death", "no_change")
      )
    )
  )
)

# Write fold change tables to files
fwrite(Dipodascales_fc_table, "./fold_change/output/Dipodascales_fc_table.tsv", sep = "\t")
fwrite(Saccharomycodales_fc_table, "./fold_change/output/Saccharomycodales_fc_table.tsv", sep = "\t")
fwrite(Trigonopsidales_fc_table, "./fold_change/output/Trigonopsidales_fc_table.tsv", sep = "\t")

# Prepare enrichment background OGs
Dipodascales_table <- yeast_01_og_count[, branch_type_dt[clade == "Dipodascales"]$id, with = F]
Dipodascales_bg_OGs <- yeast_01_og_count$Orthogroup[rowSums(Dipodascales_table) > 0]

Saccharomycodales_table <- yeast_01_og_count[, branch_type_dt[clade == "Saccharomycodales"]$id, with = F]
Saccharomycodales_bg_OGs <- yeast_01_og_count$Orthogroup[rowSums(Saccharomycodales_table) > 0]

Trigonopsidales_table <- yeast_01_og_count[, branch_type_dt[clade == "Trigonopsidales"]$id, with = F]
Trigonopsidales_bg_OGs <- yeast_01_og_count$Orthogroup[rowSums(Trigonopsidales_table) > 0]

# Write enrichment background OGs to files
fwrite(list(Dipodascales_bg_OGs), "./fold_change/output/enrichment/Dipodascales_bg_OGs.tsv", sep = "\t")
fwrite(list(Saccharomycodales_bg_OGs), "./fold_change/output/enrichment/Saccharomycodales_bg_OGs.tsv", sep = "\t")
fwrite(list(Trigonopsidales_bg_OGs), "./fold_change/output/enrichment/Trigonopsidales_bg_OGs.tsv", sep = "\t")
