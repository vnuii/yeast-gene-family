library(data.table)
library(dplyr)
library(KEGGREST)

# Read the eMapper annotations file
emapper <- fread("./y1000.final.emapper.annotations.backup", select = c("query", "Description", "Preferred_name"), key = "query") # ignore the warning

# Read the Orthogroups gene table
OG_gene_table <- fread("./orthofinder_results/yeast_Orthogroups_1154_geneIDs.tsv")

# Read the Dipodascales fold change table
Dipo_foldchange_dt <- fread("./fold_change/output/Dipodascales_fc_table.tsv")

# Get the Orthogroups with contraction or death as the significant state
Dipo_contractionlost_OGs <- Dipo_foldchange_dt[sig_state_adjust %in% c("contraction", "death")]$Orthogroup

# Read the top 60 Orthogroups from the PCA results
PC1_05_top610OGs <- fread("./PCA/threshold_05/PC1_asce_top610OGs.txt", header = FALSE) %>%
  pull(V1)

# Select relevant columns from the Orthogroups gene table for Saccharomyces cerevisiae and Candida albicans
OG_gene_table_sce_cal <- OG_gene_table[, c("Orthogroup", "saccharomyces_cerevisiae.final", "candida_albicans.final"), with = FALSE]

# Filter the Orthogroups gene table for the Orthogroups with contraction or death as the significant state
OG_gene_table_sce_cal_foldchange <- OG_gene_table_sce_cal[Orthogroup %in% Dipo_contractionlost_OGs]

# Function to get gene names based on Orthogroup and gene IDs
get_gene_name <- function(x) {
  OG <- x[[1]]
  gene_sce <- x[[2]] %>%
    str_split(",") %>%
    unlist
  gene_cal <- x[[3]] %>%
    str_split(",") %>%
    unlist
  gene_name_sce <- emapper[query %in% gene_sce]$Preferred_name %>% unique %>% .[. != "-"]
  gene_name_cal <- emapper[query %in% gene_cal]$Preferred_name %>% unique %>% .[. != "-"]
  if (length(gene_name_sce) > 0) {
    gene <- gene_sce
    name <- gene_name_sce
  } else if (length(gene_name_cal) > 0) {
    gene <- gene_cal
    name <- gene_name_cal
  } else {
    gene <- ""
    name <- ""
  }
  out <- data.table(Orthogroup = OG, gene_name = name)
  return(out)
}

# Apply the get_gene_name function to the Orthogroups gene table for fold change analysis
foldchange_OGs_GeneName_table <- apply(OG_gene_table_sce_cal_foldchange, 1, get_gene_name) %>%
  rbindlist() %>%
  .[gene_name != ""]

# Convert gene names to uppercase
foldchange_OGs_GeneName_table$gene_name <- toupper(foldchange_OGs_GeneName_table$gene_name)

# Replace "GPD" with "TDH3" in the gene names
foldchange_OGs_GeneName_table[gene_name == "GPD"]$gene_name <- "TDH3"

# Apply the get_gene_name function to the Orthogroups gene table for PCA analysis
OG_gene_table_sce_cal_PCA <- OG_gene_table_sce_cal[Orthogroup %in% PC1_05_top610OGs]

# Get gene names for the PCA analysis
PCA_OGs_GeneName_table <- apply(OG_gene_table_sce_cal_PCA, 1, get_gene_name) %>%
  rbindlist() %>%
  .[gene_name != ""]

# Convert gene names to uppercase
PCA_OGs_GeneName_table$gene_name <- toupper(PCA_OGs_GeneName_table$gene_name)

# Function to get KEGG IDs for gene names
get_kegg_id <- function(gene_name) {
  # Use KEGG API to query for KEGG IDs
  result <- keggFind("sce", gene_name)
  names_need <- names(result) %>%
    gsub("sce:", "", .)
  tmp_table <- data.table(kegg_id = names_need)
  tmp_table$gene_name <- gene_name
  return(tmp_table)
}

# Get KEGG IDs for the gene names in the fold change analysis
KEGG_id_Dipo_foldchange <- lapply(foldchange_OGs_GeneName_table$gene_name, get_kegg_id) %>%
  rbindlist()
# Match the Orthogroup column in KEGG_id_Dipo_foldchange with the Orthogroup column in foldchange_OGs_GeneName_table
KEGG_id_Dipo_foldchange$Orthogroup <- foldchange_OGs_GeneName_table$Orthogroup[match(KEGG_id_Dipo_foldchange$gene_name, foldchange_OGs_GeneName_table$gene_name)]

# Write the fold change analysis results to a TSV file
fwrite(KEGG_id_Dipo_foldchange[, c("Orthogroup", "gene_name", "kegg_id")], "./gene_digging/output/Dipodascales_foldchange_OGs_GeneName_table.tsv", sep = "\t")

# Get KEGG IDs for the gene names in the PCA analysis
KEGG_id_Dipo_PCA <- lapply(PCA_OGs_GeneName_table$gene_name, get_kegg_id) %>%
  rbindlist()

# Match the Orthogroup column in KEGG_id_Dipo_PCA with the Orthogroup column in PCA_OGs_GeneName_table
KEGG_id_Dipo_PCA$Orthogroup <- PCA_OGs_GeneName_table$Orthogroup[match(KEGG_id_Dipo_PCA$gene_name, PCA_OGs_GeneName_table$gene_name)]

# Write the PCA analysis results to a TSV file
fwrite(KEGG_id_Dipo_PCA[, c("Orthogroup", "gene_name", "kegg_id")], "./gene_digging/output/Dipodascales_PC1_05_asce_top610OGs_GeneName_table.tsv", sep = "\t")
