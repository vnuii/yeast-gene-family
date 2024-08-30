library(data.table)
library(dplyr)
library(dbscan)
library(foreach)
library(doParallel)

setwd()

yeast_01_og_count <- fread("./orthofinder_results/yeast_0.1_Orthogroups.GeneCount.tsv")
yeast_05_og_count <- fread("./orthofinder_results/yeast_0.5_Orthogroups.GeneCount.tsv")
yeast_info <- fread("./info_yeast.tsv")

# 0.5 coverage threshold PCA
yeast_05_og_count_presenceDT <- yeast_05_og_count %>%
  {
    ifelse(.[, 2:1155] > 0, 1, 0)
  } %>%
  as.data.frame()
rownames(yeast_05_og_count_presenceDT) <- yeast_05_og_count$Orthogroup
yeast_05_og_count_presenceDT_t <- t(yeast_05_og_count_presenceDT)

pca_yeast_05_presence <- prcomp(yeast_05_og_count_presenceDT_t, scale = FALSE)
pca_df_yeast_05_presence <- as.data.frame(pca_yeast_05_presence$x[, 1:10])

fwrite(pca_df_yeast_05_presence, "./PCA/pca_df_yeast_05_presence.tsv", sep = "\t", row.names = TRUE)

clade <- yeast_info$Order[match(rownames(pca_df_yeast_05_presence), yeast_info$assembly_fullID_updated)]
rotation_dt_05_presence <- data.table(OGs = rownames(pca_yeast_05_presence$rotation), pca_yeast_05_presence$rotation)

# there more negative rotation values than positive ones for PC1 and the sum of negative values is higher
rotation_dt_05_presence[PC1 > 0]$PC1 %>%
  sum() %>%
  abs()
rotation_dt_05_presence[PC1 < 0]$PC1 %>%
  sum() %>%
  abs()

# there more negative rotation values than positive ones for PC2 and the sum of negative values is higher
rotation_dt_05_presence[PC2 > 0]$PC2 %>%
  sum() %>%
  abs()
rotation_dt_05_presence[PC2 < 0]$PC2 %>%
  sum() %>%
  abs()

rotation_05_PC1_desc <- rotation_dt_05_presence[, .(OGs, PC1)] %>%
  arrange(desc(PC1))
rotation_05_PC1_asce <- rotation_dt_05_presence[, .(OGs, PC1)] %>%
  arrange(PC1)
rotation_05_PC2_desc <- rotation_dt_05_presence[, .(OGs, PC2)] %>%
  arrange(desc(PC2))
rotation_05_PC2_asce <- rotation_dt_05_presence[, .(OGs, PC2)] %>%
  arrange(PC2)

# representative OGs identification
spearman_table_05_PC1_desc <- data.table(
  PC1_coordinnates = pca_df_yeast_05_presence[, "PC1"],
  clade = clade,
  species = rownames(pca_df_yeast_05_presence)
)
spearman_result_05_PC1_desc <- data.table(
  top_features = NULL,
  spearman_rho = NULL,
  p_value = NULL
)
spearman_table_05_PC1_asce <- data.table(
  PC1_coordinnates = pca_df_yeast_05_presence[, "PC1"],
  clade = clade,
  species = rownames(pca_df_yeast_05_presence)
)
spearman_result_05_PC1_asce <- data.table(
  top_features = NULL,
  spearman_rho = NULL,
  p_value = NULL
)
spearman_table_05_PC2_desc <- data.table(
  PC2_coordinnates = pca_df_yeast_05_presence[, "PC2"],
  clade = clade,
  species = rownames(pca_df_yeast_05_presence)
)
spearman_result_05_PC2_desc <- data.table(
  top_features = NULL,
  spearman_rho = NULL,
  p_value = NULL
)
spearman_table_05_PC2_asce <- data.table(
  PC2_coordinnates = pca_df_yeast_05_presence[, "PC2"],
  clade = clade,
  species = rownames(pca_df_yeast_05_presence)
)
spearman_result_05_PC2_asce <- data.table(
  top_features = NULL,
  spearman_rho = NULL,
  p_value = NULL
)
registerDoParallel(60)

spearman_result_05_PC1_desc <- foreach(i = 1:nrow(yeast_05_og_count), .combine = rbind) %dopar% {
  spearman_table_05_PC1_desc[, paste0("PC1_top", i, "mean") := rowSums(as.data.table(yeast_05_og_count_presenceDT_t)[, rotation_05_PC1_desc$OGs[1:i], with = FALSE]) / i]
  spearman_result <- cor.test(get(paste0("PC1_top", i, "mean"), spearman_table_05_PC1_desc), spearman_table_05_PC1_desc$PC1_coordinnates, method = "spearman")
  tmp <- data.table(
    top_features = i,
    spearman_rho = spearman_result$estimate,
    p_value = spearman_result$p.value
  )
  return(tmp)
} # max spearman correlation is top which.max(spearman_result_05_PC1_desc$spearman_rho) 5
i <- which.max(spearman_result_05_PC1_desc$spearman_rho)
spearman_table_05_PC1_desc[, paste0("PC1_top", i, "mean") := rowSums(as.data.table(yeast_05_og_count_presenceDT_t)[, rotation_05_PC1_desc$OGs[1:i], with = FALSE]) / i]

spearman_result_05_PC1_asce <- foreach(i = 1:nrow(yeast_05_og_count), .combine = rbind) %dopar% {
  spearman_table_05_PC1_asce[, paste0("PC1_top", i, "mean") := rowSums(as.data.table(yeast_05_og_count_presenceDT_t)[, rotation_05_PC1_asce$OGs[1:i], with = FALSE]) / i]
  spearman_result <- cor.test(get(paste0("PC1_top", i, "mean"), spearman_table_05_PC1_asce), spearman_table_05_PC1_asce$PC1_coordinnates, method = "spearman")
  tmp <- data.table(
    top_features = i,
    spearman_rho = spearman_result$estimate,
    p_value = spearman_result$p.value
  )
  return(tmp)
} # max spearman correlation is top which.min(spearman_result_05_PC1_asce$spearman_rho) 610
i <- which.min(spearman_result_05_PC1_asce$spearman_rho)
spearman_table_05_PC1_asce[, paste0("PC1_top", i, "mean") := rowSums(as.data.table(yeast_05_og_count_presenceDT_t)[, rotation_05_PC1_asce$OGs[1:i], with = FALSE]) / i]

## PC2
spearman_result_05_PC2_desc <- foreach(i = 1:nrow(yeast_05_og_count), .combine = rbind) %dopar% {
  spearman_table_05_PC2_desc[, paste0("PC2_top", i, "mean") := rowSums(as.data.table(yeast_05_og_count_presenceDT_t)[, rotation_05_PC2_desc$OGs[1:i], with = FALSE]) / i]
  spearman_result <- cor.test(get(paste0("PC2_top", i, "mean"), spearman_table_05_PC2_desc), spearman_table_05_PC2_desc$PC2_coordinnates, method = "spearman")
  tmp <- data.table(
    top_features = i,
    spearman_rho = spearman_result$estimate,
    p_value = spearman_result$p.value
  )
  return(tmp)
} # max spearman correlation is top which.max(spearman_result_05_PC2_desc$spearman_rho) 8
i <- which.max(spearman_result_05_PC2_desc$spearman_rho)
spearman_table_05_PC2_desc[, paste0("PC2_top", i, "mean") := rowSums(as.data.table(yeast_05_og_count_presenceDT_t)[, rotation_05_PC2_desc$OGs[1:i], with = FALSE]) / i]

spearman_result_05_PC2_asce <- foreach(i = 1:nrow(yeast_05_og_count), .combine = rbind) %dopar% {
  spearman_table_05_PC2_asce[, paste0("PC2_top", i, "mean") := rowSums(as.data.table(yeast_05_og_count_presenceDT_t)[, rotation_05_PC2_asce$OGs[1:i], with = FALSE]) / i]
  spearman_result <- cor.test(get(paste0("PC2_top", i, "mean"), spearman_table_05_PC2_asce), spearman_table_05_PC2_asce$PC2_coordinnates, method = "spearman")
  tmp <- data.table(
    top_features = i,
    spearman_rho = spearman_result$estimate,
    p_value = spearman_result$p.value
  )
  return(tmp)
} # max spearman correlation is top which.min(spearman_result_05_PC2_asce$spearman_rho) 27
i <- which.min(spearman_result_05_PC2_asce$spearman_rho)
spearman_table_05_PC2_asce[, paste0("PC2_top", i, "mean") := rowSums(as.data.table(yeast_05_og_count_presenceDT_t)[, rotation_05_PC2_asce$OGs[1:i], with = FALSE]) / i]

stopImplicitCluster()
fwrite(spearman_table_05_PC1_desc, "./PCA/threshold_05/spearman_table_05_PC1_desc.tsv", sep = "\t")
fwrite(spearman_result_05_PC1_desc, "./PCA/threshold_05/spearman_result_05_PC1_desc.tsv", sep = "\t")
fwrite(spearman_table_05_PC1_asce, "./PCA/threshold_05/spearman_table_05_PC1_asce.tsv", sep = "\t")
fwrite(spearman_result_05_PC1_asce, "./PCA/threshold_05/spearman_result_05_PC1_asce.tsv", sep = "\t")
fwrite(spearman_table_05_PC2_desc, "./PCA/threshold_05/spearman_table_05_PC2_desc.tsv", sep = "\t")
fwrite(spearman_result_05_PC2_desc, "./PCA/threshold_05/spearman_result_05_PC2_desc.tsv", sep = "\t")
fwrite(spearman_table_05_PC2_asce, "./PCA/threshold_05/spearman_table_05_PC2_asce.tsv", sep = "\t")
fwrite(spearman_result_05_PC2_asce, "./PCA/threshold_05/spearman_result_05_PC2_asce.tsv", sep = "\t")


# 0.1 coverage threshold PCA
yeast_01_og_count_presenceDT <- yeast_01_og_count %>%
  {
    ifelse(.[, 2:1155] > 0, 1, 0)
  } %>%
  as.data.frame()
rownames(yeast_01_og_count_presenceDT) <- yeast_01_og_count$Orthogroup
yeast_01_og_count_presenceDT_t <- t(yeast_01_og_count_presenceDT)

pca_yeast_01_presence <- prcomp(yeast_01_og_count_presenceDT_t, scale = FALSE)
pca_df_yeast_01_presence <- as.data.frame(pca_yeast_01_presence$x[, 1:10])

fwrite(pca_df_yeast_01_presence, "./PCA/pca_df_yeast_01_presence.tsv", sep = "\t", row.names = TRUE)

clade <- yeast_info$Order[match(rownames(pca_df_yeast_01_presence), yeast_info$assembly_fullID_updated)]
rotation_dt_01_presence <- data.table(OGs = rownames(pca_yeast_01_presence$rotation), pca_yeast_01_presence$rotation)

# there more positive rotation values than negative ones for PC1 and the sum of positive values is higher
rotation_dt_01_presence[PC1 > 0]$PC1 %>%
  sum() %>%
  abs()
rotation_dt_01_presence[PC1 < 0]$PC1 %>%
  sum() %>%
  abs()

rotation_01_PC1_desc <- rotation_dt_01_presence[, .(OGs, PC1)] %>%
  arrange(desc(PC1))
rotation_01_PC1_asce <- rotation_dt_01_presence[, .(OGs, PC1)] %>%
  arrange(PC1)

# representative OGs identification
spearman_table_01_PC1_desc <- data.table(
  PC1_coordinnates = pca_df_yeast_01_presence[, "PC1"],
  clade = clade,
  species = rownames(pca_df_yeast_01_presence)
)
spearman_result_01_PC1_desc <- data.table(
  top_features = NULL,
  spearman_rho = NULL,
  p_value = NULL
)
spearman_table_01_PC1_asce <- data.table(
  PC1_coordinnates = pca_df_yeast_01_presence[, "PC1"],
  clade = clade,
  species = rownames(pca_df_yeast_01_presence)
)
spearman_result_01_PC1_asce <- data.table(
  top_features = NULL,
  spearman_rho = NULL,
  p_value = NULL
)

registerDoParallel(60)

spearman_result_01_PC1_desc <- foreach(i = 1:nrow(yeast_01_og_count), .combine = rbind) %dopar% {
  spearman_table_01_PC1_desc[, paste0("PC1_top", i, "mean") := rowSums(as.data.table(yeast_01_og_count_presenceDT_t)[, rotation_01_PC1_desc$OGs[1:i], with = FALSE]) / i]
  spearman_result <- cor.test(get(paste0("PC1_top", i, "mean"), spearman_table_01_PC1_desc), spearman_table_01_PC1_desc$PC1_coordinnates, method = "spearman")
  tmp <- data.table(
    top_features = i,
    spearman_rho = spearman_result$estimate,
    p_value = spearman_result$p.value
  )
  return(tmp)
} # max spearman correlation is top which.max(spearman_result_01_PC1_desc$spearman_rho) 421
i <- which.max(spearman_result_01_PC1_desc$spearman_rho)
spearman_table_01_PC1_desc[, paste0("PC1_top", i, "mean") := rowSums(as.data.table(yeast_01_og_count_presenceDT_t)[, rotation_01_PC1_desc$OGs[1:i], with = FALSE]) / i]

spearman_result_01_PC1_asce <- foreach(i = 1:nrow(yeast_01_og_count), .combine = rbind) %dopar% {
  spearman_table_01_PC1_asce[, paste0("PC1_top", i, "mean") := rowSums(as.data.table(yeast_01_og_count_presenceDT_t)[, rotation_01_PC1_asce$OGs[1:i], with = FALSE]) / i]
  spearman_result <- cor.test(get(paste0("PC1_top", i, "mean"), spearman_table_01_PC1_asce), spearman_table_01_PC1_asce$PC1_coordinnates, method = "spearman")
  tmp <- data.table(
    top_features = i,
    spearman_rho = spearman_result$estimate,
    p_value = spearman_result$p.value
  )
  return(tmp)
} # max spearman correlation is top which.min(spearman_result_01_PC1_asce$spearman_rho) 3
i <- which.min(spearman_result_01_PC1_asce$spearman_rho)
spearman_table_01_PC1_asce[, paste0("PC1_top", i, "mean") := rowSums(as.data.table(yeast_01_og_count_presenceDT_t)[, rotation_01_PC1_asce$OGs[1:i], with = FALSE]) / i]

stopImplicitCluster()
fwrite(spearman_table_01_PC1_desc, "./PCA/threshold_01/spearman_table_01_PC1_desc.tsv", sep = "\t")
fwrite(spearman_result_01_PC1_desc, "./PCA/threshold_01/spearman_result_01_PC1_desc.tsv", sep = "\t")
fwrite(spearman_table_01_PC1_asce, "./PCA/threshold_01/spearman_table_01_PC1_asce.tsv", sep = "\t")
fwrite(spearman_result_01_PC1_asce, "./PCA/threshold_01/spearman_result_01_PC1_asce.tsv", sep = "\t")
