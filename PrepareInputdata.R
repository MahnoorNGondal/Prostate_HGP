
## 2 - data for the LPM - RNA-seq
```{r}
## RNAseq
d.expression.tbl_tumor_rpkm_sel <- read.csv("/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_27_11_12_25/d.expression.tbl_tumor_rpkm_sel.csv")
rownames(d.expression.tbl_tumor_rpkm_sel) <- d.expression.tbl_tumor_rpkm_sel$X
d.expression.tbl_tumor_rpkm_sel$X <- NULL
rna <- as.data.frame(t(d.expression.tbl_tumor_rpkm_sel))
rownames(rna) <- gsub("\\.", "-", rownames(rna))
# write.csv(rna, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/RNAseq_data_for_LPM.csv", row.names = TRUE)

############################################
## Setup
############################################

library(dplyr)
library(limma)
library(msigdbr)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(purrr)

############################################
## 1. Align metadata and RNA-seq by case_id
############################################

# meta: metadata data.frame with columns case_id, Purity, ETS
# rna:  raw RPKM matrix/data.frame (samples x genes), rownames = case_id

# Keep only samples present in RNA and reorder RNA to match meta
meta <- meta %>%
  filter(case_id %in% rownames(rna))
rna  <- rna[meta$case_id, , drop = FALSE]

# Purity covariate
purity <- as.numeric(meta$Purity)
names(purity) <- meta$case_id

############################################
## 2. Log2 transform and filter low-expression genes
############################################

# Log2(RPKM + 1)
log_rna <- log2(as.matrix(rna) + 1)

# Filter genes expressed (RPKM > 0) in at least 5% of samples
expr_mat_raw <- 2^log_rna - 1           # back to RPKM scale
genes_keep   <- colSums(expr_mat_raw > 0) > (0.05 * nrow(expr_mat_raw))
log_rna      <- log_rna[, genes_keep, drop = FALSE]

############################################
## 3. Purity correction with limma
############################################

# BEFORE correction: PCA on top 5000 variable genes
# (use the same top genes selection as after correction for fair comparison)

# Variance on uncorrected log-rna
gene_var_pre <- apply(log_rna, 2, var, na.rm = TRUE)
top_genes    <- names(sort(gene_var_pre, decreasing = TRUE))[1:5000]

log_rna_top_pre <- log_rna[, top_genes, drop = FALSE]

# Z-score and PCA (before purity correction)
z_rna_pre <- scale(log_rna_top_pre, center = TRUE, scale = TRUE)
pca_pre   <- prcomp(z_rna_pre, center = FALSE, scale. = FALSE)
pc_pre    <- as.data.frame(pca_pre$x[, 1:2])   # PC1 & PC2 for plotting
pc_pre$case_id <- rownames(pca_pre$x)
pc_pre <- pc_pre %>%
  left_join(meta[, c("case_id", "ETS", "Purity")], by = "case_id")

# Purity correction: limma::removeBatchEffect expects genes x samples
log_rna_pc <- removeBatchEffect(t(log_rna), covariates = purity)
log_rna_pc <- t(log_rna_pc)   # back to samples x genes

############################################
## 4. Select top 5000 variable genes (after correction) and PCA
############################################

gene_var_post <- apply(log_rna_pc, 2, var, na.rm = TRUE)
top_genes     <- names(sort(gene_var_post, decreasing = TRUE))[1:5000]

log_rna_pc_top <- log_rna_pc[, top_genes, drop = FALSE]

# Z-score and PCA (after purity correction)
z_rna_post <- scale(log_rna_pc_top, center = TRUE, scale = TRUE)
pca_post   <- prcomp(z_rna_post, center = FALSE, scale. = FALSE)
pc_post    <- as.data.frame(pca_post$x[, 1:2])   # PC1 & PC2
pc_post$case_id <- rownames(pca_post$x)
pc_post <- pc_post %>%
  left_join(meta[, c("case_id", "ETS", "Purity")], by = "case_id")

############################################
## 5. PCA plots colored by ETS (before & after)
############################################

p_pca_pre <- ggplot(pc_pre, aes(x = PC1, y = PC2, color = ETS)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw() +
  labs(title = "PCA before purity correction",
       x = "PC1", y = "PC2")

p_pca_post <- ggplot(pc_post, aes(x = PC1, y = PC2, color = ETS)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw() +
  labs(title = "PCA after purity correction",
       x = "PC1", y = "PC2")

# Optionally arrange side by side
ggarrange(p_pca_pre, p_pca_post, ncol = 2)

############################################
## 6. K-means clustering and purity vs clusters
############################################

set.seed(123)  # for reproducibility
k <- 3         # choose number of clusters (adjust as needed)

# Before correction: clustering on first 10 PCs
pc_for_kmeans_pre <- pca_pre$x[, 1:10]
km_pre <- kmeans(pc_for_kmeans_pre, centers = k, nstart = 50)

meta$cluster_pre <- factor(km_pre$cluster)

# After correction: clustering on first 10 PCs
pc_for_kmeans_post <- pca_post$x[, 1:10]
km_post <- kmeans(pc_for_kmeans_post, centers = k, nstart = 50)

meta$cluster_post <- factor(km_post$cluster)

# Boxplot: purity vs k-means clusters (before)
p_box_pre <- ggboxplot(
  meta, x = "cluster_pre", y = "Purity",
  color = "cluster_pre", add = "jitter"
) +
  theme_bw() +
  labs(title = "Purity vs k-means clusters (before correction)",
       x = "k-means cluster", y = "Purity")

# Boxplot: purity vs k-means clusters (after)
p_box_post <- ggboxplot(
  meta, x = "cluster_post", y = "Purity",
  color = "cluster_post", add = "jitter"
) +
  theme_bw() +
  labs(title = "Purity vs k-means clusters (after correction)",
       x = "k-means cluster", y = "Purity")

ggarrange(p_box_pre, p_box_post, ncol = 2)

############################################
## 7. Full PCA (first 50 PCs) after correction
############################################

rna_pc_scores <- pca_post$x[, 1:50, drop = FALSE]
write.csv(rna_pc_scores, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/LPM_rna_pc_scores.csv", row.names = TRUE)
# ############################################
# ## 8. MSigDB C2 gene sets and ssGSEA
# ############################################

# # Get C2 gene sets (human)
# m_df <- msigdbr(species = "Homo sapiens", category = "C2")

# # List: pathway_name -> vector of gene symbols
# c2_list <- m_df %>%
#   split(.$gs_name) %>%
#   lapply(function(z) unique(z$gene_symbol))

# # Expression for GSVA: genes x samples
# expr_gsva <- t(log_rna_pc_top)  # same data used for PCA, corrected

# # Compute ssGSEA/GSVA scores
# c2_scores <- gsva(
#   expr = expr_gsva,
#   gset.idx.list = c2_list,
#   method = "ssgsea",     # or "gsva"
#   kcdf = "Gaussian",
#   parallel.sz = 1
# )
# # c2_scores: pathways x samples

# ############################################
# ## 9. Correlate PCs with pathway scores
# ############################################

# # Align samples between PCs and pathway scores
# common_samples <- intersect(rownames(rna_pc_scores), colnames(c2_scores))

# pc_mat   <- as.matrix(rna_pc_scores[common_samples, , drop = FALSE])       # samples x PCs
# path_mat <- t(c2_scores[, common_samples, drop = FALSE])               # samples x pathways

# # Correlation matrix: PCs x pathways
# cor_mat <- cor(pc_mat, path_mat,
#                use = "pairwise.complete.obs",
#                method = "pearson")

# # Tidy correlation data
# cor_df <- as.data.frame(as.table(cor_mat))
# colnames(cor_df) <- c("PC", "Pathway", "Correlation")

# # Optional: compute p-values for each PC–pathway pair
# get_cor_p <- function(pc_name, pw_name) {
#   x <- pc_mat[, pc_name]
#   y <- path_mat[, pw_name]
#   if (all(is.na(x)) || all(is.na(y))) return(NA_real_)
#   suppressWarnings(cor.test(x, y)$p.value)
# }

# pc_names <- rownames(cor_mat)
# pw_names <- colnames(cor_mat)

# p_mat <- map_dbl(expand.grid(pc_names, pw_names, KEEP.OUT.ATTRS = FALSE),
#                  ~ get_cor_p(.x[1], .x[2])) %>%
#   matrix(nrow = length(pc_names), byrow = TRUE,
#          dimnames = list(pc_names, pw_names))

# p_df <- as.data.frame(as.table(p_mat))
# colnames(p_df) <- c("PC", "Pathway", "Pvalue")

# # Final annotated correlation table
# cor_annot <- cor_df %>%
#   left_join(p_df, by = c("PC", "Pathway")) %>%
#   arrange(PC, desc(abs(Correlation)))

# # Example: top pathways for PC1
# pc1_top <- cor_annot %>%
#   filter(PC == "PC1") %>%
#   arrange(desc(abs(Correlation))) %>%
#   head(30)


```

## 3 - data for the LPM - Proteomics
```{r}
proteome <- read.csv("/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_30_12_17_25/prostate_CPTAC_prot_TMT_abundance_gene_MD.tsv", sep = "\t")
proteome <- proteome[,-c(1)]
rownames(proteome) <- proteome$Index
proteome$Index <- NULL


# 2. Transpose expression matrix: samples as rows, genes as columns
expr_for_pca <- as.data.frame(t(as.matrix(proteome)))

expr_for_pca$remove <- ifelse(grepl("NCI7", rownames(expr_for_pca)), "NCI7", "keep")
expr_for_pca <- dplyr::filter(expr_for_pca, remove != "NCI7")
expr_for_pca$remove <- NULL
expr_for_pca$sample <- ifelse(grepl(".T", rownames(expr_for_pca)), "tumor", "normal")
protein <- dplyr::filter(expr_for_pca, sample == "tumor")
protein$sample <- NULL
rownames(protein) <- gsub(".T", "", rownames(protein))
rownames(protein) <- gsub("\\.", "-", rownames(protein))


library(org.Hs.eg.db)
library(AnnotationDbi)

# 1) Current column names are Ensembl IDs
ens_ids <- colnames(protein)

# 2) Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys       = ens_ids,
  column     = "SYMBOL",      # what you want back
  keytype    = "ENSEMBL",     # what your IDs are
  multiVals  = "first"        # how to handle 1:many
)

# 3) Replace NAs with original ID (optional but recommended)
gene_symbols[is.na(gene_symbols)] <- ens_ids[is.na(gene_symbols)]

# 4) Assign new names
colnames(protein) <- gene_symbols
# write.csv(protein, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/Proteomics_data_for_LPM.csv", row.names = TRUE)

############################################
## Protein Analysis: 4-Stage Pipeline
############################################

library(dplyr)
library(imputeLCMD)
library(limma)
library(ggplot2)
library(ggpubr)

############################################
## 1. Setup & SAFE PCA FUNCTION
############################################

meta_clean <- meta %>% filter(case_id %in% rownames(protein))
protein <- protein[meta_clean$case_id, , drop = FALSE]
purity <- as.numeric(meta_clean$Purity)
names(purity) <- meta_clean$case_id

cat("Dataset:", nrow(protein), "samples x", ncol(protein), "proteins\n")

clean_for_pca <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if(ncol(mat) == 0 || nrow(mat) == 0) stop("Empty matrix")
  
  # Remove constants
  vars <- apply(mat, 2, var, na.rm = TRUE)
  keep_cols <- vars > 0 & is.finite(vars)
  mat <- mat[, keep_cols, drop = FALSE]
  
  if(ncol(mat) == 0) stop("No variable columns")
  
  # NA → median
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  # Final cleanup
  mat[!is.finite(mat)] <- NA
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  return(mat)
}

n_top <- min(5000, ncol(protein))

do_pca_analysis <- function(prot_mat, stage_name, meta_data) {
  prot_clean <- clean_for_pca(prot_mat)
  var_order <- apply(prot_clean, 2, var, na.rm = TRUE)
  top_prot <- names(sort(var_order, decreasing = TRUE))[1:n_top]
  
  prot_top <- prot_clean[, top_prot, drop = FALSE]
  z_prot <- scale(prot_top, center = TRUE, scale = TRUE)
  pca_obj <- prcomp(z_prot, center = FALSE, scale. = FALSE)
  
  pc_df <- data.frame(pca_obj$x[, 1:2], case_id = rownames(pca_obj$x)) %>%
    left_join(meta_data[, c("case_id", "ETS", "Purity")], by = "case_id")
  
  km <- kmeans(pca_obj$x[, 1:min(10, ncol(pca_obj$x))], centers = 3, nstart = 50)
  meta_data[[paste0("cluster_", stage_name)]] <- factor(km$cluster)
  
  list(pca = pca_obj, pc_df = pc_df, meta = meta_data)
}

set.seed(123)

############################################
## 2. STAGE 1: RAW (Before imputation, Before purity correction)
############################################

cat("\n=== STAGE 1: Raw data ===\n")
stage1 <- do_pca_analysis(protein, "raw", meta_clean)
pc1_df <- stage1$pc_df
meta_clean <- stage1$meta

############################################
## 3. STAGE 2: RAW (After purity correction, Before imputation)  
############################################

cat("\n=== STAGE 2: Raw + Purity correction ===\n")
prot_raw_t <- t(clean_for_pca(protein))
prot_raw_pc_t <- removeBatchEffect(prot_raw_t, covariates = purity)
protein_raw_pc <- t(prot_raw_pc_t)

stage2 <- do_pca_analysis(protein_raw_pc, "raw_pc", meta_clean)
pc2_df <- stage2$pc_df
meta_clean <- stage2$meta

############################################
## 4. STAGE 3: IMPUTED (After imputation, Before purity correction)
############################################

cat("\n=== STAGE 3: Imputed data ===\n")
qrilc_result <- impute.QRILC(as.matrix(protein), tune.sigma = 1)
protein_imputed <- qrilc_result[[1]]

stage3 <- do_pca_analysis(protein_imputed, "imp", meta_clean)
pc3_df <- stage3$pc_df
meta_clean <- stage3$meta

############################################
## 5. STAGE 4: IMPUTED + PURITY CORRECTED (FINAL)
############################################

cat("\n=== STAGE 4: Imputed + Purity corrected ===\n")
prot_imp_t <- t(clean_for_pca(protein_imputed))
prot_imp_pc_t <- removeBatchEffect(prot_imp_t, covariates = purity)
protein_final <- t(prot_imp_pc_t)

stage4 <- do_pca_analysis(protein_final, "final", meta_clean)
pc4_df <- stage4$pc_df
meta_clean <- stage4$meta

cat("Pipeline complete!\n")

############################################
## 6. 4-STAGE PCA PLOTS
############################################

p1 <- ggplot(pc1_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "1. RAW\n(Before imputation & purity cor)", x = "PC1", y = "PC2")

p2 <- ggplot(pc2_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "2. RAW + Purity corrected\n(Before imputation)", x = "PC1", y = "PC2")

p3 <- ggplot(pc3_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "3. IMPUTED\n(Before purity correction)", x = "PC1", y = "PC2")

p4 <- ggplot(pc4_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "4. IMPUTED + Purity corrected\n(FINAL)", x = "PC1", y = "PC2")

# 2x2 PCA grid
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 7. 4-STAGE PURITY BOXPLOTS
############################################

p_box1 <- ggboxplot(meta_clean, x = "cluster_raw", y = "Purity", 
                   color = "cluster_raw", add = "jitter") + theme_bw() +
  labs(title = "1. RAW\n(Before imputation & purity cor)")

p_box2 <- ggboxplot(meta_clean, x = "cluster_raw_pc", y = "Purity", 
                   color = "cluster_raw_pc", add = "jitter") + theme_bw() +
  labs(title = "2. RAW + Purity corrected\n(Before imputation)")

p_box3 <- ggboxplot(meta_clean, x = "cluster_imp", y = "Purity", 
                   color = "cluster_imp", add = "jitter") + theme_bw() +
  labs(title = "3. IMPUTED\n(Before purity correction)")

p_box4 <- ggboxplot(meta_clean, x = "cluster_final", y = "Purity", 
                   color = "cluster_final", add = "jitter") + theme_bw() +
  labs(title = "4. IMPUTED + Purity corrected\n(FINAL)")

# 2x2 Boxplot grid
ggarrange(p_box1, p_box2, p_box3, p_box4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 8. SUMMARY
############################################

cat("\n=== 4-STAGE SUMMARY ===\n")
cat("PC1 variances (%):")
cat(sprintf(" RAW: %.1f, RAW+PC: %.1f, IMP: %.1f, FINAL: %.1f\n",
    summary(stage1$pca)$importance[2,1]*100,
    summary(stage2$pca)$importance[2,1]*100, 
    summary(stage3$pca)$importance[2,1]*100,
    summary(stage4$pca)$importance[2,1]*100))

# Final objects
protein_pc_scores <- stage4$pca$x[, 1:min(50, ncol(stage4$pca$x))]
meta_clean$cluster_final <- meta_clean$cluster_final

print("Key outputs:")
print("- protein_pc_scores: Final PCA scores")
print("- protein_final: Final corrected matrix") 
print("- meta_clean: All 4-stage clusters")
print("- pc1_df, pc2_df, pc3_df, pc4_df: PCA dataframes")

# Add this to your pipeline to get ALL metrics:
evaluate_stage <- function(pca_obj, meta, cluster_col, purity_col = "Purity") {
  pc1_var <- summary(pca_obj)$importance[2,1]*100
  pc2_var <- summary(pca_obj)$importance[2,2]*100
  
  # ETS silhouette (biological separation)
  library(cluster)
  sil_ets <- silhouette(as.numeric(as.factor(meta$ETS)), dist(pca_obj$x[,1:2]))
  ets_sil <- mean(sil_ets[,3])
  
  # Purity ANOVA p-value (technical removal)
  purity_p <- summary(aov(meta[[purity_col]] ~ meta[[cluster_col]]))[[1]]$`Pr(>F)`[1]
  
  data.frame(
    PC1_var = pc1_var, PC2_var = pc2_var, ETS_silhouette = ets_sil, 
    Purity_pval = purity_p
  )
}

# Compare all 4 stages
metrics <- rbind(
  evaluate_stage(stage1$pca, meta_clean, "cluster_raw"),
  evaluate_stage(stage2$pca, meta_clean, "cluster_raw_pc"),
  evaluate_stage(stage3$pca, meta_clean, "cluster_imp"),
  evaluate_stage(stage4$pca, meta_clean, "cluster_final")
)
rownames(metrics) <- c("RAW", "RAW+PC", "IMP", "FINAL")
print(metrics)

############################################
## EXTRACT BEST STAGE (RAW+PC) - TOP 50 PCs
############################################

# Stage 2 (RAW+PC) was optimal per metrics
protein_pc_top50 <- stage2$pca$x[, 1:50, drop = FALSE]
write.csv(protein_pc_top50, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/LPM_protein_pc_top50.csv", row.names = TRUE)
# # Clean column names + add to metadata
# colnames(raw_pc_top50) <- paste0("PC_", 1:ncol(raw_pc_top50))
# rownames(raw_pc_top50) <- meta_clean$case_id

# # Merge with meta_clean
# meta_final <- cbind(meta_clean, raw_pc_top50)

# cat("\n🏆 RAW+PC (Stage 2) - BEST STAGE EXTRACTED\n")
# cat("Top 50 PCs:", nrow(raw_pc_top50), "samples ×", ncol(raw_pc_top50), "PCs\n")
# cat("PC1 variance:", round(summary(stage2$pca)$importance[2,1]*100,1), "%\n")
# cat("PC1+PC2 variance:", round(sum(summary(stage2$pca)$importance[2,1:2]*100),1), "%\n")
# cat("\nKey objects ready:\n")
# cat("- raw_pc_top50: standalone PC matrix\n") 
# cat("- meta_final: metadata + PCs (recommended)\n")
# cat("- protein_raw_pc: corrected protein matrix\n")

# Optional: Save results
# write.csv(raw_pc_top50, "raw_pc_top50.csv", row.names = TRUE)
# write.csv(meta_final, "meta_with_pcs.csv", row.names = FALSE)

```

## 4 - data for the LPM - PTMs:Phosphotyrosine 1
```{r}
pY <- read.csv("/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/pY_abundance_gene_MD.tsv", sep = "\t")
pY <- pY[,-c(2:5)]
rownames(pY) <- pY$Index
pY$Index <- NULL


# 2. Transpose expression matrix: samples as rows, genes as columns
pY_data <- as.data.frame(t(as.matrix(pY)))

pY_data$remove <- ifelse(grepl("NCI7", rownames(pY_data)), "NCI7", "keep")
pY_data <- dplyr::filter(pY_data, remove != "NCI7")
pY_data$remove <- NULL
pY_data$sample <- ifelse(grepl(".T", rownames(pY_data)), "tumor", "normal")
pY_data <- dplyr::filter(pY_data, sample == "tumor")
pY_data$sample <- NULL
rownames(pY_data) <- gsub(".T", "", rownames(pY_data))
rownames(pY_data) <- gsub("\\.", "-", rownames(pY_data))


library(org.Hs.eg.db)
library(AnnotationDbi)

# 1) Current column names are Ensembl IDs
ens_ids <- colnames(pY_data)

# 2) Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys       = ens_ids,
  column     = "SYMBOL",      # what you want back
  keytype    = "ENSEMBL",     # what your IDs are
  multiVals  = "first"        # how to handle 1:many
)

# 3) Replace NAs with original ID (optional but recommended)
gene_symbols[is.na(gene_symbols)] <- ens_ids[is.na(gene_symbols)]

# 4) Assign new names
colnames(pY_data) <- gene_symbols
# write.csv(pY_data, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/Phosphotyrosine_data_for_LPM.csv", row.names = TRUE)
############################################
## pY (Phosphotyrosine) Analysis: 4-Stage Pipeline
############################################

library(dplyr)
library(imputeLCMD)
library(limma)
library(ggplot2)
library(ggpubr)
library(cluster)

############################################
## 1. Setup & SAFE PCA FUNCTION (pY-specific)
############################################

meta_clean_py <- meta %>% filter(case_id %in% rownames(pY_data))
pY_data_aligned <- pY_data[meta_clean_py$case_id, , drop = FALSE]
purity_py <- as.numeric(meta_clean_py$Purity)
names(purity_py) <- meta_clean_py$case_id

cat("pY Dataset:", nrow(pY_data_aligned), "samples x", ncol(pY_data_aligned), "pY sites\n")

clean_for_pca <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if(ncol(mat) == 0 || nrow(mat) == 0) stop("Empty matrix")
  
  # Remove constants
  vars <- apply(mat, 2, var, na.rm = TRUE)
  keep_cols <- vars > 0 & is.finite(vars)
  mat <- mat[, keep_cols, drop = FALSE]
  
  if(ncol(mat) == 0) stop("No variable columns")
  
  # NA → median
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  # Final cleanup
  mat[!is.finite(mat)] <- NA
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  return(mat)
}

n_top_py <- min(5000, ncol(pY_data_aligned))

do_pca_analysis <- function(prot_mat, stage_name, meta_data) {
  prot_clean <- clean_for_pca(prot_mat)
  var_order <- apply(prot_clean, 2, var, na.rm = TRUE)
  top_prot <- names(sort(var_order, decreasing = TRUE))[1:n_top_py]
  
  prot_top <- prot_clean[, top_prot, drop = FALSE]
  z_prot <- scale(prot_top, center = TRUE, scale = TRUE)
  pca_obj <- prcomp(z_prot, center = FALSE, scale. = FALSE)
  
  pc_df <- data.frame(pca_obj$x[, 1:2], case_id = rownames(pca_obj$x)) %>%
    left_join(meta_data[, c("case_id", "ETS", "Purity")], by = "case_id")
  
  km <- kmeans(pca_obj$x[, 1:min(10, ncol(pca_obj$x))], centers = 4, nstart = 50)
  meta_data[[paste0("pY_cluster_", stage_name)]] <- factor(km$cluster)
  
  list(pca = pca_obj, pc_df = pc_df, meta = meta_data)
}

set.seed(123)

############################################
## 2. STAGE 1: RAW pY (Before imputation, Before purity correction)
############################################

cat("\n=== pY STAGE 1: Raw data ===\n")
py_stage1 <- do_pca_analysis(pY_data_aligned, "raw", meta_clean_py)
py_pc1_df <- py_stage1$pc_df
meta_clean_py <- py_stage1$meta

############################################
## 3. STAGE 2: RAW pY (After purity correction, Before imputation)  
############################################

cat("\n=== pY STAGE 2: Raw + Purity correction ===\n")
py_prot_raw_t <- t(clean_for_pca(pY_data_aligned))
py_prot_raw_pc_t <- removeBatchEffect(py_prot_raw_t, covariates = purity_py)
pY_raw_pc <- t(py_prot_raw_pc_t)

py_stage2 <- do_pca_analysis(pY_raw_pc, "raw_pc", meta_clean_py)
py_pc2_df <- py_stage2$pc_df
meta_clean_py <- py_stage2$meta

############################################
## 4. STAGE 3: IMPUTED pY (After imputation, Before purity correction)
############################################

cat("\n=== pY STAGE 3: Imputed data ===\n")
py_qrilc_result <- impute.QRILC(as.matrix(pY_data_aligned), tune.sigma = 1)
pY_imputed <- py_qrilc_result[[1]]

py_stage3 <- do_pca_analysis(pY_imputed, "imp", meta_clean_py)
py_pc3_df <- py_stage3$pc_df
meta_clean_py <- py_stage3$meta

############################################
## 5. STAGE 4: IMPUTED pY + PURITY CORRECTED (FINAL)
############################################

cat("\n=== pY STAGE 4: Imputed + Purity corrected ===\n")
py_prot_imp_t <- t(clean_for_pca(pY_imputed))
py_prot_imp_pc_t <- removeBatchEffect(py_prot_imp_t, covariates = purity_py)
pY_final <- t(py_prot_imp_pc_t)

py_stage4 <- do_pca_analysis(pY_final, "final", meta_clean_py)
py_pc4_df <- py_stage4$pc_df
meta_clean_py <- py_stage4$meta

cat("pY Pipeline complete!\n")

############################################
## 6. 4-STAGE pY PCA PLOTS
############################################

py_p1 <- ggplot(py_pc1_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "pY 1. RAW", x = "PC1", y = "PC2")

py_p2 <- ggplot(py_pc2_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "pY 2. RAW+PC", x = "PC1", y = "PC2")

py_p3 <- ggplot(py_pc3_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "pY 3. IMP", x = "PC1", y = "PC2")

py_p4 <- ggplot(py_pc4_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "pY 4. FINAL", x = "PC1", y = "PC2")

ggarrange(py_p1, py_p2, py_p3, py_p4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 7. 4-STAGE pY PURITY BOXPLOTS
############################################

py_box1 <- ggboxplot(meta_clean_py, x = "pY_cluster_raw", y = "Purity", 
                    color = "pY_cluster_raw", add = "jitter") + theme_bw() +
  labs(title = "pY 1. RAW")

py_box2 <- ggboxplot(meta_clean_py, x = "pY_cluster_raw_pc", y = "Purity", 
                    color = "pY_cluster_raw_pc", add = "jitter") + theme_bw() +
  labs(title = "pY 2. RAW+PC")

py_box3 <- ggboxplot(meta_clean_py, x = "pY_cluster_imp", y = "Purity", 
                    color = "pY_cluster_imp", add = "jitter") + theme_bw() +
  labs(title = "pY 3. IMP")

py_box4 <- ggboxplot(meta_clean_py, x = "pY_cluster_final", y = "Purity", 
                    color = "pY_cluster_final", add = "jitter") + theme_bw() +
  labs(title = "pY 4. FINAL")

ggarrange(py_box1, py_box2, py_box3, py_box4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 8. pY METRICS & SUMMARY
############################################

evaluate_stage <- function(pca_obj, meta, cluster_col, purity_col = "Purity") {
  pc1_var <- summary(pca_obj)$importance[2,1]*100
  pc2_var <- summary(pca_obj)$importance[2,2]*100
  
  sil_ets <- silhouette(as.numeric(as.factor(meta$ETS)), dist(pca_obj$x[,1:2]))
  ets_sil <- mean(sil_ets[,3])
  
  purity_p <- summary(aov(meta[[purity_col]] ~ meta[[cluster_col]]))[[1]]$`Pr(>F)`[1]
  
  data.frame(PC1_var = pc1_var, PC2_var = pc2_var, ETS_silhouette = ets_sil, Purity_pval = purity_p)
}

py_metrics <- rbind(
  evaluate_stage(py_stage1$pca, meta_clean_py, "pY_cluster_raw"),
  evaluate_stage(py_stage2$pca, meta_clean_py, "pY_cluster_raw_pc"),
  evaluate_stage(py_stage3$pca, meta_clean_py, "pY_cluster_imp"),
  evaluate_stage(py_stage4$pca, meta_clean_py, "pY_cluster_final")
)
rownames(py_metrics) <- c("pY_RAW", "pY_RAW+PC", "pY_IMP", "pY_FINAL")
print("pY Metrics:")
print(py_metrics)

cat("\npY PC1 variances (%):")
cat(sprintf(" RAW: %.1f, RAW+PC: %.1f, IMP: %.1f, FINAL: %.1f\n",
    summary(py_stage1$pca)$importance[2,1]*100,
    summary(py_stage2$pca)$importance[2,1]*100, 
    summary(py_stage3$pca)$importance[2,1]*100,
    summary(py_stage4$pca)$importance[2,1]*100))

############################################
## 9. EXTRACT BEST pY STAGE (RAW+PC) - TOP 50 PCs
############################################

pY_pc_top50 <- py_stage2$pca$x[, 1:50, drop = FALSE]
colnames(pY_pc_top50) <- paste0("pY_PC_", 1:50)
meta_clean_py <- cbind(meta_clean_py, pY_pc_top50)

cat("\n🏆 pY RAW+PC (Stage 2) - Top 50 PCs extracted\n")
cat("pY outlier check - samples with extreme PC1 (>3SD):\n")
pc1_sd <- sd(pY_pc_top50[,1])
pc1_outliers <- rownames(pY_pc_top50)[abs(pY_pc_top50[,1] - mean(pY_pc_top50[,1])) > 3*pc1_sd]
if(length(pc1_outliers) > 0) {
  cat("Outliers:", paste(pc1_outliers, collapse = ", "), "\n")
} else {
  cat("No extreme outliers detected\n")
}

print("pY Key outputs:")
print("- pY_pc_top50: Best stage PCA scores")
print("- pY_raw_pc: Best corrected pY matrix")
print("- meta_clean_py: pY metadata w/ clusters")
print("- py_metrics: Stage comparison table")

write.csv(pY_pc_top50, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/LPM_pY_pc_top50.csv", row.names = TRUE)
# phos_data <- pY_data
# metadata <- meta
# # Ensure case_id alignment
# common_ids <- intersect(rownames(phos_data), metadata$case_id)
# phos_sub   <- phos_data[common_ids, ]
# meta_sub   <- metadata %>% filter(case_id %in% common_ids)
# meta_sub   <- meta_sub[match(rownames(phos_sub), meta_sub$case_id), ]

# # Transpose for Protein-in-Rows format
# expr_raw <- t(phos_sub)
# purity_vec <- meta_sub$Purity

# # -------------------------------------------------------------------------
# # 3. ANALYZE PRE-IMPUTATION (Complete Proteins Only)
# # -------------------------------------------------------------------------
# # We use only proteins with 0 NAs to see clusters without imputation influence
# complete_cases <- rowSums(is.na(expr_raw)) == 0
# expr_complete <- expr_raw[complete_cases, ]

# # Manual Purity Correction (Residuals) for Complete Proteins
# # We subtract the linear effect of Purity from each protein
# expr_complete_corr <- apply(expr_complete, 1, function(y) {
#   residuals(lm(y ~ purity_vec))
# }) %>% t()

# # -------------------------------------------------------------------------
# # 4. ANALYZE POST-IMPUTATION (QRILC + Purity Correction)
# # -------------------------------------------------------------------------
# # Filter for proteins with < 80% missingness and impute
# expr_filtered <- expr_raw[rowMeans(is.na(expr_raw)) < 0.8, ]
# expr_imputed  <- impute.QRILC(as.matrix(expr_filtered))[[1]]

# # Batch effect correction for imputed data
# expr_imputed_corr <- removeBatchEffect(expr_imputed, covariates = purity_vec)

# # -------------------------------------------------------------------------
# # 5. VISUALIZATION & CLUSTER DETECTION
# # -------------------------------------------------------------------------

# # PCA on Complete Proteins (Pre-Imputation)
# pca_complete <- prcomp(t(expr_complete_corr), scale. = TRUE)
# p_comp <- fviz_pca_ind(pca_complete, 
#                        habillage = as.factor(meta_sub$Main_Groups), 
#                        title = "PCA: Complete Proteins Only (No Imputation, Purity Corrected)",
#                        geom = "point", pointsize = 2)

# # PCA on Full Dataset (After QRILC Imputation)
# pca_full <- prcomp(t(expr_imputed_corr), scale. = TRUE)
# p_full <- fviz_pca_ind(pca_full, 
#                        habillage = as.factor(meta_sub$Main_Groups), 
#                        title = "PCA: All Proteins (QRILC Imputed, Purity Corrected)",
#                        geom = "point", pointsize = 2)

# print(p_comp)
# print(p_full)

# # Heatmap Annotation
# anno_col <- data.frame(
#   Main_Groups = meta_sub$Main_Groups,
#   Purity = meta_sub$Purity,
#   row.names = meta_sub$case_id
# )

# # Identifying outliers/clusters using Hierarchical Clustering
# # We cluster samples based on the top 500 variable proteins
# top_vars <- names(sort(apply(expr_imputed_corr, 1, var), decreasing = TRUE))[1:min(500, nrow(expr_imputed_corr))]
# sample_dist <- dist(t(expr_imputed_corr[top_vars, ]))
# sample_hclust <- hclust(sample_dist, method = "ward.D2")

# # Generate Heatmap
# pheatmap(expr_imputed_corr[top_vars, ], 
#          annotation_col = anno_col, 
#          clustering_method = "ward.D2",
#          show_rownames = FALSE, 
#          show_colnames = FALSE,
#          main = "Clustering of Samples (After Correction & Imputation)",
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

# # -------------------------------------------------------------------------
# # 6. IDENTIFYING "DIFFERENT" SAMPLES
# # -------------------------------------------------------------------------
# # Cut the tree into 4 main clusters to see which samples group together
# clusters <- cutree(sample_hclust, k = 4)
# meta_sub$Cluster_Assignment <- as.factor(clusters)

# # Summary of which groups fall into which cluster
# cluster_summary <- table(meta_sub$Main_Groups, meta_sub$Cluster_Assignment)
# print("Distribution of Main_Groups across 4 identified clusters:")
# print(cluster_summary)
```


## 5 - data for the LPM - PTMs:Methylation 2
```{r}
methylation_beta <- read.csv("/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_04_05_23_25/betas_filtered.csv")
methylation_beta_fil <- betas_filtered 
rownames(methylation_beta_fil) <- methylation_beta_fil$X
methylation_beta_fil$X <- NULL

cn <- colnames(methylation_beta_fil)

# remove leading junk up to and including 'Methylation_Array.'
cn_clean <- sub(".*Methylation_Array\\.", "", cn)

colnames(methylation_beta_fil) <- cn_clean
colnames(methylation_beta_fil) <- gsub("\\.idat", "", colnames(methylation_beta_fil))

cn <- colnames(methylation_beta_fil)

meth_cols <- 5:length(cn)
cn[meth_cols] <- sub("_[^_]*$", "", cn[meth_cols])

colnames(methylation_beta_fil) <- cn

methylation_beta_fil <- as.data.frame(t(methylation_beta_fil))
methylation_beta_fil$sample <- ifelse(grepl("Tumor_", rownames(methylation_beta_fil)), "tumor", "normal")
methylation_beta_fil <- dplyr::filter(methylation_beta_fil, sample == "tumor")
methylation_beta_fil$sample <- NULL

rownames(methylation_beta_fil) <- gsub("Tumor_", "", rownames(methylation_beta_fil))
rownames(methylation_beta_fil) <- gsub("_8C01", "", rownames(methylation_beta_fil))
rownames(methylation_beta_fil) <- gsub("_5C01", "", rownames(methylation_beta_fil))
rownames(methylation_beta_fil) <- gsub("_6C01", "", rownames(methylation_beta_fil))
rownames(methylation_beta_fil) <- gsub("\\.", "-", rownames(methylation_beta_fil))


############################################
## Methylation Analysis Libraries
############################################

library(dplyr)
library(limma)
library(ggplot2)
library(ggpubr)
library(cluster)

############################################
## 1. Setup & SAFE PCA FUNCTION (methylation-specific)
############################################

# Align methylation data with metadata (rownames should be case_ids like C3L-07813)
meta_clean_methyl <- meta %>% filter(case_id %in% rownames(methylation_beta_fil))
methyl_aligned <- methylation_beta_fil[meta_clean_methyl$case_id, , drop = FALSE]
purity_methyl <- as.numeric(meta_clean_methyl$Purity)
names(purity_methyl) <- meta_clean_methyl$case_id

cat("Methylation Dataset:", nrow(methyl_aligned), "samples x", ncol(methyl_aligned), "CpG sites\n")

# Check for NAs (methylation typically clean)
cat("NA percentage:", round(mean(is.na(methyl_aligned))*100, 2), "%\n")

clean_for_pca <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if(ncol(mat) == 0 || nrow(mat) == 0) stop("Empty matrix")
  
  # Remove constants  
  vars <- apply(mat, 2, var, na.rm = TRUE)
  keep_cols <- vars > 0 & is.finite(vars)
  mat <- mat[, keep_cols, drop = FALSE]
  
  if(ncol(mat) == 0) stop("No variable columns")
  
  # Simple median imputation (minimal NAs expected)
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  # Final cleanup
  mat[!is.finite(mat)] <- NA
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  return(mat)
}

n_top_methyl <- min(5000, ncol(methyl_aligned))

do_pca_analysis <- function(prot_mat, stage_name, meta_data) {
  prot_clean <- clean_for_pca(prot_mat)
  var_order <- apply(prot_clean, 2, var, na.rm = TRUE)
  top_prot <- names(sort(var_order, decreasing = TRUE))[1:n_top_methyl]
  
  prot_top <- prot_clean[, top_prot, drop = FALSE]
  z_prot <- scale(prot_top, center = TRUE, scale = TRUE)
  pca_obj <- prcomp(z_prot, center = FALSE, scale. = FALSE)
  
  pc_df <- data.frame(pca_obj$x[, 1:2], case_id = rownames(pca_obj$x)) %>%
    left_join(meta_data[, c("case_id", "ETS", "Purity")], by = "case_id")
  
  km <- kmeans(pca_obj$x[, 1:min(10, ncol(pca_obj$x))], centers = 3, nstart = 50)
  meta_data[[paste0("methyl_cluster_", stage_name)]] <- factor(km$cluster)
  
  list(pca = pca_obj, pc_df = pc_df, meta = meta_data)
}

set.seed(123)

############################################
## 2. STAGE 1: RAW methylation
############################################

cat("\n=== METHYLATION STAGE 1: Raw data ===\n")
methyl_stage1 <- do_pca_analysis(methyl_aligned, "raw", meta_clean_methyl)
methyl_pc1_df <- methyl_stage1$pc_df
meta_clean_methyl <- methyl_stage1$meta

############################################
## 3. STAGE 2: RAW methylation + Purity correction (BEST)
############################################

cat("\n=== METHYLATION STAGE 2: Raw + Purity correction ===\n")
methyl_raw_t <- t(clean_for_pca(methyl_aligned))
methyl_raw_pc_t <- removeBatchEffect(methyl_raw_t, covariates = purity_methyl)
methyl_raw_pc <- t(methyl_raw_pc_t)

methyl_stage2 <- do_pca_analysis(methyl_raw_pc, "raw_pc", meta_clean_methyl)
methyl_pc2_df <- methyl_stage2$pc_df
meta_clean_methyl <- methyl_stage2$meta

############################################
## 4. STAGE 3: Median-imputed (for completeness)
############################################

cat("\n=== METHYLATION STAGE 3: Median-imputed ===\n")
methyl_imputed <- clean_for_pca(methyl_aligned)

methyl_stage3 <- do_pca_analysis(methyl_imputed, "imp", meta_clean_methyl)
methyl_pc3_df <- methyl_stage3$pc_df
meta_clean_methyl <- methyl_stage3$meta

############################################
## 5. STAGE 4: Median-imputed + Purity corrected (FINAL) 
############################################

cat("\n=== METHYLATION STAGE 4: Median-imputed + Purity corrected ===\n")
methyl_imp_t <- t(clean_for_pca(methyl_imputed))
methyl_imp_pc_t <- removeBatchEffect(methyl_imp_t, covariates = purity_methyl)
methyl_final <- t(methyl_imp_pc_t)

methyl_stage4 <- do_pca_analysis(methyl_final, "final", meta_clean_methyl)
methyl_pc4_df <- methyl_stage4$pc_df
meta_clean_methyl <- methyl_stage4$meta

cat("Methylation Pipeline complete!\n")

############################################
## 6. 4-STAGE METHYLATION PCA PLOTS
############################################

methyl_p1 <- ggplot(methyl_pc1_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Methyl 1. RAW", x = "PC1", y = "PC2")

methyl_p2 <- ggplot(methyl_pc2_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Methyl 2. RAW+PC", x = "PC1", y = "PC2")

methyl_p3 <- ggplot(methyl_pc3_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Methyl 3. Median-imp", x = "PC1", y = "PC2")

methyl_p4 <- ggplot(methyl_pc4_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Methyl 4. FINAL", x = "PC1", y = "PC2")

ggarrange(methyl_p1, methyl_p2, methyl_p3, methyl_p4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 7. 4-STAGE METHYLATION PURITY BOXPLOTS
############################################

methyl_box1 <- ggboxplot(meta_clean_methyl, x = "methyl_cluster_raw", y = "Purity", 
                        color = "methyl_cluster_raw", add = "jitter") + theme_bw() +
  labs(title = "Methyl 1. RAW")

methyl_box2 <- ggboxplot(meta_clean_methyl, x = "methyl_cluster_raw_pc", y = "Purity", 
                        color = "methyl_cluster_raw_pc", add = "jitter") + theme_bw() +
  labs(title = "Methyl 2. RAW+PC")

methyl_box3 <- ggboxplot(meta_clean_methyl, x = "methyl_cluster_imp", y = "Purity", 
                        color = "methyl_cluster_imp", add = "jitter") + theme_bw() +
  labs(title = "Methyl 3. Median-imp")

methyl_box4 <- ggboxplot(meta_clean_methyl, x = "methyl_cluster_final", y = "Purity", 
                        color = "methyl_cluster_final", add = "jitter") + theme_bw() +
  labs(title = "Methyl 4. FINAL")

ggarrange(methyl_box1, methyl_box2, methyl_box3, methyl_box4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 8. METHYLATION METRICS & SUMMARY
############################################

evaluate_stage <- function(pca_obj, meta, cluster_col, purity_col = "Purity") {
  pc1_var <- summary(pca_obj)$importance[2,1]*100
  pc2_var <- summary(pca_obj)$importance[2,2]*100
  
  sil_ets <- silhouette(as.numeric(as.factor(meta$ETS)), dist(pca_obj$x[,1:2]))
  ets_sil <- mean(sil_ets[,3])
  
  purity_p <- summary(aov(meta[[purity_col]] ~ meta[[cluster_col]]))[[1]]$`Pr(>F)`[1]
  
  data.frame(PC1_var = pc1_var, PC2_var = pc2_var, ETS_silhouette = ets_sil, Purity_pval = purity_p)
}

methyl_metrics <- rbind(
  evaluate_stage(methyl_stage1$pca, meta_clean_methyl, "methyl_cluster_raw"),
  evaluate_stage(methyl_stage2$pca, meta_clean_methyl, "methyl_cluster_raw_pc"),
  evaluate_stage(methyl_stage3$pca, meta_clean_methyl, "methyl_cluster_imp"),
  evaluate_stage(methyl_stage4$pca, meta_clean_methyl, "methyl_cluster_final")
)
rownames(methyl_metrics) <- c("METHYL_RAW", "METHYL_RAW+PC", "METHYL_MEDIAN", "METHYL_FINAL")
print("Methylation Metrics:")
print(methyl_metrics)

cat("\nMethylation PC1 variances (%):")
cat(sprintf(" RAW: %.1f, RAW+PC: %.1f, MEDIAN: %.1f, FINAL: %.1f\n",
    summary(methyl_stage1$pca)$importance[2,1]*100,
    summary(methyl_stage2$pca)$importance[2,1]*100, 
    summary(methyl_stage3$pca)$importance[2,1]*100,
    summary(methyl_stage4$pca)$importance[2,1]*100))

############################################
## 9. EXTRACT BEST METHYLATION STAGE (RAW+PC) - TOP 50 PCs
############################################

methyl_pc_top50 <- methyl_stage2$pca$x[, 1:50, drop = FALSE]
colnames(methyl_pc_top50) <- paste0("methyl_PC_", 1:50)
meta_clean_methyl <- cbind(meta_clean_methyl, methyl_pc_top50)

cat("\n🏆 METHYLATION RAW+PC (Stage 2) - Top 50 PCs extracted\n")
cat("Methylation outlier check - samples with extreme PC1 (>3SD):\n")
pc1_sd <- sd(methyl_pc_top50[,1])
pc1_outliers <- rownames(methyl_pc_top50)[abs(methyl_pc_top50[,1] - mean(methyl_pc_top50[,1])) > 3*pc1_sd]
if(length(pc1_outliers) > 0) {
  cat("Outliers:", paste(pc1_outliers, collapse = ", "), "\n")
} else {
  cat("No extreme outliers detected\n")
}

print("Methylation Key outputs:")
print("- methyl_pc_top50: Best stage PCA scores")
print("- methyl_raw_pc: Best corrected methyl matrix") 
print("- meta_clean_methyl: Methyl metadata w/ clusters")
print("- methyl_metrics: Stage comparison table")

# Save outputs
write.csv(methyl_pc_top50, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/LPM_methyl_pc_top50.csv", row.names = TRUE)

```


## 6 - data for the LPM - PTMs:Phosphorylation 3
```{r}
phospho <- read.csv("/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/phospho_abundance_gene_MD.tsv", sep = "\t")
phospho <- phospho[,-c(2:5)]
rownames(phospho) <- phospho$Index
phospho$Index <- NULL

# Transpose: samples as rows, genes as columns
phospho_data <- as.data.frame(t(as.matrix(phospho)))

# Filter NCI7 + tumor only
phospho_data$remove <- ifelse(grepl("NCI7", rownames(phospho_data)), "NCI7", "keep")
phospho_data <- dplyr::filter(phospho_data, remove != "NCI7")
phospho_data$remove <- NULL
phospho_data$sample <- ifelse(grepl(".T", rownames(phospho_data)), "tumor", "normal")
phospho_data <- dplyr::filter(phospho_data, sample == "tumor")
phospho_data$sample <- NULL
rownames(phospho_data) <- gsub(".T", "", rownames(phospho_data))
rownames(phospho_data) <- gsub("\\.", "-", rownames(phospho_data))

# Map Ensembl → Gene symbols
ens_ids <- colnames(phospho_data)
gene_symbols <- mapIds(org.Hs.eg.db, keys = ens_ids, column = "SYMBOL", 
                      keytype = "ENSEMBL", multiVals = "first")
gene_symbols[is.na(gene_symbols)] <- ens_ids[is.na(gene_symbols)]
colnames(phospho_data) <- gene_symbols

cat("Raw phospho data:", nrow(phospho_data), "samples x", ncol(phospho_data), "sites\n")

############################################
## Phospho Analysis Libraries
############################################

library(dplyr)
library(imputeLCMD)
library(limma)
library(ggplot2)
library(ggpubr)
library(cluster)

############################################
## 1. Setup & SAFE PCA FUNCTION (phospho-specific)
############################################

meta_clean_phospho <- meta %>% filter(case_id %in% rownames(phospho_data))
phospho_aligned <- phospho_data[meta_clean_phospho$case_id, , drop = FALSE]
purity_phospho <- as.numeric(meta_clean_phospho$Purity)
names(purity_phospho) <- meta_clean_phospho$case_id

cat("Phospho Dataset:", nrow(phospho_aligned), "samples x", ncol(phospho_aligned), "phospho sites\n")

clean_for_pca <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if(ncol(mat) == 0 || nrow(mat) == 0) stop("Empty matrix")
  
  vars <- apply(mat, 2, var, na.rm = TRUE)
  keep_cols <- vars > 0 & is.finite(vars)
  mat <- mat[, keep_cols, drop = FALSE]
  
  if(ncol(mat) == 0) stop("No variable columns")
  
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  mat[!is.finite(mat)] <- NA
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  return(mat)
}

n_top_phospho <- min(5000, ncol(phospho_aligned))

do_pca_analysis <- function(prot_mat, stage_name, meta_data) {
  prot_clean <- clean_for_pca(prot_mat)
  var_order <- apply(prot_clean, 2, var, na.rm = TRUE)
  top_prot <- names(sort(var_order, decreasing = TRUE))[1:n_top_phospho]
  
  prot_top <- prot_clean[, top_prot, drop = FALSE]
  z_prot <- scale(prot_top, center = TRUE, scale = TRUE)
  pca_obj <- prcomp(z_prot, center = FALSE, scale. = FALSE)
  
  pc_df <- data.frame(pca_obj$x[, 1:2], case_id = rownames(pca_obj$x)) %>%
    left_join(meta_data[, c("case_id", "ETS", "Purity")], by = "case_id")
  
  km <- kmeans(pca_obj$x[, 1:min(10, ncol(pca_obj$x))], centers = 3, nstart = 50)
  meta_data[[paste0("phospho_cluster_", stage_name)]] <- factor(km$cluster)
  
  list(pca = pca_obj, pc_df = pc_df, meta = meta_data)
}

set.seed(123)

############################################
## 2. STAGE 1: RAW phospho
############################################

cat("\n=== PHOSPHo STAGE 1: Raw data ===\n")
phospho_stage1 <- do_pca_analysis(phospho_aligned, "raw", meta_clean_phospho)
phospho_pc1_df <- phospho_stage1$pc_df
meta_clean_phospho <- phospho_stage1$meta

############################################
## 3. STAGE 2: RAW phospho + Purity correction  
############################################

cat("\n=== PHOSPHo STAGE 2: Raw + Purity correction ===\n")
phospho_raw_t <- t(clean_for_pca(phospho_aligned))
phospho_raw_pc_t <- removeBatchEffect(phospho_raw_t, covariates = purity_phospho)
phospho_raw_pc <- t(phospho_raw_pc_t)

phospho_stage2 <- do_pca_analysis(phospho_raw_pc, "raw_pc", meta_clean_phospho)
phospho_pc2_df <- phospho_stage2$pc_df
meta_clean_phospho <- phospho_stage2$meta

############################################
## 4. STAGE 3: IMPUTED phospho
############################################

cat("\n=== PHOSPHo STAGE 3: Imputed data ===\n")
phospho_qrilc <- impute.QRILC(as.matrix(phospho_aligned), tune.sigma = 1)
phospho_imputed <- phospho_qrilc[[1]]

phospho_stage3 <- do_pca_analysis(phospho_imputed, "imp", meta_clean_phospho)
phospho_pc3_df <- phospho_stage3$pc_df
meta_clean_phospho <- phospho_stage3$meta

############################################
## 5. STAGE 4: IMPUTED phospho + Purity corrected
############################################

cat("\n=== PHOSPHo STAGE 4: Imputed + Purity corrected ===\n")
phospho_imp_t <- t(clean_for_pca(phospho_imputed))
phospho_imp_pc_t <- removeBatchEffect(phospho_imp_t, covariates = purity_phospho)
phospho_final <- t(phospho_imp_pc_t)

phospho_stage4 <- do_pca_analysis(phospho_final, "final", meta_clean_phospho)
phospho_pc4_df <- phospho_stage4$pc_df
meta_clean_phospho <- phospho_stage4$meta

cat("Phospho Pipeline complete!\n")

############################################
## 6. 4-STAGE PHOSPHo PCA PLOTS
############################################

phos_p1 <- ggplot(phospho_pc1_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Phospho 1. RAW", x = "PC1", y = "PC2")

phos_p2 <- ggplot(phospho_pc2_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Phospho 2. RAW+PC", x = "PC1", y = "PC2")

phos_p3 <- ggplot(phospho_pc3_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Phospho 3. IMP", x = "PC1", y = "PC2")

phos_p4 <- ggplot(phospho_pc4_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Phospho 4. FINAL", x = "PC1", y = "PC2")

ggarrange(phos_p1, phos_p2, phos_p3, phos_p4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 7. 4-STAGE PHOSPHo PURITY BOXPLOTS
############################################

phos_box1 <- ggboxplot(meta_clean_phospho, x = "phospho_cluster_raw", y = "Purity", 
                      color = "phospho_cluster_raw", add = "jitter") + theme_bw() +
  labs(title = "Phospho 1. RAW")

phos_box2 <- ggboxplot(meta_clean_phospho, x = "phospho_cluster_raw_pc", y = "Purity", 
                      color = "phospho_cluster_raw_pc", add = "jitter") + theme_bw() +
  labs(title = "Phospho 2. RAW+PC")

phos_box3 <- ggboxplot(meta_clean_phospho, x = "phospho_cluster_imp", y = "Purity", 
                      color = "phospho_cluster_imp", add = "jitter") + theme_bw() +
  labs(title = "Phospho 3. IMP")

phos_box4 <- ggboxplot(meta_clean_phospho, x = "phospho_cluster_final", y = "Purity", 
                      color = "phospho_cluster_final", add = "jitter") + theme_bw() +
  labs(title = "Phospho 4. FINAL")

ggarrange(phos_box1, phos_box2, phos_box3, phos_box4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 8. PHOSPHo METRICS & SUMMARY
############################################

evaluate_stage <- function(pca_obj, meta, cluster_col, purity_col = "Purity") {
  pc1_var <- summary(pca_obj)$importance[2,1]*100
  pc2_var <- summary(pca_obj)$importance[2,2]*100
  
  sil_ets <- silhouette(as.numeric(as.factor(meta$ETS)), dist(pca_obj$x[,1:2]))
  ets_sil <- mean(sil_ets[,3])
  
  purity_p <- summary(aov(meta[[purity_col]] ~ meta[[cluster_col]]))[[1]]$`Pr(>F)`[1]
  
  data.frame(PC1_var = pc1_var, PC2_var = pc2_var, ETS_silhouette = ets_sil, Purity_pval = purity_p)
}

phospho_metrics <- rbind(
  evaluate_stage(phospho_stage1$pca, meta_clean_phospho, "phospho_cluster_raw"),
  evaluate_stage(phospho_stage2$pca, meta_clean_phospho, "phospho_cluster_raw_pc"),
  evaluate_stage(phospho_stage3$pca, meta_clean_phospho, "phospho_cluster_imp"),
  evaluate_stage(phospho_stage4$pca, meta_clean_phospho, "phospho_cluster_final")
)
rownames(phospho_metrics) <- c("PHOSPHo_RAW", "PHOSPHo_RAW+PC", "PHOSPHo_IMP", "PHOSPHo_FINAL")
print("Phospho Metrics:")
print(phospho_metrics)

cat("\nPhospho PC1 variances (%):")
cat(sprintf(" RAW: %.1f, RAW+PC: %.1f, IMP: %.1f, FINAL: %.1f\n",
    summary(phospho_stage1$pca)$importance[2,1]*100,
    summary(phospho_stage2$pca)$importance[2,1]*100, 
    summary(phospho_stage3$pca)$importance[2,1]*100,
    summary(phospho_stage4$pca)$importance[2,1]*100))

############################################
## 9. EXTRACT BEST PHOSPHo STAGE (RAW+PC) - TOP 50 PCs
############################################

phospho_pc_top50 <- phospho_stage2$pca$x[, 1:50, drop = FALSE]
colnames(phospho_pc_top50) <- paste0("phospho_PC_", 1:50)
meta_clean_phospho <- cbind(meta_clean_phospho, phospho_pc_top50)

cat("\n🏆 Phospho RAW+PC (Stage 2) - Top 50 PCs extracted\n")
cat("Phospho outlier check - samples with extreme PC1 (>3SD):\n")
pc1_sd <- sd(phospho_pc_top50[,1])
pc1_outliers <- rownames(phospho_pc_top50)[abs(phospho_pc_top50[,1] - mean(phospho_pc_top50[,1])) > 3*pc1_sd]
if(length(pc1_outliers) > 0) {
  cat("Outliers:", paste(pc1_outliers, collapse = ", "), "\n")
} else {
  cat("No extreme outliers detected\n")
}

print("Phospho Key outputs:")
print("- phospho_pc_top50: Best stage PCA scores")
print("- phospho_raw_pc: Best corrected phospho matrix") 
print("- meta_clean_phospho: Phospho metadata w/ clusters")
print("- phospho_metrics: Stage comparison table")

write.csv(phospho_pc_top50, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/LPM_phospho_pc_top50.csv", row.names = TRUE)
```

## 7 - data for the LPM - PTMs:Acetylation 4
```{r}
############################################
## Acetylation Analysis: 4-Stage Pipeline
############################################

# 0. PREPROCESSING (acetyl instead of pY)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

acetyl <- read.csv("/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/acetyl_abundance_gene_MD.tsv", sep = "\t")
acetyl <- acetyl[,-c(2:5)]
rownames(acetyl) <- acetyl$Index
acetyl$Index <- NULL

# Transpose: samples as rows, genes as columns
acetyl_data <- as.data.frame(t(as.matrix(acetyl)))

# Filter NCI7 + tumor only
acetyl_data$remove <- ifelse(grepl("NCI7", rownames(acetyl_data)), "NCI7", "keep")
acetyl_data <- dplyr::filter(acetyl_data, remove != "NCI7")
acetyl_data$remove <- NULL
acetyl_data$sample <- ifelse(grepl(".T", rownames(acetyl_data)), "tumor", "normal")
acetyl_data <- dplyr::filter(acetyl_data, sample == "tumor")
acetyl_data$sample <- NULL
rownames(acetyl_data) <- gsub(".T", "", rownames(acetyl_data))
rownames(acetyl_data) <- gsub("\\.", "-", rownames(acetyl_data))

# Map Ensembl → Gene symbols
ens_ids <- colnames(acetyl_data)
gene_symbols <- mapIds(org.Hs.eg.db, keys = ens_ids, column = "SYMBOL", 
                      keytype = "ENSEMBL", multiVals = "first")
gene_symbols[is.na(gene_symbols)] <- ens_ids[is.na(gene_symbols)]
colnames(acetyl_data) <- gene_symbols

cat("Raw acetyl data:", nrow(acetyl_data), "samples x", ncol(acetyl_data), "sites\n")

############################################
## Acetyl Analysis Libraries
############################################

library(imputeLCMD)
library(limma)
library(ggplot2)
library(ggpubr)
library(cluster)

############################################
## 1. Setup & SAFE PCA FUNCTION (acetyl-specific)
############################################

meta_clean_acetyl <- meta %>% filter(case_id %in% rownames(acetyl_data))
acetyl_data_aligned <- acetyl_data[meta_clean_acetyl$case_id, , drop = FALSE]
purity_acetyl <- as.numeric(meta_clean_acetyl$Purity)
names(purity_acetyl) <- meta_clean_acetyl$case_id

cat("Acetyl Dataset:", nrow(acetyl_data_aligned), "samples x", ncol(acetyl_data_aligned), "acetyl sites\n")

clean_for_pca <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if(ncol(mat) == 0 || nrow(mat) == 0) stop("Empty matrix")
  
  # Remove constants
  vars <- apply(mat, 2, var, na.rm = TRUE)
  keep_cols <- vars > 0 & is.finite(vars)
  mat <- mat[, keep_cols, drop = FALSE]
  
  if(ncol(mat) == 0) stop("No variable columns")
  
  # NA → median
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  # Final cleanup
  mat[!is.finite(mat)] <- NA
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  return(mat)
}

n_top_acetyl <- min(5000, ncol(acetyl_data_aligned))

do_pca_analysis <- function(prot_mat, stage_name, meta_data) {
  prot_clean <- clean_for_pca(prot_mat)
  var_order <- apply(prot_clean, 2, var, na.rm = TRUE)
  top_prot <- names(sort(var_order, decreasing = TRUE))[1:n_top_acetyl]
  
  prot_top <- prot_clean[, top_prot, drop = FALSE]
  z_prot <- scale(prot_top, center = TRUE, scale = TRUE)
  pca_obj <- prcomp(z_prot, center = FALSE, scale. = FALSE)
  
  pc_df <- data.frame(pca_obj$x[, 1:2], case_id = rownames(pca_obj$x)) %>%
    left_join(meta_data[, c("case_id", "ETS", "Purity")], by = "case_id")
  
  km <- kmeans(pca_obj$x[, 1:min(10, ncol(pca_obj$x))], centers = 4, nstart = 50)
  meta_data[[paste0("acetyl_cluster_", stage_name)]] <- factor(km$cluster)
  
  list(pca = pca_obj, pc_df = pc_df, meta = meta_data)
}

set.seed(123)

############################################
## 2. STAGE 1: RAW acetyl (Before imputation, Before purity correction)
############################################

cat("\n=== ACETYL STAGE 1: Raw data ===\n")
acetyl_stage1 <- do_pca_analysis(acetyl_data_aligned, "raw", meta_clean_acetyl)
acetyl_pc1_df <- acetyl_stage1$pc_df
meta_clean_acetyl <- acetyl_stage1$meta

############################################
## 3. STAGE 2: RAW acetyl (After purity correction, Before imputation)  
############################################

cat("\n=== ACETYL STAGE 2: Raw + Purity correction ===\n")
acetyl_prot_raw_t <- t(clean_for_pca(acetyl_data_aligned))
acetyl_prot_raw_pc_t <- removeBatchEffect(acetyl_prot_raw_t, covariates = purity_acetyl)
acetyl_raw_pc <- t(acetyl_prot_raw_pc_t)

acetyl_stage2 <- do_pca_analysis(acetyl_raw_pc, "raw_pc", meta_clean_acetyl)
acetyl_pc2_df <- acetyl_stage2$pc_df
meta_clean_acetyl <- acetyl_stage2$meta

############################################
## 4. STAGE 3: IMPUTED acetyl (After imputation, Before purity correction)
############################################

cat("\n=== ACETYL STAGE 3: Imputed data ===\n")
acetyl_qrilc_result <- impute.QRILC(as.matrix(acetyl_data_aligned), tune.sigma = 1)
acetyl_imputed <- acetyl_qrilc_result[[1]]

acetyl_stage3 <- do_pca_analysis(acetyl_imputed, "imp", meta_clean_acetyl)
acetyl_pc3_df <- acetyl_stage3$pc_df
meta_clean_acetyl <- acetyl_stage3$meta

############################################
## 5. STAGE 4: IMPUTED acetyl + PURITY CORRECTED (FINAL)
############################################

cat("\n=== ACETYL STAGE 4: Imputed + Purity corrected ===\n")
acetyl_prot_imp_t <- t(clean_for_pca(acetyl_imputed))
acetyl_prot_imp_pc_t <- removeBatchEffect(acetyl_prot_imp_t, covariates = purity_acetyl)
acetyl_final <- t(acetyl_prot_imp_pc_t)

acetyl_stage4 <- do_pca_analysis(acetyl_final, "final", meta_clean_acetyl)
acetyl_pc4_df <- acetyl_stage4$pc_df
meta_clean_acetyl <- acetyl_stage4$meta

cat("Acetyl Pipeline complete!\n")

############################################
## 6. 4-STAGE ACETYL PCA PLOTS
############################################

acetyl_p1 <- ggplot(acetyl_pc1_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Acetyl 1. RAW", x = "PC1", y = "PC2")

acetyl_p2 <- ggplot(acetyl_pc2_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Acetyl 2. RAW+PC", x = "PC1", y = "PC2")

acetyl_p3 <- ggplot(acetyl_pc3_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Acetyl 3. IMP", x = "PC1", y = "PC2")

acetyl_p4 <- ggplot(acetyl_pc4_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Acetyl 4. FINAL", x = "PC1", y = "PC2")

ggarrange(acetyl_p1, acetyl_p2, acetyl_p3, acetyl_p4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 7. 4-STAGE ACETYL PURITY BOXPLOTS
############################################

acetyl_box1 <- ggboxplot(meta_clean_acetyl, x = "acetyl_cluster_raw", y = "Purity",
                        color = "acetyl_cluster_raw", add = "jitter") + theme_bw() +
  labs(title = "Acetyl 1. RAW")

acetyl_box2 <- ggboxplot(meta_clean_acetyl, x = "acetyl_cluster_raw_pc", y = "Purity",
                        color = "acetyl_cluster_raw_pc", add = "jitter") + theme_bw() +
  labs(title = "Acetyl 2. RAW+PC")

acetyl_box3 <- ggboxplot(meta_clean_acetyl, x = "acetyl_cluster_imp", y = "Purity",
                        color = "acetyl_cluster_imp", add = "jitter") + theme_bw() +
  labs(title = "Acetyl 3. IMP")

acetyl_box4 <- ggboxplot(meta_clean_acetyl, x = "acetyl_cluster_final", y = "Purity",
                        color = "acetyl_cluster_final", add = "jitter") + theme_bw() +
  labs(title = "Acetyl 4. FINAL")

ggarrange(acetyl_box1, acetyl_box2, acetyl_box3, acetyl_box4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 8. ACETYL METRICS & SUMMARY
############################################

evaluate_stage <- function(pca_obj, meta, cluster_col, purity_col = "Purity") {
  pc1_var <- summary(pca_obj)$importance[2,1]*100
  pc2_var <- summary(pca_obj)$importance[2,2]*100
  
  sil_ets <- silhouette(as.numeric(as.factor(meta$ETS)), dist(pca_obj$x[,1:2]))
  ets_sil <- mean(sil_ets[,3])
  
  purity_p <- summary(aov(meta[[purity_col]] ~ meta[[cluster_col]]))[[1]]$`Pr(>F)`[1]
  
  data.frame(PC1_var = pc1_var, PC2_var = pc2_var, ETS_silhouette = ets_sil, Purity_pval = purity_p)
}

acetyl_metrics <- rbind(
  evaluate_stage(acetyl_stage1$pca, meta_clean_acetyl, "acetyl_cluster_raw"),
  evaluate_stage(acetyl_stage2$pca, meta_clean_acetyl, "acetyl_cluster_raw_pc"),
  evaluate_stage(acetyl_stage3$pca, meta_clean_acetyl, "acetyl_cluster_imp"),
  evaluate_stage(acetyl_stage4$pca, meta_clean_acetyl, "acetyl_cluster_final")
)
rownames(acetyl_metrics) <- c("ACETYL_RAW", "ACETYL_RAW+PC", "ACETYL_IMP", "ACETYL_FINAL")
print("Acetyl Metrics:")
print(acetyl_metrics)

cat("\nAcetyl PC1 variances (%):")
cat(sprintf(" RAW: %.1f, RAW+PC: %.1f, IMP: %.1f, FINAL: %.1f\n",
    summary(acetyl_stage1$pca)$importance[2,1]*100,
    summary(acetyl_stage2$pca)$importance[2,1]*100, 
    summary(acetyl_stage3$pca)$importance[2,1]*100,
    summary(acetyl_stage4$pca)$importance[2,1]*100))

############################################
## 9. EXTRACT BEST ACETYL STAGE (RAW+PC) - TOP 50 PCs
############################################

acetyl_pc_top50 <- acetyl_stage2$pca$x[, 1:50, drop = FALSE]
colnames(acetyl_pc_top50) <- paste0("acetyl_PC_", 1:50)
meta_clean_acetyl <- cbind(meta_clean_acetyl, acetyl_pc_top50)

cat("\n🏆 ACETYL RAW+PC (Stage 2) - Top 50 PCs extracted\n")
cat("Acetyl outlier check - samples with extreme PC1 (>3SD):\n")
pc1_sd <- sd(acetyl_pc_top50[,1])
pc1_outliers <- rownames(acetyl_pc_top50)[abs(acetyl_pc_top50[,1] - mean(acetyl_pc_top50[,1])) > 3*pc1_sd]
if(length(pc1_outliers) > 0) {
  cat("Outliers:", paste(pc1_outliers, collapse = ", "), "\n")
} else {
  cat("No extreme outliers detected\n")
}

print("Acetyl Key outputs:")
print("- acetyl_pc_top50: Best stage PCA scores")
print("- acetyl_raw_pc: Best corrected acetyl matrix")
print("- meta_clean_acetyl: Acetyl metadata w/ clusters")
print("- acetyl_metrics: Stage comparison table")

write.csv(acetyl_pc_top50, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/LPM_acetyl_pc_top50.csv", row.names = TRUE)
```

## 8 - data for the LPM - PTMs:Ubiquitination 5
```{r}
ubiq_DIA <- read.csv("/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_30_12_17_25/prostate_CPTAC_ubiq_DIA_abundance_single_site_MS2quant_Norm.tsv", sep = "\t")
ubiq_DIA <- as.data.frame(ubiq_DIA)
ubiq_DIA <- ubiq_DIA[,-c(1,3)]

rownames(ubiq_DIA) <- ubiq_DIA$Index
ubiq_DIA$Index <- NULL

# Transpose: samples as rows, genes as columns
ubiq_DIA_data <- as.data.frame(t(as.matrix(ubiq_DIA)))

# Filter NCI7 + tumor only
ubiq_DIA_data$remove <- ifelse(grepl("NCI7", rownames(ubiq_DIA_data)), "NCI7", "keep")
ubiq_DIA_data <- dplyr::filter(ubiq_DIA_data, remove != "NCI7")
ubiq_DIA_data$remove <- NULL
ubiq_DIA_data$sample <- ifelse(grepl(".T", rownames(ubiq_DIA_data)), "tumor", "normal")
ubiq_DIA_data <- dplyr::filter(ubiq_DIA_data, sample == "tumor")
ubiq_DIA_data$sample <- NULL
rownames(ubiq_DIA_data) <- gsub(".T", "", rownames(ubiq_DIA_data))
rownames(ubiq_DIA_data) <- gsub("\\.", "-", rownames(ubiq_DIA_data))

############################################
## 1. Setup & SAFE PCA FUNCTION (ubiq-specific)
############################################

meta_clean_ubiq <- meta %>% filter(case_id %in% rownames(ubiq_DIA_data))
ubiq_data_aligned <- ubiq_DIA_data[meta_clean_ubiq$case_id, , drop = FALSE]
purity_ubiq <- as.numeric(meta_clean_ubiq$Purity)
names(purity_ubiq) <- meta_clean_ubiq$case_id

cat("Ubiq Dataset:", nrow(ubiq_data_aligned), "samples x", ncol(ubiq_data_aligned), "ubiq sites\n")

clean_for_pca <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if(ncol(mat) == 0 || nrow(mat) == 0) stop("Empty matrix")
  
  # Remove constants
  vars <- apply(mat, 2, var, na.rm = TRUE)
  keep_cols <- vars > 0 & is.finite(vars)
  mat <- mat[, keep_cols, drop = FALSE]
  
  if(ncol(mat) == 0) stop("No variable columns")
  
  # Simple median imputation
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[, i])
    if(any(nas)) mat[nas, i] <- median(mat[, i], na.rm = TRUE)
  }
  
  # Final cleanup for any newly introduced infinities
  mat[!is.finite(mat)] <- NA
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[, i])
    if(any(nas)) mat[nas, i] <- median(mat[, i], na.rm = TRUE)
  }
  
  return(mat)
}

n_top_ubiq <- min(5000, ncol(ubiq_data_aligned))

do_pca_analysis <- function(prot_mat, stage_name, meta_data) {
  prot_clean <- clean_for_pca(prot_mat)
  var_order <- apply(prot_clean, 2, var, na.rm = TRUE)
  top_prot <- names(sort(var_order, decreasing = TRUE))[1:n_top_ubiq]
  
  prot_top <- prot_clean[, top_prot, drop = FALSE]
  z_prot <- scale(prot_top, center = TRUE, scale = TRUE)
  pca_obj <- prcomp(z_prot, center = FALSE, scale. = FALSE)
  
  pc_df <- data.frame(pca_obj$x[, 1:2], case_id = rownames(pca_obj$x)) %>%
    left_join(meta_data[, c("case_id", "ETS", "Purity")], by = "case_id")
  
  km <- kmeans(pca_obj$x[, 1:min(10, ncol(pca_obj$x))], centers = 4, nstart = 50)
  meta_data[[paste0("ubiq_cluster_", stage_name)]] <- factor(km$cluster)
  
  list(pca = pca_obj, pc_df = pc_df, meta = meta_data)
}

set.seed(123)

############################################
## 2. STAGE 1: RAW ubiq (Before purity correction)
############################################

cat("\n=== UBIQ STAGE 1: Raw data ===\n")
ubiq_stage1 <- do_pca_analysis(ubiq_data_aligned, "raw", meta_clean_ubiq)
ubiq_pc1_df <- ubiq_stage1$pc_df
meta_clean_ubiq <- ubiq_stage1$meta

############################################
## 3. STAGE 2: RAW + Purity correction
############################################

cat("\n=== UBIQ STAGE 2: Raw + Purity correction ===\n")
ubiq_prot_raw_t <- t(clean_for_pca(ubiq_data_aligned))
ubiq_prot_raw_pc_t <- removeBatchEffect(ubiq_prot_raw_t, covariates = purity_ubiq)
ubiq_raw_pc <- t(ubiq_prot_raw_pc_t)

ubiq_stage2 <- do_pca_analysis(ubiq_raw_pc, "raw_pc", meta_clean_ubiq)
ubiq_pc2_df <- ubiq_stage2$pc_df
meta_clean_ubiq <- ubiq_stage2$meta

############################################
## 4. STAGE 3: Median-imputed (for completeness, same method)
############################################

cat("\n=== UBIQ STAGE 3: Median-imputed ===\n")
ubiq_imputed <- clean_for_pca(ubiq_data_aligned)

ubiq_stage3 <- do_pca_analysis(ubiq_imputed, "imp", meta_clean_ubiq)
ubiq_pc3_df <- ubiq_stage3$pc_df
meta_clean_ubiq <- ubiq_stage3$meta

############################################
## 5. STAGE 4: Median-imputed + Purity corrected (FINAL)
############################################

cat("\n=== UBIQ STAGE 4: Median-imputed + Purity corrected ===\n")
ubiq_prot_imp_t <- t(clean_for_pca(ubiq_imputed))
ubiq_prot_imp_pc_t <- removeBatchEffect(ubiq_prot_imp_t, covariates = purity_ubiq)
ubiq_final <- t(ubiq_prot_imp_pc_t)

ubiq_stage4 <- do_pca_analysis(ubiq_final, "final", meta_clean_ubiq)
ubiq_pc4_df <- ubiq_stage4$pc_df
meta_clean_ubiq <- ubiq_stage4$meta

cat("Ubiq Pipeline complete!\n")

############################################
## 6. 4-STAGE UBIQ PCA PLOTS
############################################

ubiq_p1 <- ggplot(ubiq_pc1_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Ubiq 1. RAW", x = "PC1", y = "PC2")

ubiq_p2 <- ggplot(ubiq_pc2_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Ubiq 2. RAW+PC", x = "PC1", y = "PC2")

ubiq_p3 <- ggplot(ubiq_pc3_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Ubiq 3. Median-imputed", x = "PC1", y = "PC2")

ubiq_p4 <- ggplot(ubiq_pc4_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Ubiq 4. FINAL", x = "PC1", y = "PC2")

ggarrange(ubiq_p1, ubiq_p2, ubiq_p3, ubiq_p4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 7. 4-STAGE UBIQ PURITY BOXPLOTS
############################################

ubiq_box1 <- ggboxplot(meta_clean_ubiq, x = "ubiq_cluster_raw", y = "Purity",
                      color = "ubiq_cluster_raw", add = "jitter") + theme_bw() +
  labs(title = "Ubiq 1. RAW")

ubiq_box2 <- ggboxplot(meta_clean_ubiq, x = "ubiq_cluster_raw_pc", y = "Purity",
                      color = "ubiq_cluster_raw_pc", add = "jitter") + theme_bw() +
  labs(title = "Ubiq 2. RAW+PC")

ubiq_box3 <- ggboxplot(meta_clean_ubiq, x = "ubiq_cluster_imp", y = "Purity",
                      color = "ubiq_cluster_imp", add = "jitter") + theme_bw() +
  labs(title = "Ubiq 3. Median-imputed")

ubiq_box4 <- ggboxplot(meta_clean_ubiq, x = "ubiq_cluster_final", y = "Purity",
                      color = "ubiq_cluster_final", add = "jitter") + theme_bw() +
  labs(title = "Ubiq 4. FINAL")

ggarrange(ubiq_box1, ubiq_box2, ubiq_box3, ubiq_box4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 8. UBIQ METRICS & SUMMARY
############################################

evaluate_stage <- function(pca_obj, meta, cluster_col, purity_col = "Purity") {
  pc1_var <- summary(pca_obj)$importance[2,1]*100
  pc2_var <- summary(pca_obj)$importance[2,2]*100
  
  sil_ets <- silhouette(as.numeric(as.factor(meta$ETS)), dist(pca_obj$x[,1:2]))
  ets_sil <- mean(sil_ets[,3])
  
  purity_p <- summary(aov(meta[[purity_col]] ~ meta[[cluster_col]]))[[1]]$`Pr(>F)`[1]
  
  data.frame(PC1_var = pc1_var, PC2_var = pc2_var, ETS_silhouette = ets_sil, Purity_pval = purity_p)
}

ubiq_metrics <- rbind(
  evaluate_stage(ubiq_stage1$pca, meta_clean_ubiq, "ubiq_cluster_raw"),
  evaluate_stage(ubiq_stage2$pca, meta_clean_ubiq, "ubiq_cluster_raw_pc"),
  evaluate_stage(ubiq_stage3$pca, meta_clean_ubiq, "ubiq_cluster_imp"),
  evaluate_stage(ubiq_stage4$pca, meta_clean_ubiq, "ubiq_cluster_final")
)
rownames(ubiq_metrics) <- c("UBIQ_RAW", "UBIQ_RAW+PC", "UBIQ_MEDIAN", "UBIQ_FINAL")
print("Ubiq Metrics:")
print(ubiq_metrics)

cat("\nUbiq PC1 variances (%):")
cat(sprintf(" RAW: %.1f, RAW+PC: %.1f, MEDIAN: %.1f, FINAL: %.1f\n",
    summary(ubiq_stage1$pca)$importance[2,1]*100,
    summary(ubiq_stage2$pca)$importance[2,1]*100, 
    summary(ubiq_stage3$pca)$importance[2,1]*100,
    summary(ubiq_stage4$pca)$importance[2,1]*100))

############################################
## 9. EXTRACT BEST UBIQ STAGE (RAW+PC) - TOP 50 PCs
############################################

ubiq_pc_top50 <- ubiq_stage2$pca$x[, 1:50, drop = FALSE]
colnames(ubiq_pc_top50) <- paste0("ubiq_PC_", 1:50)
meta_clean_ubiq <- cbind(meta_clean_ubiq, ubiq_pc_top50)

cat("\n🏆 UBIQ RAW+PC (Stage 2) - Top 50 PCs extracted\n")
cat("Ubiq outlier check - samples with extreme PC1 (>3SD):\n")
pc1_sd <- sd(ubiq_pc_top50[,1])
pc1_outliers <- rownames(ubiq_pc_top50)[abs(ubiq_pc_top50[,1] - mean(ubiq_pc_top50[,1])) > 3*pc1_sd]
if(length(pc1_outliers) > 0) {
  cat("Outliers:", paste(pc1_outliers, collapse = ", "), "\n")
} else {
  cat("No extreme outliers detected\n")
}

print("Ubiq Key outputs:")
print("- ubiq_pc_top50: Best stage PCA scores")
print("- ubiq_raw_pc: Best corrected ubiq matrix")
print("- meta_clean_ubiq: Ubiq metadata w/ clusters")
print("- ubiq_metrics: Stage comparison table")

write.csv(ubiq_pc_top50, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/LPM_ubiq_pc_top50.csv", row.names = TRUE)
```


## 9 - data for the LPM - PTMs:Glycosylation 6 (not ready yet)
```{r}

```


## 10 - data for the LPM - Metabolomics
```{r}

Metabolomics_full <- read.csv("/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_30_12_17_25/Metabolomics_PRAD_unannotated-features_All-3Cs_raw-no-normalization(Compounds) (1).csv")
Metabolomics_full <- Metabolomics_full[,-c(2:30)]
rownames(Metabolomics_full) <- Metabolomics_full$Name
Metabolomics_full$remove <- ifelse(grepl("BB_", Metabolomics_full$Name), "remove", "keep")
Metabolomics_full$remove <- ifelse(grepl("ZA_", Metabolomics_full$Name), "remove", Metabolomics_full$remove)
Metabolomics_full$remove <- ifelse(grepl("F5_", Metabolomics_full$Name), "remove", Metabolomics_full$remove)
Metabolomics_full <- dplyr::filter(Metabolomics_full, remove == "keep")
Metabolomics_full$remove <- NULL
Metabolomics_full$Name <- NULL
Metabolomics_full <- Metabolomics_full[-c(1),]

metadata_cptac_prad <- read_xlsx(
  "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_30_12_17_25/SampleList_PRAD sample manifest for metabolomics_for-CPTAC_20250408.xlsx"
)

metadata_cptac_prad$`LC-MS IDs` <- as.character(metadata_cptac_prad$`LC-MS IDs`)
metadata_cptac_prad$ID <- paste0("S_0", metadata_cptac_prad$`LC-MS IDs`)
metadata_cptac_prad <- metadata_cptac_prad[-c(188),]
metadata_cptac_prad$tissue <- ifelse(grepl("_T", metadata_cptac_prad$sampleId), "Tumor", "Normal")
metadata_cptac_prad_t <- dplyr::filter(metadata_cptac_prad, tissue == "Tumor")
metadata_cptac_prad_t$sampleId <- gsub("_T", "", metadata_cptac_prad_t$sampleId)


colnames(Metabolomics_full) <- c(metadata_cptac_prad$sampleId)

Metabolomics_t <- as.data.frame(t(as.matrix(Metabolomics_full)))
Metabolomics_t$sample <- ifelse(grepl("_T", rownames(Metabolomics_t)), "tumor", "normal")
Metabolomics_t <- dplyr::filter(Metabolomics_t, sample == "tumor")
Metabolomics_t$sample <- NULL
rownames(Metabolomics_t) <- gsub("_T", "", rownames(Metabolomics_t))
rownames(Metabolomics_t) <- gsub("\\.", "-", rownames(Metabolomics_t))


# After your existing setup:
Metabolomics_t <- cbind(Metabolomics_t, metadata_cptac_prad_t[, c("aliquoted mass for metabolomics (mg)")])

# Extract weight column FIRST as numeric vector
weight_col <- as.numeric(metadata_cptac_prad_t$`aliquoted mass for metabolomics (mg)`)

cat("Weight normalization - min:", min(weight_col), "max:", max(weight_col), "mean:", mean(weight_col), "\n")
cat("Any zero/negative weights?", any(weight_col <= 0), "\n")
cat("Weight vector length:", length(weight_col), "Data cols:", ncol(Metabolomics_t)-1, "\n")

# 1. WEIGHT NORMALIZATION - Extract metabolite data FIRST (exclude weight column)
cat("\n=== Step 1: Weight normalization ===\n")

# Ensure numeric matrix before division (FIXED: column-wise numeric conversion)
metabo_data <- Metabolomics_t[, 1:(ncol(Metabolomics_t)-1)]

# Convert each column to numeric individually (avoids make.names error)
for (j in seq_len(ncol(metabo_data))) {
  metabo_data[[j]] <- as.numeric(metabo_data[[j]])
}

metabo_mat <- as.matrix(metabo_data)

Metabolomics_weight_norm <- sweep(metabo_mat, 1, weight_col, "/")

# Check for any infinite/NA values after division
n_na <- sum(is.na(Metabolomics_weight_norm))
n_inf <- sum(!is.finite(Metabolomics_weight_norm))
cat("Post-weight norm - NAs:", n_na, "Infs:", n_inf, "\n")

# Replace Infs with NA, then median impute
Metabolomics_weight_norm[!is.finite(Metabolomics_weight_norm)] <- NA
for(i in 1:ncol(Metabolomics_weight_norm)) {
  nas <- is.na(Metabolomics_weight_norm[,i])
  if(any(nas)) {
    med_val <- median(Metabolomics_weight_norm[,i], na.rm = TRUE)
    Metabolomics_weight_norm[nas, i] <- med_val
  }
}

# Convert back to data.frame for downstream steps
Metabolomics_weight_norm <- as.data.frame(Metabolomics_weight_norm)
rownames(Metabolomics_weight_norm) <- rownames(Metabolomics_t)

# 2. CHECK IF LOG TRANSFORM NEEDED
cat("\n=== Step 2: Skewness check for log transform ===\n")
skewness_raw <- apply(Metabolomics_weight_norm, 2, function(x) {
  m3 <- mean((x - mean(x))^3, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if(s == 0 || is.na(s)) return(0)
  m3 / s^3
})
cat("Median skewness (weight-normalized):", median(skewness_raw, na.rm = TRUE), "\n")
cat("Highly skewed metabolites (>1):", sum(skewness_raw > 1, na.rm = TRUE), "/", ncol(Metabolomics_weight_norm), "\n")

# Log2 transform if skewed (typical for metabolomics)
log_transform_needed <- median(skewness_raw, na.rm = TRUE) > 0.5
if(log_transform_needed) {
  cat("Applying log2 transformation...\n")
  Metabolomics_log <- log2(Metabolomics_weight_norm + 1e-9)  # Small pseudo-count
} else {
  cat("Data sufficiently symmetric, skipping log transform\n")
  Metabolomics_log <- Metabolomics_weight_norm
}

# 3. MEDIAN CENTERING (by metabolite/column)
cat("\n=== Step 3: Median centering ===\n")
Metabolomics_final <- apply(Metabolomics_log, 2, function(x) x - median(x, na.rm = TRUE))

# Convert back to data.frame and add metadata
Metabolomics_final <- as.data.frame(Metabolomics_final)
Metabolomics_final$weight_mg <- weight_col
rownames(Metabolomics_final) <- rownames(Metabolomics_t)

cat("\n=== FINAL SUMMARY ===\n")
cat("Final data shape:", nrow(Metabolomics_final), "samples x", ncol(Metabolomics_final), "metabolites + weight\n")
cat("Ready for PCA: Metabolomics_final[, 1:(ncol(Metabolomics_final)-1)]\n")

# Quick PCA check (FIXED)
metabo_pca_ready <- Metabolomics_final[, 1:(ncol(Metabolomics_final)-1)]
n_cols_pca <- min(50, ncol(metabo_pca_ready))
pca_check <- prcomp(t(metabo_pca_ready[,1:n_cols_pca]), scale. = FALSE)
cat("PC1+PC2 variance explained:", round(sum(summary(pca_check)$importance[2,1:2])*100, 1), "%\n")

cat("\n✅ Normalization complete! Ready for PCA/clustering.\n")
cat("Key outputs:\n")
cat("- Metabolomics_final: Weight-norm → Log2 → Median-centered\n")
cat("- metabo_pca_ready: PCA-ready matrix (samples x metabolites)\n")
cat("- weight_col: Original weights for reference\n")

############################################
## Metabolomics Analysis Libraries
############################################

library(dplyr)
library(limma)
library(ggplot2)
library(ggpubr)
library(cluster)

############################################
## 1. Setup & SAFE PCA FUNCTION (metabolomics-specific)
############################################

# Align metabolomics data with metadata
metabo_pca_ready <- as.data.frame(metabo_pca_ready)
meta_clean_metabo <- meta %>% filter(case_id %in% rownames(metabo_pca_ready))
metabo_aligned <- metabo_pca_ready[meta_clean_metabo$case_id, , drop = FALSE]
purity_metabo <- as.numeric(meta_clean_metabo$Purity)
names(purity_metabo) <- meta_clean_metabo$case_id

cat("Metabolomics Dataset:", nrow(metabo_aligned), "samples x", ncol(metabo_aligned), "metabolites\n")

# Check NA percentage (should be 0 after your normalization)
cat("NA percentage:", round(mean(is.na(metabo_aligned))*100, 2), "%\n")

clean_for_pca <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if(ncol(mat) == 0 || nrow(mat) == 0) stop("Empty matrix")
  
  # Remove constants (unlikely after normalization)
  vars <- apply(mat, 2, var, na.rm = TRUE)
  keep_cols <- vars > 0 & is.finite(vars)
  mat <- mat[, keep_cols, drop = FALSE]
  
  if(ncol(mat) == 0) stop("No variable columns")
  
  # Median imputation (backup only)
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  # Final cleanup
  mat[!is.finite(mat)] <- NA
  for(i in 1:ncol(mat)) {
    nas <- is.na(mat[,i])
    if(any(nas)) mat[nas, i] <- median(mat[,i], na.rm = TRUE)
  }
  
  return(mat)
}

n_top_metabo <- min(5000, ncol(metabo_aligned))

do_pca_analysis <- function(prot_mat, stage_name, meta_data) {
  prot_clean <- clean_for_pca(prot_mat)
  var_order <- apply(prot_clean, 2, var, na.rm = TRUE)
  top_prot <- names(sort(var_order, decreasing = TRUE))[1:n_top_metabo]
  
  prot_top <- prot_clean[, top_prot, drop = FALSE]
  z_prot <- scale(prot_top, center = TRUE, scale = TRUE)
  pca_obj <- prcomp(z_prot, center = FALSE, scale. = FALSE)
  
  pc_df <- data.frame(pca_obj$x[, 1:2], case_id = rownames(pca_obj$x)) %>%
    left_join(meta_data[, c("case_id", "ETS", "Purity")], by = "case_id")
  
  km <- kmeans(pca_obj$x[, 1:min(10, ncol(pca_obj$x))], centers = 3, nstart = 50)
  meta_data[[paste0("metabo_cluster_", stage_name)]] <- factor(km$cluster)
  
  list(pca = pca_obj, pc_df = pc_df, meta = meta_data)
}

set.seed(123)

############################################
## 2. STAGE 1: RAW metabolomics (already normalized)
############################################

cat("\n=== METABOLOMICS STAGE 1: Raw (weight+log+median normalized) ===\n")
metabo_stage1 <- do_pca_analysis(metabo_aligned, "raw", meta_clean_metabo)
metabo_pc1_df <- metabo_stage1$pc_df
meta_clean_metabo <- metabo_stage1$meta

############################################
## 3. STAGE 2: RAW + Purity correction (BEST expected)
############################################

cat("\n=== METABOLOMICS STAGE 2: Raw + Purity correction ===\n")
metabo_raw_t <- t(clean_for_pca(metabo_aligned))
metabo_raw_pc_t <- removeBatchEffect(metabo_raw_t, covariates = purity_metabo)
metabo_raw_pc <- t(metabo_raw_pc_t)

metabo_stage2 <- do_pca_analysis(metabo_raw_pc, "raw_pc", meta_clean_metabo)
metabo_pc2_df <- metabo_stage2$pc_df
meta_clean_metabo <- metabo_stage2$meta

############################################
## 4. STAGE 3: Median-imputed (for completeness only)
############################################

cat("\n=== METABOLOMICS STAGE 3: Median-imputed ===\n")
metabo_imputed <- clean_for_pca(metabo_aligned)

metabo_stage3 <- do_pca_analysis(metabo_imputed, "imp", meta_clean_metabo)
metabo_pc3_df <- metabo_stage3$pc_df
meta_clean_metabo <- metabo_stage3$meta

############################################
## 5. STAGE 4: Median-imputed + Purity corrected (FINAL)
############################################

cat("\n=== METABOLOMICS STAGE 4: Median-imputed + Purity corrected ===\n")
metabo_imp_t <- t(clean_for_pca(metabo_imputed))
metabo_imp_pc_t <- removeBatchEffect(metabo_imp_t, covariates = purity_metabo)
metabo_final <- t(metabo_imp_pc_t)

metabo_stage4 <- do_pca_analysis(metabo_final, "final", meta_clean_metabo)
metabo_pc4_df <- metabo_stage4$pc_df
meta_clean_metabo <- metabo_stage4$meta

cat("Metabolomics Pipeline complete!\n")

############################################
## 6. 4-STAGE METABOLOMICS PCA PLOTS
############################################

metabo_p1 <- ggplot(metabo_pc1_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Metabo 1. RAW(norm)", x = "PC1", y = "PC2")

metabo_p2 <- ggplot(metabo_pc2_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Metabo 2. RAW+PC", x = "PC1", y = "PC2")

metabo_p3 <- ggplot(metabo_pc3_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Metabo 3. Median-imp", x = "PC1", y = "PC2")

metabo_p4 <- ggplot(metabo_pc4_df, aes(PC1, PC2, color = ETS)) +
  geom_point(size = 2.5, alpha = 0.8) + theme_bw() +
  labs(title = "Metabo 4. FINAL", x = "PC1", y = "PC2")

ggarrange(metabo_p1, metabo_p2, metabo_p3, metabo_p4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 7. 4-STAGE METABOLOMICS PURITY BOXPLOTS
############################################

metabo_box1 <- ggboxplot(meta_clean_metabo, x = "metabo_cluster_raw", y = "Purity", 
                        color = "metabo_cluster_raw", add = "jitter") + theme_bw() +
  labs(title = "Metabo 1. RAW(norm)")

metabo_box2 <- ggboxplot(meta_clean_metabo, x = "metabo_cluster_raw_pc", y = "Purity", 
                        color = "metabo_cluster_raw_pc", add = "jitter") + theme_bw() +
  labs(title = "Metabo 2. RAW+PC")

metabo_box3 <- ggboxplot(meta_clean_metabo, x = "metabo_cluster_imp", y = "Purity", 
                        color = "metabo_cluster_imp", add = "jitter") + theme_bw() +
  labs(title = "Metabo 3. Median-imp")

metabo_box4 <- ggboxplot(meta_clean_metabo, x = "metabo_cluster_final", y = "Purity", 
                        color = "metabo_cluster_final", add = "jitter") + theme_bw() +
  labs(title = "Metabo 4. FINAL")

ggarrange(metabo_box1, metabo_box2, metabo_box3, metabo_box4, ncol = 2, nrow = 2, common.legend = TRUE)

############################################
## 8. METABOLOMICS METRICS & SUMMARY
############################################

evaluate_stage <- function(pca_obj, meta, cluster_col, purity_col = "Purity") {
  pc1_var <- summary(pca_obj)$importance[2,1]*100
  pc2_var <- summary(pca_obj)$importance[2,2]*100
  
  sil_ets <- silhouette(as.numeric(as.factor(meta$ETS)), dist(pca_obj$x[,1:2]))
  ets_sil <- mean(sil_ets[,3])
  
  purity_p <- summary(aov(meta[[purity_col]] ~ meta[[cluster_col]]))[[1]]$`Pr(>F)`[1]
  
  data.frame(PC1_var = pc1_var, PC2_var = pc2_var, ETS_silhouette = ets_sil, Purity_pval = purity_p)
}

metabo_metrics <- rbind(
  evaluate_stage(metabo_stage1$pca, meta_clean_metabo, "metabo_cluster_raw"),
  evaluate_stage(metabo_stage2$pca, meta_clean_metabo, "metabo_cluster_raw_pc"),
  evaluate_stage(metabo_stage3$pca, meta_clean_metabo, "metabo_cluster_imp"),
  evaluate_stage(metabo_stage4$pca, meta_clean_metabo, "metabo_cluster_final")
)
rownames(metabo_metrics) <- c("METABO_RAW(norm)", "METABO_RAW+PC", "METABO_MEDIAN", "METABO_FINAL")
print("Metabolomics Metrics:")
print(metabo_metrics)

cat("\nMetabolomics PC1 variances (%):")
cat(sprintf(" RAW: %.1f, RAW+PC: %.1f, MEDIAN: %.1f, FINAL: %.1f\n",
    summary(metabo_stage1$pca)$importance[2,1]*100,
    summary(metabo_stage2$pca)$importance[2,1]*100, 
    summary(metabo_stage3$pca)$importance[2,1]*100,
    summary(metabo_stage4$pca)$importance[2,1]*100))

############################################
## 9. EXTRACT BEST METABOLOMICS STAGE (RAW+PC) - TOP 50 PCs
############################################

metabo_pc_top50 <- metabo_stage2$pca$x[, 1:50, drop = FALSE]
colnames(metabo_pc_top50) <- paste0("metabo_PC_", 1:50)
meta_clean_metabo <- cbind(meta_clean_metabo, metabo_pc_top50)

cat("\n🏆 METABOLOMICS RAW+PC (Stage 2) - Top 50 PCs extracted\n")
cat("Metabolomics outlier check - samples with extreme PC1 (>3SD):\n")
pc1_sd <- sd(metabo_pc_top50[,1])
pc1_outliers <- rownames(metabo_pc_top50)[abs(metabo_pc_top50[,1] - mean(metabo_pc_top50[,1])) > 3*pc1_sd]
if(length(pc1_outliers) > 0) {
  cat("Outliers:", paste(pc1_outliers, collapse = ", "), "\n")
} else {
  cat("No extreme outliers detected\n")
}

print("Metabolomics Key outputs:")
print("- metabo_pc_top50: Best stage PCA scores")
print("- metabo_raw_pc: Best corrected metabolomics matrix")
print("- meta_clean_metabo: Metabolomics metadata w/ clusters")
print("- metabo_metrics: Stage comparison table")

# Save outputs
write.csv(metabo_pc_top50, "/mctp/share/users/gondal/04_CPTAC/02_Prostate/03_output/02_CPTAC/version_31_12_29_25/LPM_metabo_pc_top50.csv", row.names = TRUE)


```
