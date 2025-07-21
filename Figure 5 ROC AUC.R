#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
})

# - First, create colored deltaZ score box and whisker
# — assume `gene_mat`, `clusters`, and `meta_pat` are already defined —

# Prepare data for Cluster 1 (Immune)
genes_c1 <- names(clusters)[clusters == 1]
cols_imm <- which(meta_pat$cell_group == "Immune")
df1 <- data.frame(
  Expression = colMeans(gene_mat[genes_c1, cols_imm, drop = FALSE]),
  AlzStatus  = meta_pat$Alz_status[cols_imm]
)

# Prepare data for Cluster 2 (Neuronal_Support)
genes_c2 <- names(clusters)[clusters == 2]
cols_neu <- which(meta_pat$cell_group == "Neuronal_Support")
df2 <- data.frame(
  Expression = colMeans(gene_mat[genes_c2, cols_neu, drop = FALSE]),
  AlzStatus  = meta_pat$Alz_status[cols_neu]
)

# Shared y-axis breaks and limits
y_limits <- c(0, 1.5)
y_breaks <- seq(0, 1.5, 0.5)

# Colors and comparisons
status_cols <- c(
  Control     = "#b3c8e7",
  PreClinical = "#cabcdc",
  Clinical    = "#e5c6dd"
)
my_comparisons <- list(
  c("Control", "PreClinical"),
  c("Control", "Clinical"),
  c("PreClinical", "Clinical")
)

# Plot Cluster 1 with fixed y-axis
p1 <- ggplot(df1, aes(x = AlzStatus, y = Expression, fill = AlzStatus)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "gray20") +
  scale_fill_manual(values = status_cols) +
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +
  stat_compare_means(
    comparisons = my_comparisons,
    method      = "wilcox.test",
    label       = "p.signif",
    hide.ns     = TRUE
  ) +
  labs(title = "Cluster 1 (Immune) Δ-Z by Alz Status") +
  theme_minimal() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title      = element_text(size = 12, face = "bold")
  )

# Plot Cluster 2 with fixed y-axis
p2 <- ggplot(df2, aes(x = AlzStatus, y = Expression, fill = AlzStatus)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "gray20") +
  scale_fill_manual(values = status_cols) +
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +
  stat_compare_means(
    comparisons = my_comparisons,
    method      = "wilcox.test",
    label       = "p.signif",
    hide.ns     = TRUE
  ) +
  labs(title = "Cluster 2 (Neuronal) Δ-Z by Alz Status") +
  theme_minimal() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title      = element_text(size = 12, face = "bold")
  )

# Display side by side
grid.arrange(p1, p2, ncol = 2)

# Function to run and print pairwise Wilcoxon tests
print_pairwise_stats <- function(df, cluster_name) {
  cat("\n=== Pairwise Wilcoxon tests for", cluster_name, "===\n")
  for(comp in my_comparisons) {
    x <- df$Expression[df$AlzStatus == comp[1]]
    y <- df$Expression[df$AlzStatus == comp[2]]
    res <- wilcox.test(x, y, exact = FALSE)
    cat(sprintf("%s vs %s: W=%.1f, p=%.3g\n",
                comp[1], comp[2],
                res$statistic,
                res$p.value))
  }
}

# Print stats for Cluster 1 (Immune)
print_pairwise_stats(df1, "Cluster 1 (Immune)")

# Print stats for Cluster 2 (Neuronal) 
print_pairwise_stats(df2, "Cluster 2 (Neuronal)")


## Create ROC AUC plot 
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
  library(ggplot2)
})

# — assume you already have `patient_mat` (genes × samples),
#   `clusters` (gene → 1 or 2), and `meta_pat` (rows=samples, with Alz_status)

# 1) Compute module scores: mean Δ-Z per module per patient
immune_genes   <- names(clusters)[clusters == 1]
neuronal_genes <- names(clusters)[clusters == 2]

scores <- data.frame(
  Sample      = colnames(patient_mat),
  ImmuneScore   = colMeans(gene_mat[immune_genes, , drop=FALSE]),
  NeuronalScore = colMeans(gene_mat[neuronal_genes, , drop=FALSE])
) %>%
  mutate(
    ComboScore = ImmuneScore + NeuronalScore,
    Status     = meta_pat$Alz_status
  )

# 2) Binary outcomes: PreClinical vs Control, Clinical vs Control, etc.
#    Here we'll do PreClinical+Clinical vs Control
scores <- scores %>%
  mutate(
    DzBin = ifelse(Status %in% c("PreClinical","Clinical"), 1, 0)
  )

# 3) Compute ROC curves
roc_imm   <- roc(scores$DzBin, scores$ImmuneScore, quiet=TRUE)
roc_neur  <- roc(scores$DzBin, scores$NeuronalScore, quiet=TRUE)
roc_combo <- roc(scores$DzBin, scores$ComboScore, quiet=TRUE)

# 4) Prepare for ggplot
rocs <- list(
  Immune   = roc_imm,
  Neuronal = roc_neur,
  Combo    = roc_combo
)
df_rocs <- do.call(rbind, lapply(names(rocs), function(name) {
  coords <- coords(rocs[[name]], "all", ret=c("specificity","sensitivity"))
  data.frame(
    Module = name,
    Specificity = coords$specificity,
    Sensitivity = coords$sensitivity
  )
}))

# 5) Plot overlaid ROC curves
ggplot(df_rocs, aes(x = 1 - Specificity, y = Sensitivity, color = Module)) +
  geom_path(size = 1.2) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = c("Immune"="#e41a1c",
                                "Neuronal"="#377eb8",
                                "Combo"="#4daf4a")) +
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    title = "ROC Curves: Immune vs Neuronal vs Combined Modules",
    subtitle = sprintf("AUCs: Immune=%.3f, Neuronal=%.3f, Combo=%.3f",
                       auc(roc_imm), auc(roc_neur), auc(roc_combo))
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title   = element_text(face="bold")
  )

# box and whisker that includes combined module 
## ── 0 · pick the genes and columns for the combined module ───────────
genes_combo <- c(genes_c1, genes_c2)     # immune + neuronal genes
cols_all    <- seq_len(ncol(gene_mat))   # use every column (cells / samples)

df3 <- data.frame(
  Expression = colMeans(gene_mat[genes_combo, cols_all, drop = FALSE]),
  AlzStatus  = meta_pat$Alz_status[cols_all]
)

## ── 1 · reuse colors, breaks, comparisons defined earlier ───────────
p3 <- ggplot(df3, aes(x = AlzStatus, y = Expression, fill = AlzStatus)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, colour = "gray20") +
  scale_fill_manual(values = status_cols) +
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +
  stat_compare_means(
    comparisons = my_comparisons,
    method      = "wilcox.test",
    label       = "p.signif",
    hide.ns     = TRUE
  ) +
  labs(title = "Combined Module Δ-Z by Alz Status") +
  theme_minimal() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title      = element_text(size = 12, face = "bold")
  )

## ── 2 · display the three panels side-by-side ────────────────────────
grid.arrange(p1, p2, p3, ncol = 3)

## ── 3 · print pairwise stats for the combined module ──────
print_pairwise_stats(df3, "Combined Module")


## ── Patient-level module scores & three box-and-whisker plots ─────────

# (Assumes: gene_mat, clusters, meta_pat, status_cols, my_comparisons are in memory)

# 1) module gene lists
immune_genes   <- names(clusters)[clusters == 1]
neuronal_genes <- names(clusters)[clusters == 2]

# 2) sample-level module means
ImmuneExpr   <- colMeans(gene_mat[immune_genes,   , drop = FALSE])
NeuronalExpr <- colMeans(gene_mat[neuronal_genes, , drop = FALSE])

meta_pat <- meta_pat %>%                          # add to metadata (row-aligned)
  mutate(ImmuneExpr   = ImmuneExpr,
         NeuronalExpr = NeuronalExpr)

# 3) collapse to one row per patient (orig_patient)  
patient_df <- meta_pat %>%
  group_by(orig_patient) %>%
  summarise(
    Status        = first(Alz_status),                          # one status per patient
    ImmuneScore   = mean(ImmuneExpr  [cell_group == "Immune"],            na.rm = TRUE),
    NeuronalScore = mean(NeuronalExpr[cell_group == "Neuronal_Support"],  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(ComboScore = ImmuneScore + NeuronalScore) %>%          # combined module
  rename(Patient_ID = orig_patient)

# 4) tidy into plotting frames
df1 <- patient_df %>% transmute(Expression = ImmuneScore,   AlzStatus = Status)
df2 <- patient_df %>% transmute(Expression = NeuronalScore, AlzStatus = Status)
df3 <- patient_df %>% transmute(Expression = ComboScore,    AlzStatus = Status)

# 5) neat half- / whole-unit y-axes per panel
pretty_limits <- function(x, step_small = 0.5, step_big = 1, thresh = 4) {
  span <- diff(range(x, na.rm = TRUE)); step <- ifelse(span > thresh, step_big, step_small)
  lo <- floor(min(x)/step)*step; hi <- ceiling(max(x)/step)*step
  list(lims = c(lo, hi), brks = seq(lo, hi, by = step))
}
ylab1 <- pretty_limits(df1$Expression)
ylab2 <- pretty_limits(df2$Expression)
ylab3 <- pretty_limits(df3$Expression)

make_panel <- function(df, yinfo, ttl) {
  ggplot(df, aes(x = AlzStatus, y = Expression, fill = AlzStatus)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, colour = "gray20") +
    scale_fill_manual(values = status_cols) +
    scale_y_continuous(breaks = yinfo$brks, limits = yinfo$lims) +
    stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       label = "p.signif", hide.ns = TRUE) +
    labs(title = ttl) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(size = 12, face = "bold"))
}

p1 <- make_panel(df1, ylab1, "Immune Module Δ-Z")
p2 <- make_panel(df2, ylab2, "Neuronal-Support Module Δ-Z")
p3 <- make_panel(df3, ylab3, "Combined Module Δ-Z")

gridExtra::grid.arrange(p1, p2, p3, ncol = 3)

###################
# ── Export per‑patient Δ‑Z for the combined module ─────────────────────
suppressPackageStartupMessages({ library(readr) })

out_csv <- "~/Downloads/combined_module_deltaZ_per_patient.csv"

patient_df %>%
  select(Patient_ID,
         AlzStatus = Status,
         ComboScore) %>%          # the value plotted in the last box‑and‑whisker
  arrange(AlzStatus, Patient_ID) %>%
  write_csv(out_csv)

message("CSV written to: ", out_csv)


## bootstrapping, AUC CIs
library(pROC)
## 1 ─ ROC objects + 2 000‑bootstrap 95 % CI
roc_imm   <- roc(scores$DzBin, scores$ImmuneScore,   quiet = TRUE)
ci_imm    <- ci.auc(roc_imm,   boot.n = 2000, boot.stratified = TRUE)

roc_neur  <- roc(scores$DzBin, scores$NeuronalScore, quiet = TRUE)
ci_neur   <- ci.auc(roc_neur,  boot.n = 2000, boot.stratified = TRUE)

roc_combo <- roc(scores$DzBin, scores$ComboScore,    quiet = TRUE)
ci_combo  <- ci.auc(roc_combo, boot.n = 2000, boot.stratified = TRUE)

## 2 ─ store them in lists with matching names
rocs <- list(Immune   = roc_imm,
             Neuronal = roc_neur,
             Combo    = roc_combo)

cis  <- list(Immune   = ci_imm,
             Neuronal = ci_neur,
             Combo    = ci_combo)

## 3 ─ print AUCs and CIs
for(nm in names(rocs)) {
  cat(sprintf("%s: AUC = %.3f  (95%% CI %.3f–%.3f)\n",
              nm, auc(rocs[[nm]]), cis[[nm]][1], cis[[nm]][3]))
}
