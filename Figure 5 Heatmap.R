#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(readr)
  library(tibble)
  library(matrixStats)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(grid)
})

# ── [1] Settings ────────────────────────────────────────────────────────
WD       <- "" #add your working directory here
groups   <- c("Immune", "Neuronal_Support")
PC_tag   <- "PreClinical_vs_Control"
CL_tag   <- "Clinical_vs_Control"
bio_vars <- c("age","sex","APOE","AB42/AB40 RATIO","UPSIT")

# ── [2] Compute best logFC per gene across both comparisons ────────────
best_fc <- numeric()
update_best <- function(df) {
  idx <- which(df$logFC > 0 & df$FDR < 0.05)
  for (i in idx) {
    g <- rownames(df)[i]; v <- df$logFC[i]
    if (!g %in% names(best_fc) || v > best_fc[g]) best_fc[g] <<- v
  }
}
for (g in groups) {
  pc <- read.csv(file.path(WD, sprintf("DE_%s_%s.csv", g, PC_tag)),
                 row.names=1, check.names=FALSE)
  cl <- read.csv(file.path(WD, sprintf("DE_%s_%s.csv", g, CL_tag)),
                 row.names=1, check.names=FALSE)
  update_best(pc)
  update_best(cl)
}
keep_names <- names(best_fc)[
  !grepl("^(MT-|H2|H3)", names(best_fc)) &
    !grepl("\\.", names(best_fc)) &
    !grepl("ENSG", names(best_fc))
]

# ── [3] Identify candidate genes up in PC or CL per group ──────────────
up_by_group <- setNames(vector("list", length(groups)), groups)
for (g in groups) {
  df_pc <- read.csv(file.path(WD, sprintf("DE_%s_%s.csv", g, PC_tag)),
                    row.names=1, check.names=FALSE)
  df_cl <- read.csv(file.path(WD, sprintf("DE_%s_%s.csv", g, CL_tag)),
                    row.names=1, check.names=FALSE)
  genes_pc <- rownames(df_pc)[df_pc$logFC > 0 & df_pc$FDR < 0.05]
  genes_cl <- rownames(df_cl)[df_cl$logFC > 0 & df_cl$FDR < 0.05]
  up_by_group[[g]] <- union(genes_pc, genes_cl)
}
candidate_genes <- intersect(unique(unlist(up_by_group)), keep_names)

# ── [4] Load pseudobulk & metadata, restrict to candidate_genes ────────
expr_l <- list(); meta_l <- list()
for (g in groups) {
  counts <- as.matrix(read.csv(
    file.path(WD, sprintf("pseudobulk_counts_%s.csv", g)),
    row.names=1, check.names=FALSE))
  meta   <- read.csv(
    file.path(WD, sprintf("pseudobulk_meta_%s.csv", g)),
    row.names=1, check.names=FALSE)[colnames(counts), ]
  log2cpm <- edgeR::cpm(counts, log=TRUE, prior.count=0.5)
  expr_l[[g]] <- log2cpm[
    rownames(log2cpm) %in% candidate_genes, , drop=FALSE
  ]
  meta_l[[g]] <- meta %>% transmute(
    bulk_id    = rownames(meta),
    orig_patient,
    cell_group = g,
    Alz_status = factor(gsub("-", "", Alz_status),
                        levels=c("Control","PreClinical","Clinical")),
    across(all_of(bio_vars), identity)
  )
}
expr_mat <- do.call(cbind, expr_l)
meta_df  <- bind_rows(meta_l) %>%
  remove_rownames() %>%
  column_to_rownames("bulk_id") %>%
  .[colnames(expr_mat), ]

# ── [5] Average replicates per patient × cell_group ────────────────────
pat_cols   <- interaction(meta_df$orig_patient,
                          meta_df$cell_group, sep=".", drop=TRUE)
patient_mat <- sapply(split(seq_len(ncol(expr_mat)), pat_cols),
                      function(ix) rowMeans(expr_mat[, ix, drop=FALSE]))
colnames(patient_mat) <- names(split(seq_len(ncol(expr_mat)), pat_cols))
meta_pat <- meta_df %>%
  mutate(col_id=pat_cols) %>% group_by(col_id) %>% slice(1) %>% ungroup() %>%
  column_to_rownames("col_id") %>%
  .[colnames(patient_mat), ]

# ── [6] Compute per-gene Z-scores ──────────────────────────────────────
gene_mat <- t(scale(t(patient_mat)))

# ── [7] Filter to genes positive in all diseased samples for their DE group ─
diseased_status <- c("PreClinical","Clinical")
keep_genes <- Filter(function(gene) {
  gs   <- groups[sapply(groups, function(gr) gene %in% up_by_group[[gr]])]
  cols <- which(meta_pat$Alz_status %in% diseased_status &
                  meta_pat$cell_group %in% gs)
  length(cols)>0 && all(gene_mat[gene, cols] > 0)
}, rownames(gene_mat))
gene_mat <- gene_mat[keep_genes, , drop=FALSE]

# ── [8] Keep top 20 per group by best_fc ────────────────────────────────
immune_cand    <- keep_genes[keep_genes %in% up_by_group[["Immune"]]]
neuronal_cand  <- keep_genes[keep_genes %in% up_by_group[["Neuronal_Support"]]]
top20_immune   <- head(names(sort(best_fc[immune_cand],   decreasing=TRUE)), 20)
top20_neuronal <- head(names(sort(best_fc[neuronal_cand], decreasing=TRUE)), 20)
keep20         <- union(top20_immune, top20_neuronal)
gene_mat <- gene_mat[keep20, , drop=FALSE]

# ── [9] Re-cluster & split rows into k = 2 clusters ────────────────────
row_d    <- dist(gene_mat)
row_hc   <- hclust(row_d, method="ward.D2")
clusters <- cutree(row_hc, k=2)

# ── [10] Build annotations ──────────────────────────────────────────────
meta_pat <- meta_pat %>% rename(ABratio=`AB42/AB40 RATIO`)
sex_pal  <- c(`1`="#FB6A4A", `0`="#3182BD")
apoe_pal <- c(`2/3`="#8DD3C7", `3/3`="#80B1D3",
              `3/4`="#FDB462", `4/4`="#FB8072")
bottom_ha <- HeatmapAnnotation(
  age     = anno_barplot(meta_pat$age,     gp=gpar(fill="#9ecae1")),
  ABratio = anno_barplot(meta_pat$ABratio, gp=gpar(fill="#9ecae1")),
  UPSIT   = anno_barplot(meta_pat$UPSIT,   gp=gpar(fill="#9ecae1")),
  sex     = anno_simple(meta_pat$sex,      col=sex_pal),
  APOE    = anno_simple(meta_pat$APOE,     col=apoe_pal),
  which   = "column", annotation_name_side=NULL, gap=unit(2,"mm")
)
status_cols <- c(Control="#a6cee3", PreClinical="#1f78b4", Clinical="#08519c")
group_cols  <- c(Immune="#e41a1c", Neuronal_Support="#377eb8")
top_ha <- HeatmapAnnotation(
  Status   = meta_pat$Alz_status,
  CellType = meta_pat$cell_group,
  col      = list(Status=status_cols, CellType=group_cols),
  gp       = gpar(col=NA),
  show_annotation_name=FALSE
)
col_fun <- colorRamp2(c(-2,0,2), c("blue","white","red"))

# ── [11] Draw the heatmap ───────────────────────────────────────────────
ht <- Heatmap(
  gene_mat,
  name            = "Z-score",
  cluster_rows    = FALSE,
  row_order       = row_hc$order,
  row_split       = clusters,
  cluster_columns = FALSE,
  show_row_names    = FALSE,
  show_column_names = FALSE,
  column_split      = meta_pat$Alz_status,
  top_annotation    = top_ha,
  bottom_annotation = bottom_ha,
  col               = col_fun
)
draw(ht,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     merge_legends          = TRUE)

###z-score box and whisker plots
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# `gene_mat`, `clusters`, and `meta_pat` were created above:

# 1) Genes in each cluster
cluster_assignments <- clusters
genes_c1 <- names(cluster_assignments)[cluster_assignments == 1]
genes_c2 <- names(cluster_assignments)[cluster_assignments == 2]

# 2) Plot for Cluster 1 among Immune cells
cols_imm <- which(meta_pat$cell_group == "Immune")
expr_c1_imm <- gene_mat[genes_c1, cols_imm, drop=FALSE]
df1 <- data.frame(
  Expression = colMeans(expr_c1_imm),
  AlzStatus  = meta_pat$Alz_status[cols_imm]
)
p1 <- ggplot(df1, aes(x = AlzStatus, y = Expression)) +
  geom_boxplot() +
  labs(title = "Cluster 1 (Immune) Δ-Z by Alz Status") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# 3) Plot for Cluster 2 among Neuronal_Support cells
cols_neu <- which(meta_pat$cell_group == "Neuronal_Support")
expr_c2_neu <- gene_mat[genes_c2, cols_neu, drop=FALSE]
df2 <- data.frame(
  Expression = colMeans(expr_c2_neu),
  AlzStatus  = meta_pat$Alz_status[cols_neu]
)
p2 <- ggplot(df2, aes(x = AlzStatus, y = Expression)) +
  geom_boxplot() +
  labs(title = "Cluster 2 (Neuronal) Δ-Z by Alz Status") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display both plots
print(p1)
print(p2)

## ──  Export per‑module gene list with Z‑scores ─────────────────────
library(readr)          

# 1) pick the columns belonging to each module
cols_imm <- which(meta_pat$cell_group == "Immune")
cols_neu <- which(meta_pat$cell_group == "Neuronal_Support")

# 2) take the mean Z‑score for every gene within its module
mean_z_imm <- rowMeans(gene_mat[genes_c1, cols_imm, drop = FALSE])
mean_z_neu <- rowMeans(gene_mat[genes_c2, cols_neu, drop = FALSE])

# 3) assemble the long‑form table
out_tbl <- rbind(
  data.frame(gene = names(mean_z_imm),
             module = "Immune",
             mean_Z = mean_z_imm,
             stringsAsFactors = FALSE),
  data.frame(gene = names(mean_z_neu),
             module = "Neuronal_Support",
             mean_Z = mean_z_neu,
             stringsAsFactors = FALSE)
)

# 4) write file
write_csv(out_tbl,
          file = "~/Downloads/immune_neuronal_module_genes_zscores.csv")

message("CSV written to ~/Downloads/immune_neuronal_module_genes_zscores.csv")

