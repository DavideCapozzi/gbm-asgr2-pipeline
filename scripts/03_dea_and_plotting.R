# scripts/03_dea_and_plotting.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(jsonlite)

# 0. Ensure Output Directories Exist
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# 1. Load Clustered Data
gbm <- readRDS("data/processed/02_gbm_clustered.rds")

# 2. Dynamically Identify the Target Cluster (Max ASGR2+ proportion)
# Avoid cluster hijacking by background noise: select by max percentage, 
# but require a minimum absolute count to prevent small-cluster artifacts.
asgr2_expr <- GetAssayData(gbm, layer = "data")["ASGR2", ]
asgr2_pos_cells <- names(asgr2_expr[asgr2_expr > 0])

if (length(asgr2_pos_cells) < 10) {
  stop("Not enough ASGR2+ cells across the entire dataset to perform robust DEA.")
}

# Calculate percentages robustly
cluster_totals <- table(Idents(gbm))
asgr2_cluster_counts <- table(Idents(gbm)[asgr2_pos_cells])

# Align tables (handle clusters with 0 ASGR2+ cells)
aligned_counts <- as.numeric(asgr2_cluster_counts[names(cluster_totals)])
aligned_counts[is.na(aligned_counts)] <- 0

# Calculate percentage
asgr2_percentages <- (aligned_counts / as.numeric(cluster_totals)) * 100
names(asgr2_percentages) <- names(cluster_totals)

# Filter out clusters with less than 10 ASGR2+ cells to ensure statistical power
valid_clusters <- names(cluster_totals)[aligned_counts >= 10]

if (length(valid_clusters) == 0) {
  stop("No individual cluster has enough ASGR2+ cells (>=10) for robust intra-cluster DEA.")
}

# Select the target cluster based on the maximum percentage among valid clusters
valid_percentages <- asgr2_percentages[valid_clusters]
target_cluster <- names(valid_percentages)[which.max(valid_percentages)]

print(paste("Dynamically identified Cluster", target_cluster, "as the primary ASGR2+ niche based on max percentage."))

# 3. Subset the Target Cluster (Macrophage/Microglia population)
macrophages <- subset(gbm, idents = target_cluster)

# 4. Create Metadata for ASGR2 Status
# Cell is 'Positive' if ASGR2 expression > 0, else 'Negative'
macrophages$ASGR2_status <- ifelse(
  colnames(macrophages) %in% asgr2_pos_cells, 
  "Positive", 
  "Negative"
)

# Set identity to the new status for FindMarkers
Idents(macrophages) <- "ASGR2_status"

# 5. Perform Differential Expression Analysis (DEA)
# Comparing ASGR2 Positive vs Negative within the macrophage cluster
print("Running Wilcoxon Rank Sum test for DEA...")
dea_results <- FindMarkers(
  macrophages, 
  ident.1 = "Positive", 
  ident.2 = "Negative",
  test.use = "wilcox",
  logfc.threshold = 0.5,
  min.pct = 0.25         
)

# Add gene names as a column for easier manipulation
dea_results$gene <- rownames(dea_results)

# Remove ASGR2 to avoid "Double Dipping" effect (it is the grouping variable, not a discovery)
dea_results <- dea_results %>% filter(gene != "ASGR2")

write.csv(dea_results, "results/03_ASGR2_DEA_results.csv", row.names = FALSE)

# Filter for statistically significant genes (Adjusted P-value < 0.05)
sig_genes <- dea_results %>% filter(p_val_adj < 0.05)

# 6. Generate Volcano Plot
p_volcano <- EnhancedVolcano(
  dea_results,
  lab = dea_results$gene,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'ASGR2+ vs ASGR2- Macrophages',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 5.0,
  legendPosition = 'right'
)
ggsave("results/03_volcano_plot.pdf", p_volcano, width = 10, height = 8)

# 7. Perform Gene Ontology (GO) Enrichment Analysis
print("Running Gene Ontology (Biological Process) Analysis...")

# Extract significant upregulated genes in ASGR2+ cells
up_genes <- sig_genes %>% filter(avg_log2FC > 0) %>% pull(gene)

# Convert Gene Symbols to Entrez IDs for clusterProfiler
entrez_ids <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_results_list <- list(
  pipeline_step = "03_dea_and_go",
  target_cluster_analyzed = target_cluster,
  cells_in_target_cluster = ncol(macrophages),
  asgr2_positive_cells_analyzed = sum(macrophages$ASGR2_status == "Positive"),
  total_significant_genes = nrow(sig_genes),
  upregulated_genes = length(up_genes),
  go_terms_enriched = 0
)

if (nrow(entrez_ids) > 0) {
  go_enrich <- enrichGO(
    gene          = entrez_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP", # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  if (!is.null(go_enrich) && nrow(go_enrich@result %>% filter(p.adjust < 0.05)) > 0) {
    p_go <- dotplot(go_enrich, showCategory = 15) + ggtitle("GO Enrichment: ASGR2+ Upregulated")
    ggsave("results/03_go_dotplot.pdf", p_go, width = 10, height = 8)
    write.csv(as.data.frame(go_enrich), "results/03_GO_enrichment_results.csv", row.names = FALSE)
    
    go_results_list$go_terms_enriched <- nrow(go_enrich@result %>% filter(p.adjust < 0.05))
    print("GO Analysis successful. Plots saved.")
  } else {
    print("Warning: No significant GO terms found for the upregulated genes.")
  }
} else {
  print("Warning: Not enough upregulated genes mapped to Entrez IDs to perform GO analysis.")
}

# 8. Export Machine-Friendly Results JSON
write_json(go_results_list, "results/03_dea_summary.json", pretty = TRUE, auto_unbox = TRUE)

# Save the subsetted Seurat object for potential future granular analysis
saveRDS(macrophages, "data/processed/03_macrophages_annotated.rds")

print("Script 03 completed successfully. Pipeline execution finished.")