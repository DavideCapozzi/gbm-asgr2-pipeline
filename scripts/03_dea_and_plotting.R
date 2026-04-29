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

# 2. Dynamically Identify the Target Cluster (Max ASGR2+ cells)
# Instead of hardcoding '1', we find it programmatically based on expression
asgr2_expr <- GetAssayData(gbm, layer = "data")["ASGR2", ]
asgr2_pos_cells <- names(asgr2_expr[asgr2_expr > 0])

if (length(asgr2_pos_cells) < 10) {
  stop("Not enough ASGR2+ cells across the entire dataset to perform robust DEA.")
}

# Count ASGR2+ cells per cluster
cluster_counts <- table(Idents(gbm)[asgr2_pos_cells])
target_cluster <- names(cluster_counts)[which.max(cluster_counts)]

print(paste("Dynamically identified Cluster", target_cluster, "as the primary ASGR2+ niche."))

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
  logfc.threshold = 0.25, # Standard threshold, kept conservative to capture subtle changes
  min.pct = 0.1
)

# Add gene names as a column for easier manipulation
dea_results$gene <- rownames(dea_results)
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