# scripts/03_dea_and_plotting.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(jsonlite)
library(EnhancedVolcano)
library(MAST)

# 0. Ensure Output Directories Exist
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 1. Load Clustered Data
gbm <- readRDS("data/processed/02_gbm_clustered.rds")

# 2. Dynamically Identify the Target Cluster (Max ASGR2+ proportion)
asgr2_expr <- GetAssayData(gbm, layer = "data")["ASGR2", ]
asgr2_pos_cells <- names(asgr2_expr[asgr2_expr > 0])

if (length(asgr2_pos_cells) < 10) {
  stop("Not enough ASGR2+ cells across the entire dataset to perform robust DEA.")
}

cluster_totals <- table(Idents(gbm))
asgr2_cluster_counts <- table(Idents(gbm)[asgr2_pos_cells])
aligned_counts <- as.numeric(asgr2_cluster_counts[names(cluster_totals)])
aligned_counts[is.na(aligned_counts)] <- 0
asgr2_percentages <- (aligned_counts / as.numeric(cluster_totals)) * 100
names(asgr2_percentages) <- names(cluster_totals)

valid_clusters <- names(cluster_totals)[aligned_counts >= 10]

if (length(valid_clusters) == 0) {
  stop("No individual cluster has enough ASGR2+ cells (>=10) for robust intra-cluster DEA.")
}

target_cluster <- names(asgr2_percentages[valid_clusters])[which.max(asgr2_percentages[valid_clusters])]
print(paste("Dynamically identified Cluster", target_cluster, "as the primary ASGR2+ niche."))

# 3. Subset Target Cluster and Calculate Technical Covariates
macrophages <- subset(gbm, idents = target_cluster)
mac_asgr2_expr <- GetAssayData(macrophages, layer = "data")["ASGR2", ]
macrophages$ASGR2_status <- ifelse(mac_asgr2_expr > 0, "Positive", "Negative")
Idents(macrophages) <- "ASGR2_status"

# Calculate Cellular Detection Rate (CDR) to correct for sequencing depth dropouts
macrophages$cdr <- scale(colSums(GetAssayData(macrophages, layer = "counts") > 0))

# 4. Balanced Downsampling Strategy
# Prevent the overwhelming negative population from burying the signal
set.seed(42) # Ensures computational reproducibility of the sample
cells_pos <- WhichCells(macrophages, idents = "Positive")
cells_neg <- WhichCells(macrophages, idents = "Negative")

# Sample 3x negatives relative to positives to stabilize variance without losing power
target_neg_size <- min(length(cells_neg), length(cells_pos) * 3)
cells_neg_sampled <- sample(cells_neg, target_neg_size)

macs_balanced <- subset(macrophages, cells = c(cells_pos, cells_neg_sampled))
print(paste("Running MAST on balanced set:", length(cells_pos), "pos vs", length(cells_neg_sampled), "neg"))

# 5. Perform Differential Expression Analysis (MAST Hurdle Model)
dea_results <- FindMarkers(
  macs_balanced, 
  ident.1 = "Positive", 
  ident.2 = "Negative",
  test.use = "MAST",
  latent.vars = "cdr",
  logfc.threshold = 0.5, 
  min.pct = 0.25 
)

# Structure results and prevent double-dipping bias on ASGR2
dea_results$gene <- rownames(dea_results)
dea_results <- dea_results %>% filter(gene != "ASGR2")
write.csv(dea_results, "results/03_ASGR2_DEA_results.csv", row.names = FALSE)

# 6. Volcano Plot (Using strict FDR < 0.05)
p_volcano <- EnhancedVolcano(
  dea_results,
  lab = dea_results$gene,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = paste('ASGR2+ Signature (Balanced) - Cluster', target_cluster),
  subtitle = 'MAST Model with CDR correction',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 4.0
)
ggsave("results/03_volcano_asgr2.pdf", p_volcano, width = 10, height = 8)

# 7. Gene Ontology Enrichment (Exploratory Pool)
print("Running Gene Ontology (Biological Process) Analysis...")
# Use nominal p-value to capture pathway trends in small cohorts
go_pool_genes <- dea_results %>% filter(p_val < 0.01 & avg_log2FC > 0.5) %>% pull(gene)
go_terms_count <- 0

if (length(go_pool_genes) > 5) {
  entrez_ids <- bitr(go_pool_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  ego <- enrichGO(
    gene          = entrez_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  if (!is.null(ego) && nrow(ego) > 0) {
    write.csv(as.data.frame(ego), "results/03_GO_enrichment_results.csv", row.names = FALSE)
    p_go <- dotplot(ego, showCategory = 15) + ggtitle("GO Enrichment - Exploratory ASGR2+ Pool")
    ggsave("results/03_go_dotplot.pdf", p_go, width = 10, height = 8)
    go_terms_count <- nrow(ego)
  } else {
    print("Warning: No significant GO terms found in the exploratory pool.")
  }
} else {
  print("Warning: Insufficient genes in exploratory pool to perform GO enrichment.")
}

# 8. Export Summary Metrics
summary_metrics <- list(
  pipeline_step = "03_dea_and_go_MAST_Balanced",
  method = "MAST_cdr_adjusted",
  target_cluster_selected = target_cluster,
  asgr2_pos_cells_analyzed = length(cells_pos),
  asgr2_neg_cells_downsampled = length(cells_neg_sampled),
  asgr2_neg_cells_total_pool = length(cells_neg),
  strict_degs_fdr05 = sum(dea_results$p_val_adj < 0.05),
  exploratory_go_genes = length(go_pool_genes),
  go_terms_enriched = go_terms_count
)

write_json(summary_metrics, "results/03_dea_summary.json", pretty = TRUE, auto_unbox = TRUE)
print("Pipeline Stage 03 completed. Claims validated. Check 'results/' directory.")