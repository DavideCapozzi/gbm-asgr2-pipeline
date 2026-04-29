# scripts/02_clustering_annotation.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(jsonlite) 

# 0. Ensure Output Directories Exist
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# 1. Load Preprocessed Data
gbm <- readRDS("data/processed/01_gbm_preprocessed.rds")

# 2. Dimensionality Reduction
gbm <- RunPCA(gbm, npcs = 30, verbose = FALSE)
gbm <- RunUMAP(gbm, dims = 1:20, verbose = FALSE)
gbm <- FindNeighbors(gbm, dims = 1:20, verbose = FALSE)
gbm <- FindClusters(gbm, resolution = 0.5, verbose = FALSE)

# 3. Save Classic Seurat Colored UMAP Plot
p_umap_clusters <- DimPlot(gbm, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Global UMAP - Unannotated Clusters")
ggsave("results/02_umap_clusters.pdf", p_umap_clusters, width = 8, height = 6)

# 4. Define Canonical Markers
target_markers <- c("CD68", "CD163", "AIF1", "CD14", "P2RY12", "ASGR2", "ITGA4")
available_markers <- intersect(target_markers, rownames(gbm))
missing_markers <- setdiff(target_markers, available_markers)

if (length(missing_markers) > 0) {
  warning(paste("The following markers were not found in the dataset:", paste(missing_markers, collapse = ", ")))
}

# 5. Generate Marker FeaturePlots
if (length(available_markers) > 0) {
  p_markers <- FeaturePlot(gbm, features = available_markers, ncol = 3, pt.size = 0.5, order = TRUE)
  ggsave("results/02_umap_macrophage_markers.pdf", p_markers, width = 15, height = 10)
}

# 6. Generate Machine-Friendly DEA Feasibility Metrics (JSON)
# We calculate cells per cluster and check ASGR2 positivity to evaluate if step 3 is statistically viable
cluster_ids <- levels(Idents(gbm))
feasibility_metrics <- list(
  pipeline_step = "02_clustering",
  total_clusters = length(cluster_ids),
  ASGR2_present = "ASGR2" %in% available_markers,
  cluster_stats = list()
)

for (cluster in cluster_ids) {
  cells_in_cluster <- WhichCells(gbm, idents = cluster)
  cluster_total_cells <- length(cells_in_cluster)
  
  cluster_data <- list(
    cluster_id = cluster,
    total_cells = cluster_total_cells
  )
  
  # If ASGR2 exists, calculate how many cells express it (>0) in this cluster
  if (feasibility_metrics$ASGR2_present) {
    asgr2_expr <- GetAssayData(gbm, layer = "data")["ASGR2", cells_in_cluster]
    asgr2_pos_cells <- sum(asgr2_expr > 0)
    
    cluster_data$ASGR2_positive_cells <- asgr2_pos_cells
    cluster_data$ASGR2_positivity_pct <- round((asgr2_pos_cells / cluster_total_cells) * 100, 2)
  }
  
  feasibility_metrics$cluster_stats[[as.character(cluster)]] <- cluster_data
}

write_json(feasibility_metrics, "results/02_dea_feasibility.json", pretty = TRUE, auto_unbox = TRUE)

# Save the full clustered object to avoid re-running PCA/UMAP
saveRDS(gbm, "data/processed/02_gbm_clustered.rds")

print("ACTION REQUIRED: Inspect results/02_umap_macrophage_markers.pdf and results/02_dea_feasibility.json.")
print("Script 02 completed successfully. UMAP and Feasibility metrics exported.")