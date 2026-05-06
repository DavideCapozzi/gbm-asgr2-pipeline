# scripts/02_clustering_annotation.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(jsonlite)

# Default config for standalone execution (GSM3828672 behaviour preserved)
if (!exists("config")) {
  config <- list(
    dataset_id        = "GSM3828672",
    use_harmony       = FALSE,
    harmony_var       = NULL,
    genes_of_interest = c("ASGR2", "CLEC10A"),
    results_dir       = "results/GSM3828672",
    processed_dir     = "data/processed/GSM3828672"
  )
}

# 0. Ensure Output Directories Exist
dir.create(config$results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$processed_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load Preprocessed Data
gbm <- readRDS(file.path(config$processed_dir, "01_gbm_preprocessed.rds"))

# 2. Dimensionality Reduction
gbm <- RunPCA(gbm, npcs = 30, verbose = FALSE)

# 3. Optional Harmony Batch Correction (required for multi-sample 10X datasets)
if (isTRUE(config$use_harmony)) {
  library(harmony)
  gbm <- RunHarmony(gbm, group.by.vars = config$harmony_var, verbose = FALSE)
  reduction_use <- "harmony"
} else {
  reduction_use <- "pca"
}

gbm <- RunUMAP(gbm, dims = 1:20, reduction = reduction_use, verbose = FALSE)
gbm <- FindNeighbors(gbm, dims = 1:20, reduction = reduction_use, verbose = FALSE)
gbm <- FindClusters(gbm, resolution = 0.5, verbose = FALSE)

# 4. Save UMAP Plot
p_umap_clusters <- DimPlot(gbm, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Global UMAP - Unannotated Clusters")
ggsave(file.path(config$results_dir, "02_umap_clusters.pdf"), p_umap_clusters, width = 8, height = 6)

# 5. Define Canonical Markers
base_markers <- c("CD68", "CD163", "AIF1", "CD14", "P2RY12", "ITGA4", "CD36")
genes_of_interest <- config$genes_of_interest
target_markers <- unique(c(base_markers, genes_of_interest))
available_markers <- intersect(target_markers, rownames(gbm))
missing_markers <- setdiff(target_markers, rownames(gbm))

if (length(missing_markers) > 0) {
  warning(paste("The following markers were not found in the dataset:", paste(missing_markers, collapse = ", ")))
}

# 6. Generate Marker FeaturePlots
if (length(available_markers) > 0) {
  p_markers <- FeaturePlot(gbm, features = available_markers, ncol = 3, pt.size = 0.5, order = TRUE)
  ggsave(file.path(config$results_dir, "02_umap_macrophage_markers.pdf"), p_markers, width = 15, height = 10)
}

# 7. Generate Machine-Friendly DEA Feasibility Metrics (JSON)
cluster_ids <- levels(Idents(gbm))
genes_present <- intersect(genes_of_interest, available_markers)
feasibility_metrics <- list(
  pipeline_step   = "02_clustering",
  dataset_id      = config$dataset_id,
  batch_corrected = isTRUE(config$use_harmony),
  total_clusters  = length(cluster_ids),
  genes_present   = genes_present,
  cluster_stats   = list()
)

for (cluster in cluster_ids) {
  cells_in_cluster <- WhichCells(gbm, idents = cluster)
  cluster_total_cells <- length(cells_in_cluster)

  cluster_data <- list(
    cluster_id  = cluster,
    total_cells = cluster_total_cells
  )

  for (gene in genes_present) {
    gene_expr <- GetAssayData(gbm, layer = "data")[gene, cells_in_cluster]
    gene_pos <- sum(gene_expr > 0)
    cluster_data[[paste0(gene, "_positive_cells")]] <- gene_pos
    cluster_data[[paste0(gene, "_positivity_pct")]] <- round((gene_pos / cluster_total_cells) * 100, 2)
  }

  feasibility_metrics$cluster_stats[[as.character(cluster)]] <- cluster_data
}

write_json(feasibility_metrics, file.path(config$results_dir, "02_dea_feasibility.json"), pretty = TRUE, auto_unbox = TRUE)

# Save the full clustered object to avoid re-running PCA/UMAP
saveRDS(gbm, file.path(config$processed_dir, "02_gbm_clustered.rds"))

print("ACTION REQUIRED: Inspect 02_umap_macrophage_markers.pdf and 02_dea_feasibility.json.")
print(paste("Script 02 completed successfully for dataset:", config$dataset_id))
