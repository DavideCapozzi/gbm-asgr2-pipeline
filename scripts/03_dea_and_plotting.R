# scripts/03_dea_and_plotting.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(jsonlite)
library(EnhancedVolcano)
library(MAST)

# --- 0. Config Safety & Environment Inheritance ---
if (!exists("config")) {
  config <- list(
    dataset_id        = "GSE162631",
    genes_of_interest = c("ASGR2"),
    results_dir       = "results/GSE162631",
    processed_dir     = "data/processed/GSE162631"
  )
}

# Ensure Output Directories Exist
dir.create(config$results_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load Clustered Data
file_path_seurat <- file.path(config$processed_dir, "02_gbm_clustered.rds")
if (!file.exists(file_path_seurat)) stop(paste("Cannot find clustered data at:", file_path_seurat))
gbm <- readRDS(file_path_seurat)

# 2. Run DEA pipeline for each gene of interest
for (gene in config$genes_of_interest) {
  print(paste0("\n=== DEA for gene: ", gene, " ==="))
  
  # --- Identify Target Cluster ---
  gene_expr <- GetAssayData(gbm, layer = "data")[gene, ]
  gene_pos_cells <- names(gene_expr[gene_expr > 0])
  
  if (length(gene_pos_cells) < 10) {
    warning(paste0("Not enough ", gene, "+ cells across the dataset (n=", length(gene_pos_cells),
                   "). Skipping DEA for this gene."))
    next
  }
  
  cluster_totals <- table(Idents(gbm))
  gene_cluster_counts <- table(Idents(gbm)[gene_pos_cells])
  aligned_counts <- as.numeric(gene_cluster_counts[names(cluster_totals)])
  aligned_counts[is.na(aligned_counts)] <- 0
  gene_percentages <- (aligned_counts / as.numeric(cluster_totals)) * 100
  names(gene_percentages) <- names(cluster_totals)
  
  valid_clusters <- names(cluster_totals)[aligned_counts >= 10]
  
  if (length(valid_clusters) == 0) {
    warning(paste0("No cluster has >= 10 ", gene, "+ cells. Skipping DEA."))
    next
  }
  
  target_cluster <- names(gene_percentages[valid_clusters])[which.max(gene_percentages[valid_clusters])]
  print(paste("Dynamically identified Cluster", target_cluster, "as the primary", gene, "niche."))
  
  # --- Subset Target Cluster and Calculate Technical Covariates ---
  macrophages <- subset(gbm, idents = target_cluster)
  mac_gene_expr <- GetAssayData(macrophages, layer = "data")[gene, ]
  macrophages$gene_status <- ifelse(mac_gene_expr > 0, "Positive", "Negative")
  Idents(macrophages) <- "gene_status"
  
  macrophages$cdr <- scale(colSums(GetAssayData(macrophages, layer = "counts") > 0))
  
  # --- Balanced Downsampling ---
  set.seed(42)
  cells_pos <- WhichCells(macrophages, idents = "Positive")
  cells_neg <- WhichCells(macrophages, idents = "Negative")
  
  target_neg_size <- min(length(cells_neg), length(cells_pos) * 3)
  cells_neg_sampled <- sample(cells_neg, target_neg_size)
  
  macs_balanced <- subset(macrophages, cells = c(cells_pos, cells_neg_sampled))
  print(paste("Running MAST on balanced set:", length(cells_pos), "pos vs", length(cells_neg_sampled), "neg"))
  
  # --- MAST Hurdle Model ---
  dea_results <- FindMarkers(
    macs_balanced,
    ident.1       = "Positive",
    ident.2       = "Negative",
    test.use      = "MAST",
    latent.vars   = "cdr",
    logfc.threshold = 0.5,
    min.pct       = 0.25
  )
  
  # Structure results; exclude the gene itself to prevent double-dipping
  dea_results$gene_col <- rownames(dea_results)
  dea_results <- dea_results %>% filter(gene_col != gene)
  names(dea_results)[names(dea_results) == "gene_col"] <- "gene"
  
  write.csv(dea_results,
            file.path(config$results_dir, paste0("03_", gene, "_DEA_results.csv")),
            row.names = FALSE)
  
  # ==============================================================================
  # VISUALIZATION SUITE
  # ==============================================================================
  
  # --- 1. Clean Volcano Plot ---
  print("Generating Scientific Volcano Plot...")
  p_volcano <- EnhancedVolcano(
    dea_results,
    lab      = dea_results$gene,
    x        = 'avg_log2FC',
    y        = 'p_val_adj',
    title    = paste0("Differential Expression: ", gene, "+ Macrophages"),
    subtitle = paste0('Dataset: ', config$dataset_id, ' | Hurdle Model (MAST)'),
    caption  = "Cutoffs: FDR < 0.05 | Log2FC > 0.5",
    pCutoff  = 0.05,
    FCcutoff = 0.5,
    pointSize = 2.5,
    labSize  = 4.5,
    col      = c("grey85", "grey75", "#fdae61", "#d73027"),
    colAlpha = 0.85,
    legendPosition = 'bottom',
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    max.overlaps = 25
  ) + 
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "grey95")
    )
  ggsave(file.path(config$results_dir, paste0("03_volcano_", tolower(gene), ".pdf")), p_volcano, width = 9, height = 8)
  
  # --- 2. Top 25 DEGs DotPlot (Faceted Significance Approach) ---
  print("Generating Top 25 DEGs Dot Plot (Ordered by p_val_adj)...")
  
  # 1. Select the top 25 genes strictly by adjusted p-value (FDR)
  top_degs_df <- dea_results %>%
    filter(p_val_adj < 0.05) %>%
    arrange(p_val_adj) %>%
    slice_head(n = 25)
  
  if(nrow(top_degs_df) > 0) {
    # 2. Add Directionality to sort them logically (Up vs Down)
    top_degs_df <- top_degs_df %>%
      mutate(Direction = ifelse(avg_log2FC > 0, "Up", "Down")) %>%
      arrange(desc(Direction), p_val_adj) # Upregulated first, then Downregulated, ordered by FDR
    
    top25_genes <- top_degs_df$gene
    
    # 3. Refine labels for X-axis logic
    macrophages$plot_status <- ifelse(macrophages$gene_status == "Positive", paste0(gene, "+"), paste0(gene, "-"))
    macrophages$plot_status <- factor(macrophages$plot_status, levels = c(paste0(gene, "+"), paste0(gene, "-")))
    
    # 4. Generate the Plot leveraging the binary nature of 2-group scaling
    p_dot <- DotPlot(macrophages, features = top25_genes, group.by = "plot_status", dot.scale = 9) +
      theme_classic(base_size = 14, base_family = "sans") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "italic", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 11, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey95", linetype = "solid")
      ) +
      scale_x_discrete(limits = top25_genes) + 
      # Classic divergent palette: Red for UP, Blue for DOWN.
      # The 2-group scaling will snap values naturally to these extremities.
      scale_color_gradient2(
        low = "royalblue", mid = "lightgrey", high = "firebrick",
        name = "Scaled Expression\n(Directionality)"
      ) +
      guides(size = guide_legend(title = "Percent\nExpressed")) +
      ggtitle(paste0("Top 25 DEGs Ranked by FDR: ", gene, " + vs ", gene, " -"))
    
    ggsave(file.path(config$results_dir, paste0("03_top25_dotplot_", tolower(gene), ".pdf")), p_dot, width = 12.5, height = 4)
  } else {
    print("Not enough significant DEGs to generate DotPlot.")
  }
  
  # --- 3. Dynamic Ridge Plots (Top 4 Validated Distributions) ---
  print("Generating Ridge Plots...")
  top4_ridge_genes <- dea_results %>% 
    filter(p_val_adj < 0.05) %>% 
    slice_max(abs(avg_log2FC), n = 4) %>% 
    pull(gene)
  
  if (length(top4_ridge_genes) > 0) {
    p_ridge <- RidgePlot(macrophages, features = top4_ridge_genes, group.by = "plot_status", 
                         cols = c("#F8766D", "#00BFC4"), ncol = 2) + 
      plot_annotation(title = paste0("Expression Density Validation (Top 4 DEGs)"),
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))
    
    ggsave(file.path(config$results_dir, paste0("03_ridge_plots_", tolower(gene), ".pdf")), p_ridge, width = 10, height = 6)
  }
  
  # ==============================================================================
  # GENE ONTOLOGY & METRICS EXPORT
  # ==============================================================================
  
  print(paste("Running GO (Biological Process) analysis for", gene, "..."))
  go_pool_genes <- dea_results %>% filter(p_val < 0.01 & avg_log2FC > 0.5) %>% pull("gene")
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
      write.csv(as.data.frame(ego),
                file.path(config$results_dir, paste0("03_", gene, "_GO_enrichment_results.csv")),
                row.names = FALSE)
      p_go <- dotplot(ego, showCategory = 15) + ggtitle(paste("GO Enrichment -", gene, "Exploratory Pool"))
      ggsave(file.path(config$results_dir, paste0("03_go_dotplot_", tolower(gene), ".pdf")),
             p_go, width = 8, height = 6)
      go_terms_count <- nrow(ego)
    }
  }
  
  # --- Export Summary Metrics ---
  summary_metrics <- list(
    pipeline_step        = "03_dea_and_go_MAST_Balanced",
    dataset_id           = config$dataset_id,
    gene_of_interest     = gene,
    method               = "MAST_cdr_adjusted",
    target_cluster_selected = target_cluster,
    pos_cells_analyzed   = length(cells_pos),
    neg_cells_downsampled = length(cells_neg_sampled),
    neg_cells_total_pool = length(cells_neg),
    strict_degs_fdr05    = sum(dea_results$p_val_adj < 0.05),
    exploratory_go_genes = length(go_pool_genes),
    go_terms_enriched    = go_terms_count
  )
  
  write_json(summary_metrics,
             file.path(config$results_dir, paste0("03_", gene, "_dea_summary.json")),
             pretty = TRUE, auto_unbox = TRUE)
  
  print(paste("DEA completed for", gene, "in dataset:", config$dataset_id))
}

print(paste("Pipeline Stage 03 completed for dataset:", config$dataset_id, ". Check", config$results_dir))