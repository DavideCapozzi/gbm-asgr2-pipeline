# scripts/03_dea_and_plotting.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(jsonlite)
library(EnhancedVolcano)
library(MAST)

# Default config for standalone execution (GSM3828672 behaviour preserved)
if (!exists("config")) {
  config <- list(
    dataset_id        = "GSM3828672",
    genes_of_interest = c("ASGR2", "CLEC10A"),
    violin_genes      = NULL,
    results_dir       = "results/GSM3828672",
    processed_dir     = "data/processed/GSM3828672"
  )
}

# 0. Ensure Output Directories Exist
dir.create(config$results_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load Clustered Data
gbm <- readRDS(file.path(config$processed_dir, "02_gbm_clustered.rds"))

# 2. Run DEA pipeline for each gene of interest
for (gene in config$genes_of_interest) {
  print(paste0("\n=== DEA for gene: ", gene, " ==="))

  # --- Identify Target Cluster (max gene+ proportion) ---
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
    warning(paste0("No cluster has >= 10 ", gene, "+ cells. Skipping DEA for this gene."))
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

  # --- Volcano Plot ---
  p_volcano <- EnhancedVolcano(
    dea_results,
    lab      = dea_results$gene,
    x        = 'avg_log2FC',
    y        = 'p_val_adj',
    title    = paste0(gene, "+ Signature (Balanced) - Cluster ", target_cluster),
    subtitle = 'MAST Model with CDR correction',
    pCutoff  = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize  = 4.0
  )
  ggsave(file.path(config$results_dir, paste0("03_volcano_", tolower(gene), ".pdf")),
         p_volcano, width = 10, height = 8)

  # --- Split Expression Plot: Frequency | Intensity ---
  # Reflects the two components of the MAST hurdle model:
  #   LEFT  panel: P(expression > 0) — the zero/non-zero component
  #   RIGHT panel: E[expression | expressed] — the continuous component
  # P-value annotation uses p_val_adj from dea_results directly (no recalculation).

  # Resolve which DEG(s) to plot
  violin_candidates <- if (!is.null(config$violin_genes) && length(config$violin_genes) > 0) {
    found <- intersect(config$violin_genes, dea_results$gene)
    if (length(found) == 0)
      warning(paste("None of config$violin_genes found in DEA results for", gene, "— using top DEG instead."))
    found
  } else {
    character(0)
  }

  if (length(violin_candidates) == 0 && nrow(dea_results) > 0) {
    # Base R ordering avoids S4Vectors::slice masking dplyr::slice after MAST loads
    ordered_idx <- order(dea_results$p_val_adj, -dea_results$avg_log2FC)
    violin_candidates <- dea_results$gene[ordered_idx[1]]
  }

  if (length(violin_candidates) > 0) {
    n_vln <- length(violin_candidates)

    # Extract normalised expression for all candidates at once
    expr_mat <- as.matrix(GetAssayData(macs_balanced, layer = "data")[violin_candidates, , drop = FALSE])
    plot_df <- do.call(rbind, lapply(violin_candidates, function(g) {
      data.frame(
        feature    = g,
        status     = factor(macs_balanced$gene_status, levels = c("Positive", "Negative")),
        expression = as.numeric(expr_mat[g, ]),
        stringsAsFactors = FALSE
      )
    }))

    # Map MAST p_val_adj to significance stars (reads dea_results, no recalculation)
    pval_to_stars <- function(p) {
      if (is.na(p) || p >= 0.05) "ns"
      else if (p < 0.001) "***"
      else if (p < 0.01)  "**"
      else                 "*"
    }

    fill_colors <- c("Positive" = "#F8766D", "Negative" = "#00BFC4")
    # ASCII-only title: avoids PDF font issues with Unicode dashes
    plot_title  <- paste0(gene, "+ vs ", gene, "- cells | Cluster ", target_cluster)

    base_theme <- theme_classic() +
      theme(
        legend.position = "none",
        axis.text.x     = element_text(size = 11),
        axis.title.y    = element_text(size = 9),
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 10)
      )

    # Build one frequency | intensity pair per candidate DEG
    plot_pairs <- lapply(violin_candidates, function(g) {

      gene_df <- plot_df[plot_df$feature == g, ]

      # --- Frequency data: % cells with expression > 0 per group ---
      pcts <- tapply(gene_df$expression > 0, gene_df$status,
                     function(x) round(mean(x) * 100))
      freq_df <- data.frame(
        status = factor(names(pcts), levels = c("Positive", "Negative")),
        pct    = as.numeric(pcts),
        stringsAsFactors = FALSE
      )

      # --- Intensity data: only cells with expression > 0 ---
      intens_df <- gene_df[gene_df$expression > 0, ]

      # --- P-value annotation: pulled from MAST output, never recalculated ---
      pval  <- dea_results$p_val_adj[match(g, dea_results$gene)]
      if (length(pval) == 0 || is.na(pval)) pval <- 1
      stars    <- pval_to_stars(pval)
      pval_lab <- paste0("p.adj = ", formatC(pval, format = "e", digits = 2))

      # Annotation y positions scaled to data range
      y_max   <- if (nrow(intens_df) > 0) max(intens_df$expression, na.rm = TRUE) else 1
      y_seg   <- y_max * 1.05   # bracket line just above violins
      y_pval  <- y_max * 1.12   # p-value text above bracket
      y_stars <- y_max * 1.20   # asterisks at the top

      # LEFT panel: frequency barplot (hurdle component 1)
      p_freq <- ggplot(freq_df, aes(x = status, y = pct, fill = status)) +
        geom_col(width = 0.55, color = "black", linewidth = 0.3) +
        geom_text(aes(label = paste0(pct, "%")),
                  vjust = -0.4, size = 3.5, fontface = "bold") +
        scale_fill_manual(values = fill_colors) +
        scale_y_continuous(
          limits = c(0, max(freq_df$pct, na.rm = TRUE) * 1.30),
          expand = expansion(mult = c(0, 0))
        ) +
        labs(x = NULL, y = "% Expressing cells", title = g) +
        base_theme

      # RIGHT panel base: violin + boxplot when enough data, strip chart otherwise
      if (nrow(intens_df) > 4) {
        p_intens <- ggplot(intens_df, aes(x = status, y = expression, fill = status)) +
          geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
          stat_summary(fun = median, geom = "crossbar",
                       width = 0.5, fatten = 0, color = "black", linewidth = 0.7)
      } else {
        # Sparse expressed cells: show individual points only
        p_intens <- ggplot(intens_df, aes(x = status, y = expression, fill = status)) +
          geom_jitter(shape = 21, width = 0.1, size = 2.5, alpha = 0.8, color = "grey30")
      }

      # Add MAST annotation, scales and theme to RIGHT panel
      p_intens <- p_intens +
        annotate("segment",
                 x = 1.05, xend = 1.95, y = y_seg, yend = y_seg,
                 color = "grey40", linewidth = 0.5) +
        annotate("text",
                 x = 1.5, y = y_pval, label = pval_lab,
                 size = 2.8, color = "grey35") +
        annotate("text",
                 x = 1.5, y = y_stars, label = stars,
                 size = 5.5, fontface = "bold") +
        scale_fill_manual(values = fill_colors) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.30))) +
        labs(x = NULL, y = "Expression (expressed cells only)",
             title = paste0(g, " (expressed cells)")) +
        base_theme

      # Combine horizontally: frequency 40%, intensity 60%
      p_freq + p_intens + plot_layout(widths = c(2, 3))
    })

    # Stack gene pairs vertically; Reduce handles n_vln == 1 correctly
    p_final <- Reduce(`/`, plot_pairs) +
      plot_annotation(
        title = plot_title,
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold",
                                    size = 12, margin = margin(b = 8))
        )
      )

    ggsave(
      file.path(config$results_dir, paste0("03_", gene, "_expression.pdf")),
      p_final, width = 7, height = 5 * n_vln
    )
    print(paste("Split expression plot saved for:", paste(violin_candidates, collapse = ", ")))
  } else {
    print(paste("No genes available for expression plot for", gene, "— skipping."))
  }

  # --- Gene Ontology Enrichment ---
  print(paste("Running GO (Biological Process) analysis for", gene, "..."))
  # Use string literal to avoid ambiguity with loop variable 'gene'
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
      ggsave(file.path(config$results_dir, paste0("03_", gene, "_go_dotplot.pdf")),
             p_go, width = 7, height = 6)
      go_terms_count <- nrow(ego)
    } else {
      print(paste("Warning: No significant GO terms found for", gene, "exploratory pool."))
    }
  } else {
    print(paste("Warning: Insufficient genes in exploratory pool for", gene, "GO enrichment."))
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
