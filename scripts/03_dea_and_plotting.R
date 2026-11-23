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

# Functional module table used to group/order the genes on the top-25 dot plot.
# Sourced unconditionally: it is a tiny file, and guarding on exists("GENE_MODULES")
# would silently ignore edits to it when re-sourcing this script in a live session.
source("scripts/gene_modules.R")

# --- 0. Config Safety & Environment Inheritance ---
if (!exists("config")) {
  config <- list(
    dataset_id        = "GSE162631",
    genes_of_interest = c("ASGR2"),
    violin_genes      = NULL,
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

# Fold-change mean function. Seurat's default for the `data` layer is
# log((rowSums(expm1(x)) + 1) / n, base = 2), i.e. it assumes NATURAL-log data.
# The log2tpm layer is log2(TPM + 1), so undo it with 2^x - 1 instead (written as
# expm1(x * log(2)) to keep the matrix sparse). Same formula otherwise, so the
# LogNormalize datasets keep Seurat's default (NULL).
fc_mean_fxn <- if (identical(config$norm_method, "log2tpm")) {
  function(x) log2((rowSums(expm1(x * log(2))) + 1) / NCOL(x))
} else NULL

# Multiple-testing label. FindMarkers' p_val_adj is a BONFERRONI correction over all
# genes in the object (Seurat default), not BH/FDR. Every plot caption and the JSON
# summary say so explicitly; the 0.05 threshold itself is unchanged.
P_ADJ_METHOD <- "bonferroni"
P_ADJ_LABEL  <- "adj. p (Bonferroni)"

# Module-faceted DotPlot shared by the top-25 and top-10 plots.
# `features` is the named list from build_module_feature_list(); genes in
# `sig_genes` get a "*" appended to their x-axis label.
module_dotplot <- function(obj, features, title, sig_genes = NULL) {
  # Passing `features` as a NAMED LIST makes Seurat facet the x-axis by module
  # (facet_grid(~feature.groups, scales = "free_x", space = "free_x")).
  p <- DotPlot(obj, features = features, group.by = "plot_status", dot.scale = 9) +
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
      panel.grid.major.x = element_line(color = "grey95", linetype = "solid"),
      # theme_classic() is applied after DotPlot, so it resets the blank strip
      # background DotPlot set - re-assert the facet styling here.
      strip.background = element_blank(),
      strip.text.x = element_text(size = 7.5, face = "bold", color = "black", lineheight = 0.9),
      # A module holding a single gene gets a one-gene-wide panel; let its label
      # overhang rather than be truncated mid-word.
      strip.clip = "off",
      panel.spacing = unit(0.7, "lines") # DotPlot's default of 1 line is too wide for ~10 panels
    ) +
    # NB: no scale_x_discrete(limits = ...) here - DotPlot already sets the factor
    # levels from `features`, and an explicit global `limits` would force every
    # gene into every facet panel under scales = "free_x".
    # Classic divergent palette: Red for UP, Blue for DOWN.
    # The 2-group scaling will snap values naturally to these extremities.
    scale_color_gradient2(low = "royalblue", mid = "lightgrey", high = "firebrick") +
    # DotPlot already added guides(color = guide_colorbar(title = "Average Expression")),
    # which overrides scale_color_gradient2(name = ...) - so set both titles here.
    guides(
      size  = guide_legend(title = "Percent\nExpressed"),
      color = guide_colorbar(title = "Scaled Expression\n(Directionality)")
    ) +
    ggtitle(title)

  if (length(sig_genes) > 0) {
    # Labels only (no limits), so each facet panel keeps just its own genes
    p <- p + scale_x_discrete(labels = function(x) ifelse(x %in% sig_genes, paste0(x, "*"), x))
  }
  p
}

# 2. Run DEA pipeline for each gene of interest
for (gene in config$genes_of_interest) {
  print(paste0("\n=== DEA for gene: ", gene, " ==="))

  # Per-gene output folder: results/<DATASET>/<gene>/
  gene_dir <- file.path(config$results_dir, tolower(gene))
  dir.create(gene_dir, showWarnings = FALSE, recursive = TRUE)

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

  # Display labels for the x/y axes of the dot, ridge and heatmap plots
  macrophages$plot_status <- ifelse(macrophages$gene_status == "Positive", paste0(gene, "+"), paste0(gene, "-"))
  macrophages$plot_status <- factor(macrophages$plot_status, levels = c(paste0(gene, "+"), paste0(gene, "-")))

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
    min.pct       = 0.25,
    mean.fxn      = fc_mean_fxn
  )
  
  # Structure results; exclude the gene itself to prevent double-dipping
  dea_results$gene_col <- rownames(dea_results)
  dea_results <- dea_results %>% filter(gene_col != gene)
  names(dea_results)[names(dea_results) == "gene_col"] <- "gene"
  
  write.csv(dea_results,
            file.path(gene_dir, paste0("03_", gene, "_DEA_results.csv")),
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
    caption  = paste0("Cutoffs: ", P_ADJ_LABEL, " < 0.05 | Log2FC > 0.5"),
    ylab     = bquote(~-Log[10] ~ "adj. p (Bonferroni)"),
    legendLabels = c("NS", "Log2 FC", "adj. p", "adj. p and Log2 FC"),
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
  ggsave(file.path(gene_dir, paste0("03_volcano_", tolower(gene), ".pdf")), p_volcano, width = 9, height = 8)

  # --- 2. Split Expression Plot: Frequency | Intensity ---
  # Reflects the two components of the MAST hurdle model:
  #   LEFT  panel: P(expression > 0) — the zero/non-zero component
  #   RIGHT panel: E[expression | expressed] — the continuous component
  # P-value annotation uses p_val_adj from dea_results directly (no recalculation).

  # Resolve which DEG(s) to plot: config$violin_genes overrides, else top DEG by adj. p (Bonferroni)
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
      pval_lab <- paste0("p.adj (Bonf.) = ", formatC(pval, format = "e", digits = 2))

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
      file.path(gene_dir, paste0("03_", gene, "_expression.pdf")),
      p_final, width = 7, height = 5 * n_vln
    )
    print(paste("Split expression plot saved for:", paste(violin_candidates, collapse = ", ")))
  } else {
    print(paste("No genes available for expression plot for", gene, "— skipping."))
  }

  # --- 3. Top 25 DEGs DotPlot (Faceted Significance Approach) ---
  print("Generating Top 25 DEGs Dot Plot (Ordered by p_val_adj)...")
  
  # 1. Select the top 25 genes strictly by adjusted p-value (Bonferroni)
  top_degs_df <- dea_results %>%
    filter(p_val_adj < 0.05) %>%
    arrange(p_val_adj) %>%
    slice_head(n = 25)
  
  if(nrow(top_degs_df) > 0) {
    # 2. Group genes by functional module (deterministic, see scripts/gene_modules.R).
    #    Direction is NOT used for sorting: with only two groups Seurat's per-gene
    #    scaling snaps every value to +/-0.707, so the red/blue colour already
    #    encodes direction completely. Grouping is the more informative axis order.
    top_degs_df <- top_degs_df %>%
      mutate(module = assign_gene_modules(gene)) %>%
      arrange(p_val_adj) # within each module panel: adj.-p ascending

    top25_features <- build_module_feature_list(top_degs_df$gene, top_degs_df$module)
    print("Functional module assignment for the top 25 DEGs:")
    print(table(top_degs_df$module))

    # 3. Generate the Plot leveraging the binary nature of 2-group scaling
    p_dot <- module_dotplot(macrophages, top25_features,
                            paste0("Top 25 DEGs by ", P_ADJ_LABEL, ", grouped by function: ", gene, " + vs ", gene, " -"))

    ggsave(file.path(gene_dir, paste0("03_top25_dotplot_", tolower(gene), ".pdf")), p_dot, width = 14, height = 5)
  } else {
    print("Not enough significant DEGs to generate DotPlot.")
  }

  # Nominal p-value ranking shared by the top-10 and top-100 plots. Unlike the
  # top-25 plot these do not filter on adj. p (so they also render for low-power
  # datasets); genes with Bonferroni adj. p < 0.05 are marked with "*" instead.
  # Base R ordering avoids S4Vectors::slice masking dplyr::slice after MAST loads
  ranked_df <- dea_results[order(dea_results$p_val, -abs(dea_results$avg_log2FC)), ]
  sig_genes <- ranked_df$gene[ranked_df$p_val_adj < 0.05]

  # --- 4. Top 10 Up / Top 10 Down DotPlot (module-faceted, one row per direction) ---
  print("Generating Top 10 Up / Down DEGs Dot Plot (Ordered by p_val)...")
  direction_sets <- list(
    Up   = head(ranked_df[ranked_df$avg_log2FC > 0, ], 10),
    Down = head(ranked_df[ranked_df$avg_log2FC < 0, ], 10)
  )

  updown_plots <- list()
  for (dir_name in names(direction_sets)) {
    dir_df <- direction_sets[[dir_name]]
    if (nrow(dir_df) == 0) next
    # Same deterministic module grouping as the top-25 plot; p-ascending within panels
    dir_features <- build_module_feature_list(dir_df$gene, assign_gene_modules(dir_df$gene))
    updown_plots[[dir_name]] <- module_dotplot(
      macrophages, dir_features,
      paste0(dir_name, " in ", gene, "+ (top ", nrow(dir_df), ")"),
      sig_genes = sig_genes
    )
  }

  if (length(updown_plots) > 0) {
    missing_dirs <- setdiff(names(direction_sets), names(updown_plots))
    updown_caption <- paste0("Top 10 per direction by nominal p-value; * = ", P_ADJ_LABEL, " < 0.05",
                             if (length(missing_dirs) > 0)
                               paste0(" | no genes with ", paste(missing_dirs, collapse = "/"), " fold change"))
    p_top10 <- Reduce(`/`, updown_plots) +
      plot_annotation(
        title   = paste0("Top 10 Up / Down DEGs, grouped by function: ", gene, " + vs ", gene, " -"),
        caption = updown_caption,
        theme   = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      )
    ggsave(file.path(gene_dir, paste0("03_top10_updown_dotplot_", tolower(gene), ".pdf")),
           p_top10, width = 12, height = 4.5 * length(updown_plots))
  } else {
    print("No DEGs available for Top 10 Up / Down DotPlot.")
  }

  # --- 5. Top 100 log2FC Heatmap (p-value order, no grouping) ---
  # Column 1: pooled MAST log2FC (balanced set, as in the DEA). Multi-sample
  # datasets add one column per patient: Positive vs Negative log2FC within that
  # sample, over all target-cluster cells (not the downsampled set) so each
  # patient keeps as many cells as possible. Too few cells -> NA (grey tile).
  print("Generating Top 100 log2FC Heatmap...")
  top100_df <- head(ranked_df, 100)

  if (nrow(top100_df) > 0) {
    fc_cols <- list()
    fc_cols[[paste0("Pooled (MAST)\n(+", length(cells_pos), " / -", length(cells_neg_sampled), ")")]] <-
      top100_df$avg_log2FC

    min_cells_patient <- 3
    if (length(config$samples) > 1 && "sample" %in% colnames(macrophages@meta.data)) {
      for (s in config$samples) {
        cells_s <- colnames(macrophages)[macrophages$sample == s]
        n_pos_s <- sum(macrophages$gene_status[cells_s] == "Positive")
        n_neg_s <- sum(macrophages$gene_status[cells_s] == "Negative")
        col_label <- paste0(s, "\n(+", n_pos_s, " / -", n_neg_s, ")")
        if (n_pos_s >= min_cells_patient && n_neg_s >= min_cells_patient) {
          fc_s <- FoldChange(subset(macrophages, cells = cells_s),
                             ident.1 = "Positive", ident.2 = "Negative",
                             features = top100_df$gene, mean.fxn = fc_mean_fxn)
          fc_cols[[col_label]] <- fc_s[top100_df$gene, "avg_log2FC"]
        } else {
          fc_cols[[col_label]] <- rep(NA_real_, nrow(top100_df))
        }
      }
    }

    # Long format for geom_tile; first gene (lowest p) at the top
    heat_df <- do.call(rbind, lapply(names(fc_cols), function(col) {
      data.frame(gene = top100_df$gene, column = col, log2FC = fc_cols[[col]],
                 stringsAsFactors = FALSE)
    }))
    heat_df$gene   <- factor(heat_df$gene, levels = rev(top100_df$gene))
    heat_df$column <- factor(heat_df$column, levels = names(fc_cols))
    # Symmetric colour limits capped at the 98th percentile of |log2FC|: a few
    # extreme per-patient values (small groups) would otherwise wash out every
    # other tile. Values beyond the cap are squished to the end colour; the
    # matrix CSV keeps the raw values.
    fc_abs <- abs(heat_df$log2FC)
    fc_lim <- as.numeric(quantile(fc_abs, 0.98, na.rm = TRUE))
    fc_capped <- any(fc_abs > fc_lim, na.rm = TRUE)

    p_heat <- ggplot(heat_df, aes(x = column, y = gene, fill = log2FC)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "royalblue", mid = "white", high = "firebrick", midpoint = 0,
                           limits = c(-fc_lim, fc_lim), oob = scales::squish, na.value = "grey85",
                           name = "log2FC") +
      scale_y_discrete(labels = function(x) ifelse(x %in% sig_genes, paste0(x, "*"), x)) +
      scale_x_discrete(position = "top", expand = c(0, 0)) +
      labs(x = NULL, y = NULL,
           title = paste0("Top ", nrow(top100_df), " DEGs by p-value\n", gene, "+ vs ", gene, "-"),
           caption = paste0("Rows ordered by nominal p-value\n* = ", P_ADJ_LABEL, " < 0.05",
                            if (fc_capped)
                              paste0("\nColour capped at |log2FC| = ", round(fc_lim, 2), " (98th pct)"),
                            if (length(fc_cols) > 1)
                              paste0("\nPatient columns: all cluster ", target_cluster, " cells\n",
                                     "grey = < ", min_cells_patient, " cells in a group"))) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 6, face = "italic", color = "black"),
        axis.text.x = element_text(size = 8, face = "bold", color = "black"),
        panel.grid = element_blank(),
        # Centre title/caption on the whole figure, not the (narrow) tile panel
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.caption = element_text(hjust = 0, size = 7, color = "grey30")
      )
    if (length(fc_cols) > 1) {
      # Separate the pooled column from the per-patient ones
      p_heat <- p_heat + geom_vline(xintercept = 1.5, linewidth = 0.6, color = "black")
    }

    ggsave(file.path(gene_dir, paste0("03_top100_fc_heatmap_", tolower(gene), ".pdf")),
           p_heat, width = max(4.5, 2.5 + 1 * length(fc_cols)), height = 14)

    # Data behind the heatmap: one log2FC column per heatmap column
    fc_matrix <- data.frame(gene = top100_df$gene, p_val = top100_df$p_val,
                            p_val_adj = top100_df$p_val_adj, stringsAsFactors = FALSE)
    for (col in names(fc_cols)) {
      fc_matrix[[paste0("log2FC_", sub("\n.*", "", sub(" \\(MAST\\)", "", col)))]] <- fc_cols[[col]]
    }
    write.csv(fc_matrix,
              file.path(gene_dir, paste0("03_", gene, "_top100_fc_matrix.csv")),
              row.names = FALSE)
  } else {
    print("No DEGs available for Top 100 heatmap.")
  }

  # --- 6. Dynamic Ridge Plots (Top 4 Validated Distributions) ---
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
    
    ggsave(file.path(gene_dir, paste0("03_ridge_plots_", tolower(gene), ".pdf")), p_ridge, width = 10, height = 6)
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
                file.path(gene_dir, paste0("03_", gene, "_GO_enrichment_results.csv")),
                row.names = FALSE)
      p_go <- dotplot(ego, showCategory = 15, label_format = 40) + ggtitle(paste("GO Enrichment -", gene, "Exploratory Pool"))
      ggsave(file.path(gene_dir, paste0("03_go_dotplot_", tolower(gene), ".pdf")),
             p_go, width = 8, height = 10)
      go_terms_count <- nrow(ego)
    }
  }
  
  # --- Export Summary Metrics ---
  summary_metrics <- list(
    pipeline_step        = "03_dea_and_go_MAST_Balanced",
    dataset_id           = config$dataset_id,
    gene_of_interest     = gene,
    method               = "MAST_cdr_adjusted",
    fc_mean_fxn          = if (is.null(fc_mean_fxn)) "seurat_default" else "log2_data_corrected",
    target_cluster_selected = target_cluster,
    pos_cells_analyzed   = length(cells_pos),
    neg_cells_downsampled = length(cells_neg_sampled),
    neg_cells_total_pool = length(cells_neg),
    p_adjust_method      = P_ADJ_METHOD,
    strict_degs_padj05   = sum(dea_results$p_val_adj < 0.05),
    exploratory_go_genes = length(go_pool_genes),
    go_terms_enriched    = go_terms_count
  )
  
  write_json(summary_metrics,
             file.path(gene_dir, paste0("03_", gene, "_dea_summary.json")),
             pretty = TRUE, auto_unbox = TRUE)
  
  print(paste("DEA completed for", gene, "in dataset:", config$dataset_id))
}

print(paste("Pipeline Stage 03 completed for dataset:", config$dataset_id, ". Check", config$results_dir))