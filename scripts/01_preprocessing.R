# scripts/01_preprocessing.R
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(Matrix)
library(jsonlite)

# Default config for standalone execution (GSM3828672 behaviour preserved)
if (!exists("config")) {
  config <- list(
    dataset_id     = "GSM3828672",
    format         = "tpm_tsv",
    data_path      = "data/raw/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz",
    project_name   = "GBM_SmartSeq2",
    norm_method    = "log2tpm",
    feature_method = "dispersion",
    nfeatures      = 2000,
    qc_min_feat    = 1000,
    qc_max_feat    = 8000,
    qc_max_mt      = 15,
    scale_all      = TRUE,
    use_harmony    = FALSE,
    harmony_var    = NULL,
    samples        = NULL,
    results_dir    = "results/GSM3828672",
    processed_dir  = "data/processed/GSM3828672"
  )
}

# 0. Ensure Output Directories Exist
dir.create(config$results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$processed_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: find the directory containing the 10X files (matrix.mtx.gz), excluding macOS
# artifacts (__MACOSX) which unzip leaves behind and whose filenames also match the pattern.
find_10x_dir <- function(base_dir) {
  mtx_files <- list.files(base_dir, pattern = "matrix\\.mtx(\\.gz)?$",
                          recursive = TRUE, full.names = TRUE)
  mtx_files <- mtx_files[!grepl("__MACOSX", mtx_files)]
  if (length(mtx_files) == 0) stop("No matrix.mtx(.gz) found under: ", base_dir)
  dirname(mtx_files[1])
}

# 1. Load Data
if (config$format == "tpm_tsv") {
  if (!file.exists(config$data_path)) stop("Data file not found. Check path.")

  # Read table and convert to sparse matrix immediately to suppress Seurat coercion warnings
  tpm_df <- read.table(gzfile(config$data_path), header = TRUE, row.names = 1, sep = "\t")
  tpm_matrix <- as(as.matrix(tpm_df), "dgCMatrix")
  rm(tpm_df)

  gbm <- CreateSeuratObject(counts = tpm_matrix, project = config$project_name)
  rm(tpm_matrix)

} else if (config$format == "10x_nested_zip") {
  if (!file.exists(config$data_path)) stop("Data file not found: ", config$data_path)

  extract_dir <- file.path(dirname(config$data_path), "GSE162631")

  # Extract archive on first run only (idempotent check on first sample directory)
  if (!dir.exists(file.path(extract_dir, config$samples[1]))) {
    message("Extracting GSE162631 archive (one-time operation, may take several minutes)...")
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

    # Decompress outer gz → zip (stream via system gunzip, avoids full in-memory load of 800MB)
    tmp_zip <- tempfile(fileext = ".zip")
    ret <- system2("gunzip", args = c("-c", shQuote(config$data_path)), stdout = tmp_zip)
    if (ret != 0) { unlink(tmp_zip); stop("gunzip failed for: ", config$data_path) }

    # Discover the prefix directory inside the outer zip (e.g. "raw_counts_matrix")
    # sub("/.*", "") strips everything from the first "/" onward → extracts the top-level dir name
    avail <- unzip(tmp_zip, list = TRUE)$Name
    zip_dirs <- unique(sub("/.*", "", avail[grepl("/", avail)]))
    zip_prefix <- if (length(zip_dirs) == 1) zip_dirs else ""

    for (s in config$samples) {
      sample_zip_entry <- if (zip_prefix == "") paste0(s, ".zip") else paste0(zip_prefix, "/", s, ".zip")
      if (!sample_zip_entry %in% avail) {
        unlink(tmp_zip)
        stop("Expected '", sample_zip_entry, "' not found inside archive. Available: ",
             paste(avail, collapse = ", "))
      }
      unzip(tmp_zip, files = sample_zip_entry, exdir = extract_dir)

      inner_zip <- if (zip_prefix == "")
        file.path(extract_dir, paste0(s, ".zip"))
      else
        file.path(extract_dir, zip_prefix, paste0(s, ".zip"))

      sample_dir <- file.path(extract_dir, s)
      dir.create(sample_dir, showWarnings = FALSE)
      unzip(inner_zip, exdir = sample_dir)
      unlink(inner_zip)

      # Remove macOS artifacts injected by zip inside each sample directory
      mac_in_sample <- file.path(sample_dir, "__MACOSX")
      if (dir.exists(mac_in_sample)) unlink(mac_in_sample, recursive = TRUE)
    }

    unlink(tmp_zip)
    if (zip_prefix != "") {
      prefix_dir <- file.path(extract_dir, zip_prefix)
      if (dir.exists(prefix_dir)) unlink(prefix_dir, recursive = TRUE)
    }
    macosx_dir <- file.path(extract_dir, "__MACOSX")
    if (dir.exists(macosx_dir)) unlink(macosx_dir, recursive = TRUE)
    message("Extraction complete.")
  }

  # Load each sample into a Seurat object
  seurat_list <- lapply(config$samples, function(s) {
    sample_dir <- file.path(extract_dir, s)

    # Clean up any macOS artifacts that survived a previous extraction (idempotent)
    mac_in_sample <- file.path(sample_dir, "__MACOSX")
    if (dir.exists(mac_in_sample)) unlink(mac_in_sample, recursive = TRUE)

    tenx_dir <- find_10x_dir(sample_dir)

    # Cell Ranger v3 names the gene annotation file genes.tsv.gz; Seurat 5 requires features.tsv.gz.
    # Rename in-place so subsequent runs are also correct (file.rename is a no-op if already done).
    genes_f    <- file.path(tenx_dir, "genes.tsv.gz")
    features_f <- file.path(tenx_dir, "features.tsv.gz")
    if (file.exists(genes_f) && !file.exists(features_f)) {
      file.rename(genes_f, features_f)
    }

    counts <- Read10X(data.dir = tenx_dir)
    obj <- CreateSeuratObject(counts = counts, project = config$project_name,
                              min.cells = 3, min.features = 200)
    obj$sample <- s
    obj
  })

  # Merge all samples; add.cell.ids prevents barcode collisions across samples
  gbm <- merge(seurat_list[[1]], y = seurat_list[-1],
               add.cell.ids = config$samples, project = config$project_name)
  rm(seurat_list)

  # Seurat 5: join per-sample layers into a single counts layer before downstream operations
  gbm <- JoinLayers(gbm)

} else {
  stop("Unknown config$format: '", config$format, "'. Expected 'tpm_tsv' or '10x_nested_zip'.")
}

# Track pre-filter cell count for QC metrics
cells_pre_filter <- ncol(gbm)

# 3. Calculate QC Metrics Dynamically
mt_genes <- grep(pattern = "^MT-", x = rownames(gbm), value = TRUE)
vln_features <- c("nFeature_RNA", "nCount_RNA")

if (length(mt_genes) > 0) {
  gbm[["percent.mt"]] <- PercentageFeatureSet(gbm, pattern = "^MT-")
  vln_features <- c(vln_features, "percent.mt")
} else {
  warning("No mitochondrial genes found matching '^MT-'. Setting percent.mt to 0.")
  gbm[["percent.mt"]] <- 0
}

# 4. Save pre-filtering QC plots
p_qc_pre <- VlnPlot(
  gbm,
  features = vln_features,
  ncol = length(vln_features),
  pt.size = 0.1,
  layer = "counts"
)
ggsave(file.path(config$results_dir, "01_qc_pre_filter_vln.pdf"), p_qc_pre, width = 12, height = 6)

# 5. Apply QC Filters
gbm <- subset(gbm,
              subset = nFeature_RNA > config$qc_min_feat &
                       nFeature_RNA < config$qc_max_feat &
                       percent.mt < config$qc_max_mt)

# 6. Normalization
if (config$norm_method == "log2tpm") {
  # TPM is already a relative measure; apply log2(x+1) directly on the counts layer
  counts_mat <- GetAssayData(gbm, layer = "counts")
  counts_mat@x <- log2(counts_mat@x + 1)
  gbm <- SetAssayData(gbm, layer = "data", new.data = counts_mat)
  rm(counts_mat)
} else {
  # Standard log-normalization for raw UMI counts (10X)
  gbm <- NormalizeData(gbm, normalization.method = "LogNormalize", scale.factor = 10000)
}

# 7. Feature Selection and Scaling
# dispersion is appropriate for TPM; vst is appropriate for raw UMI counts
gbm <- FindVariableFeatures(gbm, selection.method = config$feature_method, nfeatures = config$nfeatures)

# Scale all genes for small datasets (SmartSeq2); scale only HVGs for large 10X to avoid OOM
features_to_scale <- if (isTRUE(config$scale_all)) rownames(gbm) else VariableFeatures(gbm)
gbm <- ScaleData(gbm, features = features_to_scale)

# 8. Generate Machine-Friendly QC Metrics (JSON)
cells_post_filter <- ncol(gbm)
qc_metrics <- list(
  pipeline_step        = "01_preprocessing",
  dataset_id           = config$dataset_id,
  cells_total_raw      = cells_pre_filter,
  cells_total_filtered = cells_post_filter,
  survival_rate_pct    = round((cells_post_filter / cells_pre_filter) * 100, 2),
  median_nFeature_RNA  = median(gbm$nFeature_RNA),
  median_nCount_RNA    = median(gbm$nCount_RNA),
  median_percent_mt    = median(gbm$percent.mt)
)

write_json(qc_metrics, file.path(config$results_dir, "01_qc_metrics.json"), pretty = TRUE, auto_unbox = TRUE)

# Save checkpoint
saveRDS(gbm, file.path(config$processed_dir, "01_gbm_preprocessed.rds"))
print(paste("Script 01 completed successfully for dataset:", config$dataset_id))
