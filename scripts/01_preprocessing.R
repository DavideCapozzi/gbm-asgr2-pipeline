# scripts/01_preprocessing.R
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(Matrix)
library(jsonlite) 

# 0. Ensure Output Directories Exist
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# 1. Load Data
data_path <- "data/raw/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz"
if (!file.exists(data_path)) stop("Data file not found. Check path.")

# Read table and convert to sparse matrix immediately to suppress Seurat coercion warnings
tpm_df <- read.table(gzfile(data_path), header = TRUE, row.names = 1, sep = "\t")
tpm_matrix <- as(as.matrix(tpm_df), "dgCMatrix")
rm(tpm_df) # Free memory

# 2. Initialize Seurat Object
gbm <- CreateSeuratObject(counts = tpm_matrix, project = "GBM_SmartSeq2")

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
ggsave("results/01_qc_pre_filter_vln.pdf", p_qc_pre, width = 12, height = 6)

# 5. Apply QC Filters
gbm <- subset(gbm, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 15)

# 6. Normalization (Log2 TPM + 1)
counts_mat <- GetAssayData(gbm, layer = "counts")
counts_mat@x <- log2(counts_mat@x + 1)
gbm <- SetAssayData(gbm, layer = "data", new.data = counts_mat)
rm(counts_mat) # Free memory

# 7. Feature Selection and Scaling
# Use dispersion instead of vst as vst expects raw counts, while we have TPM
gbm <- FindVariableFeatures(gbm, selection.method = "dispersion", nfeatures = 2000)
gbm <- ScaleData(gbm, features = rownames(gbm))

# 8. Generate Machine-Friendly QC Metrics (JSON)
cells_post_filter <- ncol(gbm)
qc_metrics <- list(
  pipeline_step = "01_preprocessing",
  cells_total_raw = cells_pre_filter,
  cells_total_filtered = cells_post_filter,
  survival_rate_pct = round((cells_post_filter / cells_pre_filter) * 100, 2),
  median_nFeature_RNA = median(gbm$nFeature_RNA),
  median_nCount_RNA = median(gbm$nCount_RNA),
  median_percent_mt = median(gbm$percent.mt)
)

write_json(qc_metrics, "results/01_qc_metrics.json", pretty = TRUE, auto_unbox = TRUE)

# Save checkpoint
saveRDS(gbm, "data/processed/01_gbm_preprocessed.rds")
print("Script 01 completed successfully. QC metrics exported to JSON.")