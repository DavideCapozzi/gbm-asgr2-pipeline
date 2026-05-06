# scripts/main.R
# Usage:
#   Rscript scripts/main.R GSM3828672
#   Rscript scripts/main.R GSE162631

DATASET <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(DATASET)) DATASET <- "GSM3828672"

configs <- list(
  GSM3828672 = list(
    dataset_id        = "GSM3828672",
    format            = "tpm_tsv",
    data_path         = "data/raw/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz",
    project_name      = "GBM_SmartSeq2",
    norm_method       = "log2tpm",
    feature_method    = "dispersion",
    nfeatures         = 2000,
    qc_min_feat       = 1000,
    qc_max_feat       = 8000,
    qc_max_mt         = 15,
    scale_all         = TRUE,
    use_harmony       = FALSE,
    harmony_var       = NULL,
    samples           = NULL,
    genes_of_interest = c("ASGR2", "CLEC10A"),
    violin_genes      = NULL,
    results_dir       = "results/GSM3828672",
    processed_dir     = "data/processed/GSM3828672"
  ),
  GSE162631 = list(
    dataset_id        = "GSE162631",
    format            = "10x_nested_zip",
    data_path         = "data/raw/GSE162631_raw_counts_matrix.zip.gz",
    project_name      = "GBM_10X_GSE162631",
    norm_method       = "lognorm",
    feature_method    = "vst",
    nfeatures         = 2000,
    qc_min_feat       = 200,
    qc_max_feat       = 6000,
    qc_max_mt         = 20,
    scale_all         = FALSE,
    use_harmony       = TRUE,
    harmony_var       = "sample",
    samples           = c("R1_T", "R2_T", "R3_T", "R4_T"),
    genes_of_interest = c("ASGR2", "CLEC10A"),
    violin_genes      = NULL,
    results_dir       = "results/GSE162631",
    processed_dir     = "data/processed/GSE162631"
  )
)

if (!DATASET %in% names(configs)) {
  stop(paste0("Unknown dataset: '", DATASET, "'. Valid options: ", paste(names(configs), collapse = ", ")))
}

config <- configs[[DATASET]]
assign("config", config, envir = .GlobalEnv)

# Pre-flight checks — fail early with a clear message rather than mid-pipeline
if (isTRUE(config$use_harmony) && !requireNamespace("harmony", quietly = TRUE)) {
  stop("Package 'harmony' is required for dataset '", DATASET, "' but is not installed.\n",
       "Install it with: conda install -n gbm-asgr2-pipeline r-harmony\n",
       "Or from within R: install.packages('harmony')")
}

message(paste0("\n=== Running pipeline for dataset: ", DATASET, " ===\n"))

source("scripts/01_preprocessing.R")
source("scripts/02_clustering_annotation.R")
source("scripts/03_dea_and_plotting.R")

message(paste0("\n=== Pipeline completed for dataset: ", DATASET, " ==="))
message(paste0("Results: ", config$results_dir))
message(paste0("Processed data: ", config$processed_dir))
