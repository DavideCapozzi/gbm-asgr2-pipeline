# Candidate GBM single-cell datasets for the CLEC10A / lectin extension

Screened 2026-11-23. Accessions, file names and sizes were read from GEO
(`acc.cgi` SOFT text and `suppl/filelist.txt`) unless tagged `[UNVERIFIED]`.

## Selection criteria

1. **GBM** (IDH-wildtype adult, or with an IDH-wt subset that can be selected).
2. **Comparable** with the existing analysis: human, droplet or plate scRNA-seq with enough
   myeloid cells to form CLEC10A+ / CLEC10A− groups (≥ 20 positive cells in one compartment).
3. **Pipeline-friendly**: a count (or TPM) matrix that `01_preprocessing.R` can load with at
   most a small new `format` branch. Controlled-access raw data (EGA) only is a no.

## Already used

| Accession | Study | Tech | Notes |
|---|---|---|---|
| GSM3828672 (in GSE131928) | Neftel 2019, Cell | Smart-seq2, TPM | 7,427 cells; includes CD45− sorted samples (MGH128–152); CLEC10A target cluster 14 was astrocytic → invalid, see `PLAN_clec10a_consensus.md` |
| GSE162631 | Xie 2021 | 10X, raw counts | R1_T–R4_T tumour, 59,717 cells; CLEC10A target cluster 13 = cDC2 |

## Ranking

| # | Accession | Study | Cells / patients (GBM) | Sorting | Files on GEO | Pipeline effort | Priority |
|---|---|---|---|---|---|---|---|
| 1 | **GSE163120** | Pombo Antunes 2021, Nat Neurosci | ND: ~22k cells, 7 pts; R: ~43k cells, 5 samples; CITE-seq: ~17k cells, 3 pts (R2, R5, ND8) | **CD45+** | dense CSV gene×cell + author annotations (41 MB, 83 MB, 53 MB) | new `csv_dense` loader | **High — do first** |
| 2 | **GSE182109** | Abdelfattah 2022, Nat Commun | 201,986 cells, 44 samples, 18 pts (11 ndGBM, 5 rGBM, 2 LGG; multi-region) | none (whole TME) | per-sample 10X mtx triplets, RAW.tar 2.46 GB | new `10x_flat_mtx` loader (Read10X after regrouping files) | **High** |
| 3 | GSM3828673 (in GSE131928) | Neftel 2019 | 16,201 cells, 9 pts; 268 CLEC10A+ cells, 2,818 PTPRC+ | none | TPM TSV, **already in `data/raw/GSE131928_RAW.tar`** | none (`tpm_tsv`) | Medium — not independent |
| 4 | GSE235676 | "Siglec-9 acts as an immune checkpoint molecule on macrophages in glioblastoma" (2023) `[authors/journal UNVERIFIED]` | 24 pts (Pt1–18 10X, Pt19–24 Singleron) | `[UNVERIFIED]` | per-pt mtx, RAW.tar 1.15 GB; integrated RDS 6.2 GB | `10x_flat_mtx`; IDH status to check | Medium — lectin-network angle |
| 5 | GSE154795 | Lee 2021, Nat Commun | 40 transcriptomes, ndGBM + rGBM ± neoadjuvant anti-PD-1 | CD45+ | processed TXT 1.26 GB; RDS 5.95 GB | format of TXT to check | Medium — immunotherapy angle |
| 6 | GSE222520 | MD Anderson, IDH-stratified gliomas | IDH-wt primary 4, IDH-wt recurrent 4 (+ IDH-mut 10, non-glioma brain NGB 3) | CD45+ | per-sample tar.gz, 1.1 GB | `[UNVERIFIED]` content of tarballs | Medium — GBM vs IDH-mut vs non-tumour |
| 7 | GSE84465 | Darmanis 2017, Cell Rep | 3,589 cells, 4 pts, core + periphery | none | single counts CSV, 20.5 MB | new `counts_csv` loader | Low — same platform as GSM3828672, few myeloid cells |

### Not recommended now

| Accession / resource | Reason |
|---|---|
| GBmap (Ruiz-Moreno 2025, Neuro-Oncology; own data GSE211376) — 1.1 M cells, 240 pts, 26 datasets | Best for a final atlas-level check, but h5ad + harmonised labels, outside the pipeline. Overlaps with most datasets above. CELLxGENE availability `[UNVERIFIED]` |
| Ravi 2022, Cancer Cell — Visium, Dryad doi:10.5061/dryad.h70rxwdmj (~10 GB) | Spatial modality; good later for CLEC10A+ myeloid niches vs hypoxia/Tn-synthesis, not for DEA replication |
| Klemm 2020, Cell — sorted immune populations, bulk RNA-seq, 48 pts | Bulk, not sc; quick sanity check of CLEC10A by sorted population if data are accessible `[UNVERIFIED access]` |
| GSE138794 (Wang 2019) | Mixed gliomas, sc + sn, raw data restricted |
| GSE117891 (Yu 2020) | STRT-seq, mixed gliomas incl. metastasis, ~6k cells |
| GSE135437 (Sankowski 2019) | mCEL-Seq2, odd per-plate format, mostly microglia |
| GSE173278 | Mostly tumour cells, organoids, cell lines |
| Mouse datasets (e.g. GL261 in GSE163120) | Human MGL ≠ mouse Mgl1/Mgl2; not comparable |

## Dataset notes

### GSE163120 — first choice

- Files (human): `GSM4972211_Human.GBM.ND1_2_3_4_5_6_7.filtered.gene.bc.matrix.csv.gz`,
  `GSM4972210_Human.GBM.R1_2_3_4_4nc...`, `GSM4972212_Citeseq_Human.GBM.R2_5_ND8.filtered.{RNA,ADT}...`
  plus `annot.*.csv.gz` (UMAP, cluster label, sample per barcode). Human raw reads are on EGA
  (EGAS00001004871) — not needed.
- Author labels `[VERIFIED from annot files]`: ND — TAM 1, TAM 2, prol. TAM, Monocytes, DC,
  T, NK, B. R — adds DC 1–4, Treg, Mast, Plasma B. CITE-seq — **cDC2 (761 cells)**, cDC1, pDC,
  TAM 1–5, Mono 1–2, …
  → independent check of our marker-based annotation (02b).
- **CITE-seq panel includes CD301 (MGL/CLEC10A protein)** `[VERIFIED from ADT matrix rownames]`,
  plus **CD371 (CLEC12A)**, CD163, CD206, CD1c, CD49d, TMEM119, LOX-1, HLA-DR, CD273 (PD-L2),
  CD370 (CLEC9A), CD209 (DC-SIGN), FcεRIα. This lets us
  1. check CLEC10A mRNA against CD301 protein in the same cell;
  2. reproduce the Pace 2026 flow gating (macrophage CD49d+ vs microglia TMEM119+) on
     sequenced cells.
  PD-L1 (CD274) is **not** in the panel.
- To check before use: whether the matrices are raw UMI counts (expected: integers).

### GSE182109

- 10X v3 per sample: `GSM55185xx_<sample>_{barcodes,features,matrix}` (e.g. `ndGBM-01-A`,
  `rGBM-01-B`); one ndGBM sample (`ndGBM-01-A`) has a different features file size — check gene
  annotation consistency.
- Whole TME (tumour + immune) → the only high-priority dataset where **tumour-cell
  Tn-synthesis genes** and CLEC10A+ myeloid cells can be related **per patient**.
- Multi-region: integrate by patient (Harmony on `patient`, not on sample), exclude LGG.
- ~200k cells: within 62 GB RAM if `scale_all = FALSE` as for GSE162631.

### GSM3828673

- Patients 102, 105, 114, 115, 118, 124, 125, 126, 143. **Six of nine are also in the
  Smart-seq2 set** (102, 105, 115, 124, 125, 143) → technical replicate, not an independent
  cohort. Cheapest possible run (format exists, file on disk), useful as a sanity check.
