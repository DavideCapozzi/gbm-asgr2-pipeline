# Pace, Alfano et al. 2026 — ASGR2 and CLEC12A as C-type lectin hubs in GBM

`[VERIFIED 2026-11-23 against the PMC full text; supplementary tables NOT read (MDPI returned 403)]`

| | |
|---|---|
| Citation | Pace A, Alfano C, D'Angelo L, Napoletano C, Zizzari IG, Santoro A, Nuti M, Farina L, Petti M, Rughetti A. *ASGR2 and CLEC12A as Prognostically Relevant C-Type Lectin Hubs in Glioblastoma.* Int J Mol Sci. 2026;27(6):2626 |
| IDs | DOI 10.3390/ijms27062626 · PMID 41898489 · PMCID PMC13026556 |
| Type | Original research: bulk transcriptomics network analysis + flow cytometry |

**Important:** this paper is about **ASGR2 and CLEC12A**, not MGL. CLEC10A is not among its
hubs, and it is mentioned only in the discussion (citing Dusoswa 2020). Whether CLEC10A was in
the 39-lectin network is unknown: the lectin list is in Supplementary Table S2, not read.

## Data and methods

| Item | Detail |
|---|---|
| Tumour cohort | TCGA-GBM, 152 primary → **145 IDH-wt** (RNA-seq + clinical) |
| Validation cohort | CGGA mRNAseq_693, **109** non-IDH-mutant primary GBM (abstract says 133) |
| Control | GTEx v8 brain cortex, **255** samples |
| Lectin panel | 67 genes (14 Siglecs, 15 galectins, 38 C-type lectins); 57 in TCGA → **39** after filter (sum < 500 counts removed); 36 in CGGA (SIGLEC5, CLEC5A, CD207 dropped) |
| Normalisation | DESeq2, log2 |
| DE | t-test, BH, \|log2FC\| ≥ 1.5, adj-p < 0.05 |
| Differential co-expression | Spearman per pair, positive links only; Fisher Z; Z = (zT − zC)/√(1/(nT−3) + 1/(nC−3)); threshold Z > 3; hubs = top 20% by degree |
| Survival | Fuzzy C-means (ppclust) on ASGR2 + CLEC12A; Kaplan–Meier (survival, survminer) |
| Immune deconvolution | TIMER2.0, purity-adjusted Spearman |
| Wet lab | Flow cytometry on **20** IDH-wt GBM patients, tumour + matched PBMC |
| Flow panel | CD45, CD15, CD14, CD66b, P2RY12, CD49d, CLEC12A, LOX1, CD11c, ASGR2 (macrophage vs microglia separated by CD49d / P2RY12) |

## Results (as reported)

- TCGA-vs-GTEx hubs: **ASGR2, LGALS8, CLEC4A, SIGLEC7, CLEC7A, SIGLEC10, SIGLEC11, CLEC12A**.
  ASGR2, LGALS8, CLEC4A are hubs **without** being differentially expressed.
- Hub degree TCGA vs CGGA: Spearman 0.6. ASGR2 and CLEC12A are hubs in both networks.
- TCGA survival: neither lectin alone is prognostic; **low ASGR2 + low CLEC12A** → better OS
  (76 low vs 66 high patients, log-rank p = 0.039).
- TIMER2.0: both correlate positively with macrophages, neutrophils, myeloid DCs; CLEC12A also
  with CD8+ T cells.
- Flow, tumour: ASGR2 higher on infiltrating **macrophages than microglia**
  (0.34 ± 0.06 vs 0.06 ± 0.01, p = 0.0094); CLEC12A on all myeloid subsets.
  ASGR2^low CLEC12A^low on CD45+ → better OS (HR 6.678, 95% CI 1.092–40.84, p = 0.0249).
- Flow, blood: ASGR2 **undetectable** on circulating cells (tumour-restricted); CLEC12A higher
  on circulating CD45+, monocytes, M-MDSC than in tumour. High CLEC12A+ CD45+ (> 37%) and
  CLEC12A+ M-MDSC (> 7.4%) → worse OS (p = 0.0302, p = 0.0169).
- Stated working hypothesis: **lectin synergy** in myeloid cells sustains immunosuppression.

## Limitations (ours to address)

- Bulk RNA: co-expression across patients cannot say whether two lectins are in the **same
  cell**; TIMER is a deconvolution estimate.
- Small flow cohort (n = 20; wide HR CIs).
- No single-cell transcriptomics.

## What our single-cell data can test

| Claim / hypothesis in the paper | sc-RNA test | Best dataset |
|---|---|---|
| ASGR2 higher in macrophages than microglia | ASGR2 by annotated macrophage (CD49d/ITGA4+, P2RY12−) vs microglia (P2RY12+, TMEM119+) | GSE163120 (CITE-seq has CD49d and TMEM119 protein), GSE182109 |
| Lectin "synergy" in myeloid cells | single-cell co-expression of hub lectins within myeloid compartments, per patient | all |
| CLEC12A on monocytes / M-MDSC | CLEC12A by monocyte vs TAM cluster; CD371 protein | GSE163120 CITE-seq (CD371 = CLEC12A) |
| Hubs drive immunosuppression | hub-lectin+ vs − cells in the same compartment against a pre-defined suppression signature | all |

The sc data can also place **CLEC10A** relative to these hubs (same cells as ASGR2? as CLEC12A?),
which is the link between this paper and the MGL story.
