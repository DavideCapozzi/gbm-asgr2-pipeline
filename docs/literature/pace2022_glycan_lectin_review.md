# Pace et al. 2022 — Glycan–lectin interactions as immunosuppression drivers in GBM

`[VERIFIED 2026-11-23 against the PMC full text]`

| | |
|---|---|
| Citation | Pace A, Scirocchi F, Napoletano C, Zizzari IG, D'Angelo L, Santoro A, Nuti M, Rahimi H, Rughetti A. *Glycan-Lectin Interactions as Novel Immunosuppression Drivers in Glioblastoma.* Int J Mol Sci. 2022;23(11):6312 |
| IDs | DOI 10.3390/ijms23116312 · PMID 35682991 · PMCID PMC9181495 |
| Type | **Narrative review** — no new data, no datasets |
| Group | Lab. of Tumor Immunology and Cell Therapy, Sapienza University of Rome (Rughetti / Napoletano / Zizzari) |

## What it says (facts as stated by the review)

- Myeloid cells dominate the GBM microenvironment; TAMs (microglia + bone-marrow-derived
  macrophages, BMDM) are 30–50% of the tumour mass. Markers used to separate them:
  microglia CD45^low CD11b+ CX3CR1+ P2RY12+ CD49d− (TMEM119+); BMDM CD45^high CD49d+ P2RY12−.
- CD163+ "M2-like" TAMs are abundant in tumour core and peritumoral area; in the core they
  co-occur with IDO and PD-L1. MDSCs (M-MDSC, PMN-MDSC; LOX1+ PMN-MDSC) are immunosuppressive
  via ARG1, ROS, iNOS.
- Aberrant glycosylation in tumours: increased N-glycan branching, **truncated O-glycans**
  (Tn = GalNAcα1-O-Ser/Thr; sialyl-Tn), more sialylation (ST3GAL1 up in glioma, poor prognosis).
- Lectin families discussed: C-type lectins, Siglecs (5, 7, 9), galectins.

### MGL / CLEC10A specifically (section 4.1 and Table 1)

- **MGL (CLEC10A, CD301)** is a type II C-type lectin specific for the **Tn O-glycan**.
- Expression (Table 1): **DC and cDC2**, **CD163+ macrophages**, activated microglia.
- In EAE (mouse), MGL on macrophages/microglia induces apoptosis of autoreactive T cells and
  IL-10 release → MGL as a negative regulator of neuroinflammation (Ilarregui 2019).
- MGL engagement changes APC phenotype and metabolism; MGL+ APCs regulate naïve T cells and
  Tregs via Tn-glycosylated **CD45RA** (PTPRC isoform).
- In human GBM tissue, MGL and its Tn ligand are found, preferentially on tumour-infiltrating
  **CD163+** cells; in a Tn-high syngeneic mouse model, MGL+ areas had strong infiltration of
  **PD-L1+ TAMs** (Dusoswa 2020 — see [dusoswa2020_mgl_tn_gbm.md](dusoswa2020_mgl_tn_gbm.md)).
- Table 1 mechanisms: ERK phosphorylation and NF-κB; enhanced IL-10 secretion; TLR signalling;
  DC and TAM/CD163+-mediated immunosuppression.
- Table 1 MGL ligand carriers (from ovarian-cancer glycoproteomics, Napoletano 2020):
  matrix/adhesion VCAN, SDC3, PODXL, NID2, FN1, DAG1, APP, AGRN; metabolism ERP44, LAMP1,
  LAMP2, QSOX1, SEL1L, LRRC8D; mucins MUC1, MUC16, MUC24.
- Therapeutic precedent: the Tn–MGL axis was targeted with glycomimetic peptides in an ovarian
  cancer mouse model (Eggink 2018).

## Gaps the review leaves open (where single-cell data adds value)

1. **Which cells** express CLEC10A in human GBM, quantitatively: cDC2 vs CD163+ TAM vs microglia.
2. Whether CLEC10A+ cells, **within the same cell type**, carry an immunosuppressive programme
   (IL10, CD274/PD-L1, CD163, MRC1, TGFB1, …) versus an antigen-presenting one (MHC-II, CD1C).
3. Whether tumour cells express the **Tn-synthesis machinery** (GALNTs high, C1GALT1 /
   C1GALT1C1 low) in the same patients where CLEC10A+ myeloid cells are abundant.

## Caveats for our use

- MGL ligand carriers in Table 1 come from ovarian cancer; Tn status depends on
  glycosylation enzymes, not on carrier mRNA. Carrier expression is **not** a Tn proxy.
- The review is from the same group whose work we are extending: treat its framing
  (immunosuppression) as the hypothesis to test, not as a prior result.
