# Dusoswa et al. 2020 — GBM exploits truncated O-glycans via MGL

`[VERIFIED 2026-11-27 against the full text: publisher's version, CC BY 4.0, Leiden repository
hdl.handle.net/1887/3626991. Europe PMC fullTextXML returns HTTP 500 for this PMCID. SI Appendix not read.]`

| | |
|---|---|
| Citation | Dusoswa SA, Verhoeff J, Abels E, Méndez-Huergo SP, Croci DO, Kuijper LH, de Miguel E, Wouters VMCJ, Best MG, Rodriguez E, et al. *Glioblastomas exploit truncated O-linked glycans for local and distant immune modulation via the macrophage galactose-type lectin.* PNAS. 2020;117(7):3693–3703 |
| IDs | DOI 10.1073/pnas.1907921117 · PMID 32019882 · PMCID PMC7035608 |
| Why it matters | The **only** primary study on MGL in GBM; cited by both Pace papers as the MGL evidence |

## Findings, by type of evidence

**Human tissue (glycan and protein level)**
- Tn glycans (HPA lectin staining, cryosections) higher in GBM than in lower-grade glioma and
  epilepsy tissue; images representative of **5 patients per condition** (Fig. 1E–F).
- MGL-binding ELISA on tissue lysates: GBM n = 12, LGG n = 5, epilepsy n = 8 (Fig. 1G).
- Immunofluorescence: GBM tissue highly infiltrated by **MGL+ cells that are also CD163+** (Fig. 2A–B).

**Human bulk mRNA (reanalysis of Gravendeel et al., microarray via the R2 platform)**
- Healthy brain n = 8, astrocytoma WHO II n = 13, WHO III n = 60, GBM n = 159.
- **GALNT1, GALNT2, GALNT7, GALNT10** up-regulated in GBM, from a volcano plot of 270
  glycosylation genes (Fig. 1D).
- MGL and CD163 mRNA higher in GBM; MGL vs CD163 in GBM: Pearson R = 0.29, P ≤ 0.0001 (Fig. 2C–E).
- Survival: low MGL (n = 84) vs high MGL (n = 85), median cutoff 3.25, P = 0.032 (Fig. 2F).

**Mouse (not transferable to human MGL)**
- GL261 with *C1galt1c1* (Cosmc) knocked out → Tn-high tumours. In vitro, MGL-ligand-high
  glioma cells induce tolerogenic markers (incl. PD-L1) on macrophages and IL-10 on microglia (SI).
- Mass cytometry (37 markers) on tumour, spleen and bone marrow of 20 mice (10 vs 10):
  more **CCR2+ PD-L1+ macrophages** in Tn-high tumours, plus distant bone-marrow changes.
- The authors name mouse **MGL2** as the homolog of human MGL by glycan specificity; no
  *Mgl2*−/− experiment was done.

## Limits that matter for us

- The MGL–CD163 link in humans rests on IF images and a bulk correlation (R = 0.29); there is
  **no single-cell evidence** of which myeloid cells carry MGL.
- The PD-L1+ TAM result is mouse-only.
- GALNT up-regulation is bulk mRNA: it does not tell which cells make Tn, and mRNA ≠ glycan.

## What our single-cell data can test

- Which myeloid populations express CLEC10A (cDC2 vs CD163+ TAMs vs microglia), per patient.
- CLEC10A+ vs CLEC10A− TAMs for CD274 (PD-L1), CD163, IL10.
- At patient level: tumour-cell Tn-synthesis genes (GALNT1/2/7/10, C1GALT1, C1GALT1C1)
  against the abundance of CLEC10A+ TAMs. Gene expression is only a **proxy** for Tn display.
- GSE163120 CITE-seq has **CD301** (MGL protein): mRNA vs protein in the same cell.
- Use human datasets only (mouse has Mgl1 and Mgl2 with different specificities).
