# Dusoswa et al. 2020 — GBM exploits truncated O-glycans via MGL

`[ABSTRACT ONLY — full text not retrieved (Europe PMC full text unavailable, 2026-11-23). Details marked UNVERIFIED come from a search summary.]`

| | |
|---|---|
| Citation | Dusoswa SA, Verhoeff J, Abels E, Méndez-Huergo SP, Croci DO, Kuijper LH, de Miguel E, Wouters VMCJ, Best MG, Rodriguez E, et al. *Glioblastomas exploit truncated O-linked glycans for local and distant immune modulation via the macrophage galactose-type lectin.* PNAS. 2020;117(7):3693–3703 |
| IDs | DOI 10.1073/pnas.1907921117 · PMID 32019882 · PMCID PMC7035608 |
| Why it matters | The **only** primary study on MGL in GBM; cited by both Pace papers as the MGL evidence |

## Findings (from the abstract)

- Human GBM tissue: increased **truncated O-linked glycans (Tn)** and of their receptor **MGL**,
  on **CD163+ TAMs**.
- Orthotopic immunocompetent mouse glioma overexpressing MGL ligands: mass cytometry shows
  heterogeneous myeloid infiltration, more **PD-L1+ TAMs**, and distant bone-marrow changes.
- Conclusion: GBM uses cell-surface O-glycans for local and distant immune modulation.

`[UNVERIFIED]` Truncated glycans higher in GBM than lower-grade glioma and epilepsy tissue;
GALNT1, GALNT2, GALNT7, GALNT10 up-regulated in GBM.

## Implications for our analysis

- The paper places MGL on **CD163+ TAMs**, so the immunosuppression hypothesis should be
  tested **in the TAM compartment**, not only in cDC2 (where CLEC10A is a lineage marker).
- Testable in sc data: CLEC10A+ vs CLEC10A− TAMs for CD274 (PD-L1), CD163, IL10;
  and, at the patient level, tumour-cell Tn-synthesis genes (GALNT1/2/7/10 vs C1GALT1,
  C1GALT1C1) against the abundance of CLEC10A+ TAMs. Gene expression is only a **proxy**
  for Tn glycan display.
- Mouse data is not a substitute: mouse has two Mgl genes (Mgl1/Clec10a, Mgl2) whose
  ligand specificity differs from human MGL. Use human datasets only.
