# Motevasseli et al. 2024 — tumour–TAM interactions in IDH-stratified glioma (sc + spatial)

`[VERIFIED 2026-11-27 against the PMC full text, sentence-level grep for hub lectins]`

| | |
|---|---|
| Citation | Motevasseli M, et al. *Distinct tumor-TAM interactions in IDH-stratified glioma microenvironments unveiled by single-cell and spatial transcriptomics.* Acta Neuropathol Commun. 2024 |
| IDs | DOI 10.1186/s40478-024-01837-5 · PMID 39148129 · PMCID PMC11328419 |
| Data | Re-analysis of public data: Johnson et al. 2021 human glioma scRNA-seq (Synapse `singlecellglioma`), Ravi et al. spatial (Scanpy) |

## Facts reported
- CD14+ myeloid cluster, 12,167 cells. MG-TAM = TMEM119, P2RY12, CX3CR1; classical monocytes = SELL, FCN1, VCAN; **BMD-TAM = TGFBI, CLEC12A** (Fig. 1B).
- The only hub lectin mentioned is CLEC12A, and only as a marker. ASGR2 and CLEC10A are absent from the text.

## Related check (Müller 2017, the usual source for BMD vs MG markers)
Müller S et al. Genome Biol 2017, PMID 29262845, PMC5738907, *Single-cell profiling of human gliomas reveals macrophage ontogeny…*: the main text names TGFBI and BIN1, but **not CLEC12A, ASGR2 or CLEC10A** (supplementary Table S3 not checked).

## What our sc data can test
- CLEC12A as a BMD-TAM marker is **already known**, so "who carries CLEC12A" is not new on its own. What is new: whether ASGR2 and CLEC10A sit in the same BMD-TAM cells as CLEC12A.
- Johnson 2021 (Synapse) is a further IDH-stratified sc resource for Phase 3.
