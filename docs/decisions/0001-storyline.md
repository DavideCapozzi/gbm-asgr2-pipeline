# ADR 0001 — Storyline: lectin hubs (ASGR2, CLEC10A) at single-cell resolution

Date: 2026-11-27 · Status: accepted, **conditional on the detectability pilot** (not yet run)

## Context
Three storylines were written and scored 1–5 (novelty vs lab, testability, feasibility, risk) with evidence ids from the project evidence log (E1–E24):
(a) lectin hubs at sc resolution = 14; (b) tumour glycosylation ↔ lectin+ TAM per patient = 11; (c) CITE-seq protein validation = 14.
A blind reviewer (no scores) ranked (a) > (c) > (b).

## Decision
- **Chosen: (a).** ASGR2 and CLEC10A have no published sc description in GBM (E17, E18), while SIGLEC9 (E20) and CLEC12A (E14, E16) give built-in positive controls. Needs only public data and per-patient descriptive statistics.
- **Reserve: (c)**, also usable as the protein module of (a) for CLEC10A/CLEC12A (CD301/CD371 ADT; E11, E12). Cannot test ASGR2.
- **Dropped for now: (b)**: mRNA ≠ glycan, ≥ 15 unsorted IDH-wt patients needed, confounded by purity and myeloid abundance.

## Falsifiers for (a)
1. Depth-matched BMD vs MG ASGR2 effect has a CI including 0 across ≥ 5 patients.
2. ASGR2/CLEC12A same-cell co-occurrence ≈ expected after conditioning on BMD label and UMI depth.
3. Stop rule: ASGR2 > 0 in < 1% of macrophages in every dataset → storyline re-scoped (ASGR2 = stated mRNA limit, focus CLEC10A + CITE-seq).

## Dataset requirements for (a)
| requirement | threshold | why |
|---|---|---|
| BMD-TAM vs MG-TAM separable | TGFBI/CLEC12A/ITGA4 vs P2RY12/TMEM119 give two clusters in each patient | core contrast (E4, E14) |
| cDC2 and monocytes present | ≥ 30 cells each (CD1C/FCER1A; CD14/FCN1) | CLEC10A cDC2 (E13); lineage control for ASGR2 |
| ASGR2 detectable | > 0 in ≥ 5% of BMD-TAM in ≥ 3 patients | below this, % > 0 is noise |
| CLEC10A detectable | > 0 in ≥ 2% of CD163+ TAM in ≥ 3 patients | E7 at mRNA level |
| patients | ≥ 5 IDH-wt, ≥ 200 myeloid cells each | per-patient paired test (skeptic) |
| tumour type | IDH-wt; primary and recurrent analysed apart | biology differs |
| chemistry | 10x v3 preferred; v2 / Smart-seq2 as separate strata | detection sensitivity |
| depth | per-cell UMI available | depth matching |
| protein (optional) | CITE-seq with CD301, CD371, CD49d, TMEM119 | reserve (c) |

## Consequences
Dataset selection (session 3, step 3) is provisional until `scripts/pilot/census_detectability.py` runs (laptop needs `python3.10-venv`). If the stop rule fires, revise this ADR before any download.
