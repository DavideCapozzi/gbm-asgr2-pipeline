# docs/ — research context for the MGL (CLEC10A) extension

Reference material for anyone (human or agent) extending the pipeline to new datasets.
Each file has a single purpose; read only the one you need.

| File | Purpose | Read it when… |
|---|---|---|
| [literature/pace2022_glycan_lectin_review.md](literature/pace2022_glycan_lectin_review.md) | The lab's 2022 review: the MGL–Tn axis and lectins in GBM immunosuppression | framing hypotheses, choosing marker/ligand genes |
| [literature/pace2026_asgr2_clec12a_hubs.md](literature/pace2026_asgr2_clec12a_hubs.md) | The lab's 2026 paper: ASGR2 and CLEC12A as hubs of the lectin co-expression network (bulk + flow) | linking sc results to the published lectin network |
| [literature/dusoswa2020_mgl_tn_gbm.md](literature/dusoswa2020_mgl_tn_gbm.md) | Main external MGL-in-GBM paper (Tn glycans → MGL on CD163+ TAMs → PD-L1+ TAMs) | defining the immunosuppression hypothesis |
| [datasets/candidate_datasets.md](datasets/candidate_datasets.md) | Public GBM scRNA-seq datasets screened for replication, with format, size and priority | adding a dataset to `scripts/main.R` |

## Conventions used in these files

- **Status tags.** `[VERIFIED <date>]` = checked against the primary source or the actual
  file on that date. `[UNVERIFIED]` = from a search summary or memory; check before relying on it.
- **Facts vs hypotheses.** Each literature note keeps *what the paper shows* separate from
  *what our single-cell data can test*. Never move a hypothesis into the facts section.
- **Gene → claim labels** must be checked gene by gene (a wrong label in a figure is a
  falsifiable error). Approximate members go to "Other".
- Dates are absolute (YYYY-MM-DD). Accessions are given exactly as on GEO/PMC.
- The execution plan and phase roadmap are local-only (gitignored) files in the repo root:
  `PLAN_clec10a_consensus.md`, `ROADMAP_clec10a_extension.md`.
