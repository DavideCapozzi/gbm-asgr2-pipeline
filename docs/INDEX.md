# docs/ — research context for the MGL (CLEC10A) extension

Reference material for anyone (human or agent) extending the pipeline to new datasets.
Each file has a single purpose; read only the one you need.

| File | Purpose | Read it when… |
|---|---|---|
| [literature/pace2022_glycan_lectin_review.md](literature/pace2022_glycan_lectin_review.md) | The lab's 2022 review: the MGL–Tn axis and lectins in GBM immunosuppression | framing hypotheses, choosing marker/ligand genes |
| [literature/pace2026_asgr2_clec12a_hubs.md](literature/pace2026_asgr2_clec12a_hubs.md) | The lab's 2026 paper: ASGR2 and CLEC12A as hubs of the lectin co-expression network (bulk + flow) | linking sc results to the published lectin network |
| [literature/dusoswa2020_mgl_tn_gbm.md](literature/dusoswa2020_mgl_tn_gbm.md) | Main external MGL-in-GBM paper (Tn glycans → MGL on CD163+ TAMs → PD-L1+ TAMs) | defining the immunosuppression hypothesis |
| [literature/motevasseli2024_idh_tam_sc_clec12a.md](literature/motevasseli2024_idh_tam_sc_clec12a.md) | Human glioma sc: CLEC12A is a known BMD-TAM marker (with TGFBI); Müller 2017 negative check | judging novelty of CLEC12A in sc |
| [literature/mei2023_siglec9_tam_gbm.md](literature/mei2023_siglec9_tam_gbm.md) | SIGLEC9-high monocyte-derived TAMs in human GBM sc (abstract only) | Siglec positive control, novelty of storyline (a) |
| [literature/schmassmann2023_siglec_sialic_gbm.md](literature/schmassmann2023_siglec_sialic_gbm.md) | Siglec-9/E–sialic acid axis, survival + mouse KO (abstract only) | Siglec positive control |
| [literature/santegoets2019_siglec_mdsc_glioma.md](literature/santegoets2019_siglec_mdsc_glioma.md) | Siglec-3/5/7/9 on glioma MDSCs; glioma cells carry Siglec-7/9 ligands (flow) | Siglec ligand–receptor framing |
| [literature/cornelissen2024_sialoglycan_myeloid_gbm.md](literature/cornelissen2024_sialoglycan_myeloid_gbm.md) | Desialylation of GBM cells lowers CD163 induction on monocytes (in vitro) | storyline (b), Siglec version |
| [literature/trindade2026_transmembrane_lectins_review.md](literature/trindade2026_transmembrane_lectins_review.md) | 2026 review: CLEC12A as ITIM receptor, MGL as modulatory; no sc ASGR2 data cited | Discussion framing |
| [literature/napoletano2020_mgl_glycoproteome_ovarian.md](literature/napoletano2020_mgl_glycoproteome_ovarian.md) | Lab paper: MGL ligand carriers in ovarian cancer O-glycoproteome | storyline (b) carrier genes |
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
- `docs/decisions/0001-storyline.md` — ADR: storyline (a) chosen, (c) reserve, dataset requirements
- `docs/datasets/detectability_pilot.md` — CELLxGENE pilot status (GBmap ids; not measured yet)
- `scripts/pilot/census_detectability.py` — Census 14-gene detectability pilot
