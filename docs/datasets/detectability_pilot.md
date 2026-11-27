# Detectability pilot — lectin genes in GBM myeloid cells (CELLxGENE)

Status 2026-11-27 session 4 (laptop): **MEASURED on Core GBmap** (Census `2025-11-08`, raw counts, 14 genes). Glioma cohort atlas `8f3cb9c4`: **0 cells, not in Census 2025-11-08** → needs h5ad or a newer release. Extended GBmap, SL040/SL057, Darmanis: not run (cannot change the selection: Extended contains Core; others small). Run: 95 s, 10.8 GB peak RAM per dataset; `soma.init_buffer_bytes=2**28` stalled > 20 min (do not use).

## What was verified
- GBmap is on CELLxGENE: collection `999f2a15…`, **Core GBmap** `c888b684-6c51-431f-972a-6c963044cef0` (338,564 cells), **Extended GBmap** `56c4912d-2bae-4b64-98f2-af8a84389208` (1,135,677 cells); assays 10x 3' v2/v3, 10x 5' v1, BD Rhapsody (curation API `/curation/v1/datasets`).
- Other human glioma datasets on CELLxGENE: glioma cohort atlas `8f3cb9c4…` (638,930 cells, 66 samples, IDH-wt GBM + IDH-mut, 10x 3' v3, doi 10.1002/mco2.71022); SL040 `21722308…` and SL057 `145dcf6a…` (GBM, 10x v3, bioRxiv 2023.09.01.555882); Darmanis 2017 `774c18c5…` (Smart-seq2, 3,589 cells); pediatric HGG snRNA (out of scope).
- Ensembl ids of the 14 genes resolved via WMG `primary_filter_dimensions`.

## What failed
- Laptop: no `pip`, no `ensurepip` (`python3.10-venv` not installed) → no `~/.venvs/census`. Bootstrapping pip from a downloaded script was blocked by the permission classifier. Fix (simulated safe: 3 packages, 0 removals, 0 upgrades): `sudo apt install --no-install-recommends python3.10-venv`, then `python3 -m venv ~/.venvs/census && ~/.venvs/census/bin/pip install cellxgene-census`.
- Fallback WMG API (`/wmg/v2/query`): dataset-filtered queries (Core GBmap, glioma atlas) fail with a server-side schema error (empty result); disease-filtered (MONDO:0018177) queries return HTTP 500, with and without `compare`.

## Ready to run
`scripts/pilot/census_detectability.py` (obs + 14 genes, raw counts): per dataset × cell type × donor → n cells, % > 0, mean in positives; myeloid summary with donors ≥ 50 cells and donors with ≥ 5% positive; prints the stop-rule check (ASGR2 < 1% everywhere).

## Caveats to apply when reading the output
- Census `cell_type` is the CL term only (macrophage / microglial cell / monocyte / dendritic cell); GBmap's finer labels (Mo-TAM vs Mg-TAM) are not Census columns → BMD vs MG = CL macrophage vs microglial cell, checked with sanity markers.
- Sanity checks (session-3 prompt): SIGLEC9 and CLEC12A higher in macrophage than microglia (E14, E20); TMEM119/P2RY12 higher in microglia. If they fail, annotation is not trusted.
- Skeptic point: % > 0 tracks UMI depth; compare BMD vs MG also after depth matching (needs `nCount` → add `raw_sum` to obs columns).
- GBmap pools many studies: source study is not a Census column; donor_id prefix may identify it `[UNVERIFIED]`.

| dataset | cell type | n cells | donors ≥ 50 cells | ASGR2 % > 0 | CLEC10A % > 0 | SIGLEC9 | CLEC12A | TMEM119 |
|---|---|---|---|---|---|---|---|---|
| Core GBmap | macrophage | 61,134 | 65 | 3.3 | 6.2 | 24.7 | 12.1 | 20.9 (P2RY12 8.6) |
| Core GBmap | microglial cell | 54,722 | 69 | 1.5 | 0.9 | 21.3 | 5.3 | 40.5 (P2RY12 24.8) |
| Core GBmap | monocyte | 14,215 | 26 | 9.9 | 15.5 | 18.5 | 30.8 | 3.1 |
| Core GBmap | dendritic cell | 3,961 | 14 | 9.7 | 35.7 | 9.5 | 23.5 | 2.2 |
| Glioma cohort atlas | — | 0 (not in Census 2025-11-08) | | | | | | |

Source: `results/pilot/myeloid_cells_c888b684.tsv.gz` + `scripts/pilot/pilot_checks.py` (103 donors, 134,032 myeloid cells).

## Results (Core GBmap)
- **Sanity checks: all pass** → CL annotation trusted for macrophage vs microglia: SIGLEC9 24.7 vs 21.3, CLEC12A 12.1 vs 5.3 (macrophage higher); TMEM119 40.5 vs 20.9, P2RY12 24.8 vs 8.6 (microglia higher). SIGLEC9 gap small.
- **ADR ASGR2 requirement met:** 8/65 donors have ≥ 5% ASGR2+ macrophages (need ≥ 3); median donor 1.8% → low, sparse signal; highest in monocytes (9.9%) and DC (9.7%), not macrophages.
- **ADR CLEC10A requirement met:** 25/55 donors have ≥ 2% CLEC10A+ among CD163+ macrophages; median 1.5%. CLEC10A is mainly DC (35.7%) and monocytes (15.5%), consistent with E13 (cDC2).
- **Depth confound (skeptic):** ASGR2 % macrophage vs microglia by UMI tercile: low 0.7 vs 0.5, mid 2.5 vs 1.1, high 5.7 vs 3.6 → the gap persists within terciles (not explained by depth), but detection is depth-limited (< 1% in the low tercile). Median UMI 6,690 vs 5,008.
- **Stop rule not hit** (max 9.9% in a myeloid type, 3.3% in macrophages).
- **Overlap:** Core donor-id prefixes include `R` (28,857 cells) and `ND` (22,375) = GSE163120 naming (R1–5, ND1–8) → GSE163120 RNA is likely already inside Core GBmap `[prefix match, UNVERIFIED]`; also `MGH` (Neftel), `BT`, `SF`, `G`, `SM`.
