# Detectability pilot — lectin genes in GBM myeloid cells (CELLxGENE)

Status 2026-11-27 (laptop): **NOT MEASURED. Blocked on environment.** No numbers below are measured yet.

## What was verified
- GBmap is on CELLxGENE: collection `999f2a15…`, **Core GBmap** `c888b684-6c51-431f-972a-6c963044cef0` (338,564 cells), **Extended GBmap** `56c4912d-2bae-4b64-98f2-af8a84389208` (1,135,677 cells); assays 10x 3' v2/v3, 10x 5' v1, BD Rhapsody (curation API `/curation/v1/datasets`).
- Other human glioma datasets on CELLxGENE: glioma cohort atlas `8f3cb9c4…` (638,930 cells, 66 samples, IDH-wt GBM + IDH-mut, 10x 3' v3, doi 10.1002/mco2.71022); SL040 `21722308…` and SL057 `145dcf6a…` (GBM, 10x v3, bioRxiv 2023.09.01.555882); Darmanis 2017 `774c18c5…` (Smart-seq2, 3,589 cells); pediatric HGG snRNA (out of scope).
- Ensembl ids of the 14 genes resolved via WMG `primary_filter_dimensions`.

## What failed
- Laptop: no `pip`, no `ensurepip` (`python3.10-venv` not installed) → no `~/.venvs/census`. Bootstrapping pip from a downloaded script was blocked by the permission classifier. Fix: `sudo apt install python3.10-venv python3-pip`, then `python3 -m venv ~/.venvs/census && ~/.venvs/census/bin/pip install cellxgene-census`.
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
| (to fill from `results/pilot/detectability_myeloid_agg.tsv`) | | | | | | | | |
