"""Detectability pilot: % cells > 0 and mean in positives for 14 lectin/marker genes,
per source dataset x cell type x donor, from CELLxGENE Census (obs + 14 genes only).

Run:  ~/.venvs/census/bin/python scripts/pilot/census_detectability.py [census_version]
Out:  results/pilot/detectability_long.tsv, results/pilot/detectability_myeloid_wide.tsv
Needs: pip install cellxgene-census (venv needs python3.10-venv: `sudo apt install python3.10-venv`).
Never loads the full matrix: var filter keeps 14 genes (~1.1 M cells x 14 sparse ints).
"""
import sys, pathlib
import numpy as np, pandas as pd
import cellxgene_census

GENES = ["ASGR2", "CLEC12A", "CLEC10A", "SIGLEC9", "SIGLEC7", "CD163", "TGFBI",
         "P2RY12", "TMEM119", "ITGA4", "CD1C", "FCER1A", "C1QA", "PTPRC"]
# Verified on CELLxGENE curation API 2026-11-27.
DATASETS = {
    "c888b684-6c51-431f-972a-6c963044cef0": "Core GBmap",
    "56c4912d-2bae-4b64-98f2-af8a84389208": "Extended GBmap",
    "8f3cb9c4-38d5-4b9c-b309-bcccd76d0075": "Glioma cohort atlas (66 samples, MedComm 2026)",
    "21722308-5091-4b63-9c07-4f116ab9a7b3": "SL040 (bioRxiv 2023.09.01.555882)",
    "145dcf6a-2461-4fa3-a0af-7fc56db0bd33": "SL057 (bioRxiv 2023.09.01.555882)",
    "774c18c5-efa1-4dc5-9e5e-2c824bab2e34": "Darmanis 2017 (Smart-seq2)",
}
MYELOID = ("macrophage", "microglial cell", "monocyte", "dendritic cell", "myeloid",
           "mononuclear phagocyte", "neutrophil", "conventional dendritic cell")
OBS = ["dataset_id", "donor_id", "cell_type", "disease", "assay", "tissue", "is_primary_data", "raw_sum"]  # raw_sum = UMI depth (skeptic: depth confound)

def fetch(census, dsid, name):
    """One dataset at a time (Core GBmap alone peaked at 10.7 GB RSS in a single query)."""
    ad = cellxgene_census.get_anndata(
        census, organism="Homo sapiens", X_name="raw",
        obs_value_filter=f"dataset_id == '{dsid}'",
        var_value_filter=f"feature_name in {GENES}", obs_column_names=OBS)
    X = ad.X.tocsc(); genes = ad.var["feature_name"].tolist()
    obs = ad.obs.reset_index(drop=True); obs["dataset"] = name
    rows = []
    for key, idx in obs.groupby(["dataset", "cell_type", "donor_id"], observed=True).indices.items():
        sub = X[idx]
        for j, g in enumerate(genes):
            col = sub[:, j].toarray().ravel(); pos = col > 0
            rows.append((*key, g, len(idx), 100 * pos.mean(), col[pos].mean() if pos.any() else 0.0))
    # per-cell myeloid table (obs + 14 genes) for depth-confound and CD163+ checks
    m = obs.cell_type.astype(str).str.contains("|".join(MYELOID), case=False).values
    cells = obs.loc[m, ["dataset", "donor_id", "cell_type", "assay", "is_primary_data", "raw_sum"]].copy()
    cells[genes] = X.tocsr()[np.where(m)[0]].toarray()
    print(f"{name}: {len(obs):,} cells, {m.sum():,} myeloid", flush=True)
    return rows, cells

def main(version="2025-11-08"):
    import os, tiledbsoma
    out = pathlib.Path("results/pilot"); out.mkdir(parents=True, exist_ok=True)
    only = os.environ.get("PILOT_ONLY")  # comma-separated dataset_id prefixes, for a quick test
    ds = {k: v for k, v in DATASETS.items() if not only or any(k.startswith(o) for o in only.split(","))}
    # PILOT_BUFFER (bytes) shrinks TileDB read buffers; 2**28 stalled remote reads >20 min (2026-11-27). Default off.
    buf = os.environ.get("PILOT_BUFFER")
    ctx = tiledbsoma.SOMATileDBContext(tiledb_config={"soma.init_buffer_bytes": int(buf)}) if buf else None
    rows, cells = [], []
    with cellxgene_census.open_soma(census_version=version, **({"context": ctx} if ctx else {})) as census:
        for dsid, name in ds.items():
            r, c = fetch(census, dsid, name); rows += r; cells.append(c)
    pd.concat(cells).to_csv(out / f"myeloid_cells{'_' + only.replace(',', '_') if only else ''}.tsv.gz", sep="\t", index=False)
    long = pd.DataFrame(rows, columns=["dataset", "cell_type", "donor_id", "gene", "n_cells", "pct_gt0", "mean_in_pos"])
    tag = "_" + only.replace(",", "_") if only else ""
    long.to_csv(out / f"detectability_long{tag}.tsv", sep="\t", index=False)
    mye = long[long.cell_type.str.contains("|".join(MYELOID), case=False)]
    # per dataset x cell type: pooled %, n donors with >=50 cells, n donors with pct>=5
    agg = (mye.assign(w=mye.n_cells * mye.pct_gt0 / 100)
              .groupby(["dataset", "cell_type", "gene"])[["n_cells", "pct_gt0", "w"]]
              .apply(lambda d: pd.Series({
                  "n_cells": d.n_cells.sum(),
                  "pct_gt0": 100 * d.w.sum() / d.n_cells.sum(),
                  "donors_ge50": (d.n_cells >= 50).sum(),
                  "donors_pct_ge5": ((d.n_cells >= 50) & (d.pct_gt0 >= 5)).sum()})))
    wide = agg["pct_gt0"].unstack("gene").round(2)
    wide.to_csv(out / f"detectability_myeloid_wide{tag}.tsv", sep="\t")
    agg.to_csv(out / f"detectability_myeloid_agg{tag}.tsv", sep="\t")
    print(wide.to_string())
    # stop rule (session 3): ASGR2 < 1% in all datasets' macrophages
    asgr2 = agg.xs("ASGR2", level="gene")
    print("\nASGR2 max pooled % in myeloid types:", asgr2.pct_gt0.max().round(2),
          "-> STOP RULE HIT" if asgr2.pct_gt0.max() < 1 else "-> continue")

if __name__ == "__main__":
    main(*sys.argv[1:2])
