"""Checks on the pilot per-cell myeloid tables (results/pilot/myeloid_cells_<id>.tsv.gz):
table for docs, sanity markers, ADR requirements per donor, ASGR2 depth confound, donor-id prefixes.
Run: ~/.venvs/census/bin/python scripts/pilot/pilot_checks.py results/pilot/myeloid_cells_*.tsv.gz
Rank: ... pilot_checks.py --rank results/pilot/myeloid_cells_<tag>.tsv.gz [--genes scripts/pilot/<list>.txt]
      -> results/pilot/panel_rank.tsv (one row per gene; genes asked but absent from Census are listed)
"""
import sys, re
import pandas as pd
OBSC = ["dataset", "donor_id", "cell_type", "assay", "is_primary_data", "raw_sum"]
MAC, MIC = "macrophage", "microglial cell"

def rank(f, asked=None):
    c = pd.read_csv(f, sep="\t"); genes = [g for g in c.columns if g not in OBSC]
    mm = c[c.cell_type.isin([MAC, MIC])].copy()
    mm["terc"] = pd.qcut(mm.raw_sum, 3, labels=["low", "mid", "high"])
    big = c.groupby(["cell_type", "donor_id"]).size().ge(50)
    out = []
    for g in genes:
        pct = (c[g] > 0).groupby(c.cell_type).mean().mul(100)
        best = pct.idxmax()
        d = c[c.cell_type == best].groupby("donor_id")[g].apply(lambda x: 100 * (x > 0).mean())
        d = d[d.index.isin(big[best][big[best]].index)]
        t = (mm[g] > 0).groupby([mm.terc, mm.cell_type], observed=True).mean().unstack()
        sign = (t[MAC] - t[MIC]).apply(lambda v: 1 if v > 0 else -1)
        ratio = pct.get(MAC, 0) / max(pct.get(MIC, 0), 0.1)
        out.append(dict(gene=g, best_type=best, best_pct=round(pct.max(), 1),
                        donors_ge5_best=int((d >= 5).sum()), donors_best=len(d),
                        frac_donors_ge5=round((d >= 5).mean(), 2) if len(d) else 0.0,
                        mac_pct=round(pct.get(MAC, 0), 1), mic_pct=round(pct.get(MIC, 0), 1),
                        mac_mic_ratio=round(ratio, 2), depth_sign_consistent=f"{max((sign > 0).sum(), (sign < 0).sum())}/3",
                        dc_pct=round(pct.get("dendritic cell", 0), 1), mono_pct=round(pct.get("monocyte", 0), 1)))
    r = pd.DataFrame(out).sort_values(["frac_donors_ge5", "best_pct"], ascending=False)
    r.to_csv("results/pilot/panel_rank.tsv", sep="\t", index=False); print(r.to_string(index=False))
    if asked:
        miss = [g for g in asked if g not in genes]; print("not in Census var:", ", ".join(miss) or "none")

if len(sys.argv) > 1 and sys.argv[1] == "--rank":
    a = sys.argv[2:]; asked = None
    if "--genes" in a:
        i = a.index("--genes"); gf = a[i + 1]; a = a[:i] + a[i + 2:]
        asked = [x.split("#")[0].strip() for x in open(gf) if x.split("#")[0].strip()]
    rank(a[0], asked); sys.exit(0)

G = ["ASGR2", "CLEC10A", "SIGLEC9", "CLEC12A", "TMEM119", "P2RY12", "CD163"]
for f in sys.argv[1:]:
    c = pd.read_csv(f, sep="\t")
    if c.empty: print(f"\n##### {f}: 0 cells (dataset not in this Census release)"); continue
    c = c[c.is_primary_data.astype(str) == "True"] if "Extended" in f else c
    ds = c.dataset.iloc[0]; print(f"\n##### {ds}  ({len(c):,} myeloid cells, {c.donor_id.nunique()} donors)")
    pos = c[G] > 0
    t = pos.groupby(c.cell_type).mean().mul(100).round(1)
    nd = c.groupby(["cell_type", "donor_id"]).size().ge(50).groupby("cell_type").sum()
    t.insert(0, "donors_ge50", nd); t.insert(0, "n", c.cell_type.value_counts())
    print(t.to_string())
    mac, mic = "macrophage", "microglial cell"
    if mac in t.index and mic in t.index:
        ok = lambda g, hi, lo: f"{g} {t.at[hi, g]} vs {t.at[lo, g]} {'OK' if t.at[hi, g] > t.at[lo, g] else 'FAIL'}"
        print("sanity:", "; ".join([ok("SIGLEC9", mac, mic), ok("CLEC12A", mac, mic), ok("TMEM119", mic, mac), ok("P2RY12", mic, mac)]))
    # ADR requirements, per donor with >= 50 cells
    m = c[c.cell_type == mac]
    d = m.groupby("donor_id").filter(lambda x: len(x) >= 50).groupby("donor_id")
    a = d.ASGR2.apply(lambda x: 100 * (x > 0).mean())
    print(f"ADR ASGR2: donors with >=5% of macrophages >0: {(a >= 5).sum()}/{len(a)} (need >=3); median {a.median():.1f}%")
    cd = m[m.CD163 > 0].groupby("donor_id").filter(lambda x: len(x) >= 50).groupby("donor_id")
    b = cd.CLEC10A.apply(lambda x: 100 * (x > 0).mean())
    print(f"ADR CLEC10A: donors with >=2% of CD163+ macrophages >0: {(b >= 2).sum()}/{len(b)} (need >=3); median {b.median():.1f}%")
    # depth confound: ASGR2 % macrophage vs microglia within raw_sum terciles (cut on the pooled pair)
    mm = c[c.cell_type.isin([mac, mic])].copy()
    if mm.cell_type.nunique() == 2:
        mm["terc"] = pd.qcut(mm.raw_sum, 3, labels=["low", "mid", "high"])
        dep = (mm.assign(p=mm.ASGR2 > 0).groupby(["terc", "cell_type"], observed=True).p.mean().mul(100).round(1).unstack())
        dep["gap"] = (dep[mac] - dep[mic]).round(1)
        print("depth (ASGR2 % >0 by UMI tercile):", "; ".join(f"{k}: {r[mac]} vs {r[mic]} (gap {r.gap})" for k, r in dep.iterrows()))
        print("median UMI mac/mic:", int(mm[mm.cell_type == mac].raw_sum.median()), "/", int(mm[mm.cell_type == mic].raw_sum.median()))
    pre = c.donor_id.astype(str).str.replace(r"[-_ ]?\d+[A-Za-z]?$", "", regex=True).value_counts()
    print("donor-id prefixes (top 15):", ", ".join(f"{k}({v})" for k, v in pre.head(15).items()))
