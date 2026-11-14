"""step14_bench_score.py - score the friend's two EvSpark benchmark files with the kept metrics.
  fulllength_bench_L40S.json : phiX174 (Microviridae, in-domain), 3 seeds x {plain Evo 2 ('nati'), EvSpark ('spec')}, 64-letter prompt (the first 64 letters of NC_001422.1) + 5,322 generated letters = 5,386.
  longrun_bench_L40S.json    : E. coli K-12 region, 3 seeds x {plain, EvSpark}, 8,192 generated letters, plus the real region; OUT OF DOMAIN for phage-calibrated metrics.
Scores use exactly the step4 definitions (higher = more natural). 'Calibration range' = 5th-95th percentile of the 30 calibration genomes.
Important: the calibration genomes start at the F2 gene; phiX174 and the benchmark sequences start at the NCBI origin, which lies inside gene A, so a linear gene finder
cuts a gene at the ends. The honest baseline is therefore the real phiX174 genome in the SAME frame (row 'phiX174_natural'), shown next to the calibration range.
Inputs: bench/full/*.csv, bench/long/*.csv (metric tables from the scripts), bench/checkv_full/quality_summary.tsv, optional bench/phold/phold_per_cds_predictions.tsv."""
import os, numpy as np, pandas as pd
M = "data/metrics_"
cal_tabs = [pd.read_csv(M + f + ".csv") for f in ("quick", "phanotate", "pfam", "phrog", "checkv")]
cal = cal_tabs[0]
for t in cal_tabs[1:]:
    cal = cal.merge(t, on="id")
cal = cal[cal.id.str.startswith("CAL|")]
def z(col, c):
    iqr = c.quantile(.75) - c.quantile(.25); sd = iqr if iqr > 0 else (c.std() or 1.0)
    return (col - c.median()) / sd
def scores(q):
    S = pd.DataFrame({"id": q.id})
    zmax = lambda cols: -pd.concat([z(q[c], cal[c]).abs() for c in cols], axis=1).max(axis=1)
    S["M1"] = zmax(["length", "m1_gc", "m1_maxrun"])
    if "m2_completeness" in q: S["M2"] = q.m2_completeness
    S["M4"] = zmax(["m4_genes", "m4_density", "m4_meanlen"]); S["M5"] = zmax(["m5_genes", "m5_density", "m5_meanlen"])
    S["M6"] = q.m6_core; S["M8"] = q.m8_order; S["M9"] = (q.m9_gap + q.m9_codon) / 2; S["M10"] = q.m10_rbs
    S["M11"] = q.m11_pfam_share
    if "m12_phrog_share" in q: S["M12"] = q.m12_phrog_share
    if "m13_share" in q: S["M13"] = q.m13_share
    return S
calS = scores(cal)
if os.path.exists("data/metrics_phold.csv"): pass
rng = {m: (calS[m].quantile(.05), calS[m].quantile(.95)) for m in calS.columns if m != "id"}
# M13 calibration range comes from Phold run on the 30 calibration genomes (phold/out2, ids start with CAL_)
_p2 = pd.read_csv("phold/out2/phold_per_cds_predictions.tsv", sep="	"); _p2["l"] = _p2.phrog != "No_PHROG"
_c = _p2[_p2.contig_id.str.startswith("CAL_")].groupby("contig_id").l.mean()
rng["M13"] = (_c.quantile(.05), _c.quantile(.95))
def load(kind):
    d = f"bench/{kind}"
    t = [pd.read_csv(f"{d}/metrics_{f}.csv") for f in ("quick", "phanotate", "pfam") if os.path.exists(f"{d}/metrics_{f}.csv")]
    q = t[0]
    for x in t[1:]: q = q.merge(x, on="id")
    if os.path.exists(f"{d}/metrics_phrog.csv"): q = q.merge(pd.read_csv(f"{d}/metrics_phrog.csv"), on="id")
    q = q[~q.id.str.startswith("CAL|")]
    return q
pd.set_option("display.width", 250)
out = []
for kind, label in (("full", "phiX174 full-length (in-domain)"), ("long", "E. coli 8 kb (out of domain)")):
    q = load(kind)
    ck = f"bench/checkv_{kind if kind == 'long' else 'full'}/quality_summary.tsv"
    if os.path.exists(ck):
        c = pd.read_csv(ck, sep="	"); q = q.merge(pd.DataFrame({"id": c.contig_id, "m2_completeness": pd.to_numeric(c.completeness, errors="coerce").fillna(0)}), on="id")
    ph = "bench/phold/phold_per_cds_predictions.tsv"
    if os.path.exists(ph):
        p = pd.read_csv(ph, sep="	"); p["l"] = p.phrog != "No_PHROG"
        s13 = p.groupby("contig_id").l.mean().rename("m13_share").reset_index().rename(columns={"contig_id": "id"})
        s13 = s13[s13.id.isin(q.id)]; q = q.merge(s13, on="id")
    S = scores(q); S["group"] = S.id.str.replace(r"_s\d$", "", regex=True)
    cols = [c for c in S.columns if c.startswith("M")]
    print(f"\n=== {label}: per-sequence scores (higher = more natural) ===")
    print(S.round(2).to_string(index=False))
    print("\n  calibration range (5th-95th):", {m: (round(r[0], 2), round(r[1], 2)) for m, r in rng.items() if m in cols})
    base = S[S.group.str.contains("natural|real")].iloc[0]
    print("\n  group means (3 seeds) and how many of 3 fall inside the calibration range; real baseline in the first row")
    rows = []
    for g, d in S.groupby("group"):
        r = {"group": g}
        for m in cols:
            ins = ((d[m] >= rng[m][0]) & (d[m] <= rng[m][1])).sum() if m in rng else 0
            r[m] = f"{d[m].mean():.2f} ({int(ins)}/{len(d)})"
        rows.append(r)
    print(pd.DataFrame(rows).to_string(index=False))
    S.to_csv(f"results/bench_{kind}_scores.csv", index=False)
