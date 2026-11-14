"""step16_kmer_arm_score.py - score the k-mer arm against plain Evo 2 on the SAME 23 prompts, with the kept metrics.
Input : kmerarm/work/data/metrics_{quick,phanotate,pfam,phrog}.csv  (same step3 scripts, one run, ids ACC__natural / ACC__BASE / ACC__KMER)
        kmerarm/checkv_out/quality_summary.tsv (M2), data/metrics_*.csv (30 calibration genomes -> calibration ranges)
Scores use the step14/step4 definitions (higher = more natural). 'Inside' = within the 5th-95th percentile of the 30 calibration genomes.
M13 (Phold) is not run here (needs Modal, costs money)."""
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
M = "data/metrics_"
cal_tabs = [pd.read_csv(M + f + ".csv") for f in ("quick", "phanotate", "pfam", "phrog", "checkv")]
cal = cal_tabs[0]
for t in cal_tabs[1:]: cal = cal.merge(t, on="id")
cal = cal[cal.id.str.startswith("CAL|")]
def z(col, c):
    iqr = c.quantile(.75) - c.quantile(.25); sd = iqr if iqr > 0 else (c.std() or 1.0)
    return (col - c.median()) / sd
def scores(q):
    S = pd.DataFrame({"id": q.id})
    zmax = lambda cols: -pd.concat([z(q[c], cal[c]).abs() for c in cols], axis=1).max(axis=1)
    S["M1"] = zmax(["length", "m1_gc", "m1_maxrun"]); S["M2"] = q.m2_completeness
    S["M4"] = zmax(["m4_genes", "m4_density", "m4_meanlen"]); S["M5"] = zmax(["m5_genes", "m5_density", "m5_meanlen"])
    S["M6"] = q.m6_core; S["M8"] = q.m8_order; S["M9"] = (q.m9_gap + q.m9_codon) / 2; S["M10"] = q.m10_rbs
    S["M11"] = q.m11_pfam_share; S["M12"] = q.m12_phrog_share
    return S
calS = scores(cal)
rng = {m: (calS[m].quantile(.05), calS[m].quantile(.95)) for m in calS.columns if m != "id"}
W = "kmerarm/work/data/metrics_"
q = pd.read_csv(W + "quick.csv")
for f in ("phanotate", "pfam", "phrog"): q = q.merge(pd.read_csv(W + f + ".csv"), on="id")
c = pd.read_csv("kmerarm/checkv_out/quality_summary.tsv", sep="\t")
q = q.merge(pd.DataFrame({"id": c.contig_id, "m2_completeness": pd.to_numeric(c.completeness, errors="coerce").fillna(0)}), on="id")
q = q[~q.id.str.startswith("CAL|")]
S = scores(q); S["acc"] = S.id.str.split("__").str[0]; S["arm"] = S.id.str.split("__").str[1]
S.to_csv("results/kmer_arm_scores.csv", index=False)
raw = q.assign(acc=q.id.str.split("__").str[0], arm=q.id.str.split("__").str[1])
raw.to_csv("results/kmer_arm_raw_metrics.csv", index=False)
P = {a: S[S.arm == a].set_index("acc").sort_index() for a in ("natural", "BASE", "KMER")}
Ms = [m for m in S.columns if m.startswith("M")]
rows = []
for m in Ms:
    lo, hi = rng[m]; b, k, n = P["BASE"][m], P["KMER"][m], P["natural"][m]
    d = (k - b)
    try: p = wilcoxon(k, b, zero_method="wilcox").pvalue if (d != 0).any() else 1.0
    except ValueError: p = np.nan
    rows.append(dict(metric=m, natural_median=n.median(), BASE_median=b.median(), KMER_median=k.median(), median_diff=d.median(),
                     KMER_better=int((d > 0).sum()), KMER_worse=int((d < 0).sum()), tie=int((d == 0).sum()), wilcoxon_p=p,
                     BASE_inside=((b >= lo) & (b <= hi)).mean(), KMER_inside=((k >= lo) & (k <= hi)).mean(), natural_inside=((n >= lo) & (n <= hi)).mean()))
T = pd.DataFrame(rows).round(3); T.to_csv("results/kmer_arm_vs_base.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(T.to_string(index=False))
r = raw.pivot(index="acc", columns="arm")
print("\nraw core genes (M6, of 5): natural / BASE / KMER medians", [r["m6_core"][a].median() for a in ("natural", "BASE", "KMER")])
print("length: BASE / KMER median", r["length"]["BASE"].median(), r["length"]["KMER"].median(), " KMER shorter than natural in", int((r["length"]["KMER"] < r["length"]["natural"]).sum()), "of", len(r))
print("KMER length / natural length: median", (r["length"]["KMER"] / r["length"]["natural"]).median().round(3), " min", (r["length"]["KMER"] / r["length"]["natural"]).min().round(3))
