"""step6_evo2_analysis.py - score the 80 Evo 2 outputs (design section 6).
(1) Distribution of each metric vs the natural range: share of outputs inside the calibration 5th-95th percentile
    (scores from step4: higher = more natural), plus medians for natural genomes, G1 (1,000-letter prompt), G2 (300-letter prompt).
(2) 'Equivalent damage': for each output, the damage level whose mean metric profile is closest (Euclidean distance on
    standardised scores of M1, M2, M4-M12; M3 excluded because it never varies). Candidates: natural and every D2-D9 level
    (D1/D10 are harmless and look natural). A crude summary, not a claim that Evo 2 'applies' that damage."""
import numpy as np, pandas as pd
S = pd.read_csv("results/scores_all.csv")
man = pd.read_csv("data/ladder_manifest.tsv", sep="\t", dtype={"severity": str}, quoting=3)
M = ["M1", "M2", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12"]
cal = S[S.id.str.startswith("CAL|")]
ev = S[S.id.str.contains(r"__G[12]$")].copy(); ev["set"] = ev.id.str[-2:]
nat = S[S.id.str.endswith("__natural")]
rows = []
for m in M:
    lo, hi = cal[m].quantile(.05), cal[m].quantile(.95)
    r = dict(metric=m, natural_median=nat[m].median(), G1_median=ev[ev.set == "G1"][m].median(), G2_median=ev[ev.set == "G2"][m].median(),
             natural_inside=((nat[m] >= lo) & (nat[m] <= hi)).mean(), G1_inside=((ev[ev.set == "G1"][m] >= lo) & (ev[ev.set == "G1"][m] <= hi)).mean(),
             G2_inside=((ev[ev.set == "G2"][m] >= lo) & (ev[ev.set == "G2"][m] <= hi)).mean())
    rows.append(r)
T = pd.DataFrame(rows).round(2); T.to_csv("results/evo2_vs_natural.csv", index=False)
pd.set_option("display.width", 200); print(T.to_string(index=False))
# equivalent damage
lad = man.merge(S, on="id"); lad["level"] = lad.damage + "_" + lad.severity.astype(str)
keep = lad[lad.damage.isin(["natural", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9"])]
sd = S[~S.id.str.contains(r"__G[12]$") & ~S.id.str.startswith("CAL|")][M].std().replace(0, 1)
prof = keep.groupby("level")[M].mean()
Z = lambda df: df[M] / sd
out = []
for r in ev.itertuples():
    v = (pd.Series(r._asdict())[M].astype(float)) / sd
    d = ((prof / sd) - v).pow(2).sum(axis=1) ** .5
    out.append(dict(id=r.id, set=r.set, closest=d.idxmin(), dist=round(d.min(), 2), natural_dist=round(d["natural_nan" if False else "natural_-"], 2)))
O = pd.DataFrame(out); O["type"] = O.closest.str.split("_").str[0]; O.to_csv("results/evo2_equivalent_damage.csv", index=False)
print("\nclosest damage type (count):"); print(pd.crosstab(O.type, O.set).to_string())
print("\nmost common closest levels:"); print(O.groupby("set").closest.agg(lambda s: s.value_counts().head(5).to_dict()).to_string())
print("\nmedian distance to closest vs to natural:", O.dist.median(), O.natural_dist.median())
