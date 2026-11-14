"""step9b_evo2_restricted.py - do the Evo 2 headline numbers shift when restricted to outputs with a findable F2? (user request, item 3)
Groups: all 80; intact single F2 (38); intact F2 reaching generated letters (30); no intact F2 (42).
Share inside the calibration 5th-95th percentile per metric (original scores from step4), medians, and closest damage type (step6 method)."""
import pandas as pd
S = pd.read_csv("results/scores_all.csv")
f = pd.read_csv("evo2/f2_findable.tsv", sep="\t")
man = pd.read_csv("data/ladder_manifest.tsv", sep="\t", dtype={"severity": str}, quoting=3)
M = ["M1", "M2", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12"]
cal = S[S.id.str.startswith("CAL|")]
ev = S[S.id.isin(f.id)].merge(f, on="id"); ev["set"] = ev.id.str[-2:]
groups = {"all 80": ev, "intact F2 (38)": ev[ev.intact_F2 == 1], "intact F2 reaching generated (30)": ev[ev.intact_F2_reaching_generated == 1],
          "no intact F2 (42)": ev[ev.intact_F2 == 0]}
rows = []
for m in M:
    lo, hi = cal[m].quantile(.05), cal[m].quantile(.95)
    r = {"metric": m}
    for g, d in groups.items():
        r[g + " inside%"] = round(100 * ((d[m] >= lo) & (d[m] <= hi)).mean())
    rows.append(r)
T = pd.DataFrame(rows); pd.set_option("display.width", 220); print("share (%) of outputs inside the natural (calibration 5-95%) range:"); print(T.to_string(index=False))
print("\nmedian M6 (core genes of 5) and M2 (CheckV completeness %):")
for g, d in groups.items(): print(f"  {g:36s} n={len(d):2d}  M6 {d.M6.median():.0f}   M2 {d.M2.median():.0f}   M11 {d.M11.median():.2f}   M12 {d.M12.median():.2f}")
# closest damage type, same method as step6
lad = man.merge(S, on="id"); lad["level"] = lad.damage + "_" + lad.severity.astype(str)
keep = lad[lad.damage.isin(["natural", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9"])]
sd = S[~S.id.str.contains(r"__G[12]$") & ~S.id.str.startswith("CAL|")][M].std().replace(0, 1)
prof = keep.groupby("level")[M].mean() / sd
def closest(row): return ((prof - row[M].astype(float) / sd) ** 2).sum(axis=1).pow(.5).idxmin().split("_")[0]
ev["type"] = ev.apply(closest, axis=1)
print("\nclosest damage type (count):")
for g, d in groups.items(): print(f"  {g:36s}", d.apply(closest, axis=1).value_counts().to_dict())
