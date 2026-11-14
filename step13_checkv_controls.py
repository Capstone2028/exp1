"""step13_checkv_controls.py - is the CheckV completeness gap between G1 (82%) and G2 (17%) just the verbatim prompt?
Controls (all CheckV 1.1.1, db v1.5, per Evo 2 output, 40 G1 + 40 G2):
  orig : the full Evo 2 output (as reported)
  PO   : the prompt alone (1,000 or 300 letters) = also the length-matched natural truncation of the first P letters
  GO   : the generated letters alone (prompt removed)
  NR   : the natural parent's matching region (same letters, same length as GO): length-matched natural
  NW   : a random window of the natural parent, same length as GO
CheckV 'Not-determined' completeness is counted as 0 and reported separately. Differences are paired by genome; 95% intervals from
2,000 bootstrap resamples of genomes (seed 1337)."""
import numpy as np, pandas as pd
rng = np.random.default_rng(1337)
c = pd.read_csv("evo2/checkv_ctrl/out/quality_summary.tsv", sep="\t")
o = pd.read_csv("evo2/checkv_out/quality_summary.tsv", sep="\t")
o["contig_id"] = "OR@" + o.contig_id
c = pd.concat([c, o])
c["kind"] = c.contig_id.str.split("@").str[0]; c["id"] = c.contig_id.str.split("@").str[1]
c["comp"] = pd.to_numeric(c.completeness, errors="coerce"); c["nd"] = c.comp.isna()
c["comp"] = c.comp.fillna(0); c["set"] = c.id.str[-2:]; c["parent"] = c.id.str.split("__").str[0]
W = c.pivot(index="id", columns="kind", values="comp")
ND = c.pivot(index="id", columns="kind", values="nd")
names = {"OR": "original (full output)", "PO": "prompt only", "GO": "generated only", "NR": "natural, same region+length as GO", "NW": "natural, random window, same length"}
print("CheckV completeness %, median [IQR] (and number of 'not-determined' out of 40):")
for k in ("OR", "PO", "GO", "NR", "NW"):
    row = []
    for g in ("G1", "G2"):
        v = W[k][W.index.str.endswith(g)]
        row.append(f"{v.median():5.1f} [{v.quantile(.25):.0f}-{v.quantile(.75):.0f}] nd={int(ND[k][ND.index.str.endswith(g)].sum())}")
    print(f"  {names[k]:38s} G1 {row[0]:26s} G2 {row[1]}")
def paired(a, b):                                   # G1-minus-G2 on (a - b) per output, paired by genome
    ids1 = [i for i in W.index if i.endswith("G1")]; par = [i.split("__")[0] for i in ids1]
    x = np.array([(W.loc[p + "__G1", a] - (W.loc[p + "__G1", b] if b else 0)) - (W.loc[p + "__G2", a] - (W.loc[p + "__G2", b] if b else 0)) for p in par])
    bs = [x[rng.integers(0, len(x), len(x))].mean() for _ in range(2000)]
    return x.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5)
print("\nG1 minus G2 gap in completeness (percentage points, mean over 40 paired genomes) [95% CI]:")
for lab, a, b in (("original", "OR", None), ("prompt floor (prompt only)", "PO", None), ("original minus prompt-only", "OR", "PO"),
                  ("generated only", "GO", None), ("natural, matched to generated-only", "NR", None)):
    m, lo, hi = paired(a, b); print(f"  {lab:38s} {m:+6.1f} [{lo:+.1f}, {hi:+.1f}]")
print("\nper-genome: generated-only completeness vs its length-matched natural (NR), G1 / G2 medians:",
      [(g, round(W['GO'][W.index.str.endswith(g)].median(), 1), round(W['NR'][W.index.str.endswith(g)].median(), 1)) for g in ("G1", "G2")])
W.to_csv("results/evo2_checkv_controls.csv")
