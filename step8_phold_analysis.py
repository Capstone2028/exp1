"""step8_phold_analysis.py - M13 (Phold structure labels): share of called genes labelled, per sequence, on the 220-sequence
subset (20 test genomes x 11 variants). Same AUROC rules as step4 (natural vs damaged, 2,000 bootstrap resamples of genomes, seed 1337),
but only 20 genomes -> intervals are wide, so this is indicative. Higher share = more natural.
Output: data/metrics_phold.csv, results/phold_auroc.csv"""
import numpy as np, pandas as pd
d = pd.read_csv("phold/out/phold_per_cds_predictions.tsv", sep="\t")
d["labelled"] = d.phrog != "No_PHROG"
s = d.groupby("contig_id").agg(m13_genes=("cds_id", "size"), m13_labelled=("labelled", "sum")).reset_index()
s["m13_share"] = (s.m13_labelled / s.m13_genes).round(4)
s = s.rename(columns={"contig_id": "id"}); s.to_csv("data/metrics_phold.csv", index=False)
s["parent"] = s.id.str.split("__").str[0]; s["var"] = s.id.str.split("__").str[1]
genomes = sorted(s.parent.unique())
piv = s.pivot(index="parent", columns="var", values="m13_share").loc[genomes]
rng = np.random.default_rng(1337)
auc = lambda a, b: float((a[:, None] > b[None, :]).mean() + 0.5 * (a[:, None] == b[None, :]).mean())
rows = []
for v in piv.columns:
    if v == "natural": continue
    a, b = piv["natural"].values, piv[v].values
    bs = [auc(a[i], b[i]) for i in (rng.integers(0, len(a), len(a)) for _ in range(2000))]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    rows.append(dict(level=v, auroc=round(auc(a, b), 2), lo=round(lo, 2), hi=round(hi, 2)))
R = pd.DataFrame(rows); R.to_csv("results/phold_auroc.csv", index=False)
print(R.to_string(index=False))
print("natural M13 median", piv["natural"].median(), "| M11 on same genomes:")
q = pd.read_csv("data/metrics_pfam.csv"); q = q[q.id.isin(s.id)]
print("Pfam share natural median", q[q.id.str.endswith("natural")].m11_pfam_share.median())
