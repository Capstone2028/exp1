"""step8b_phold_rotated.py - M13 (Phold) on the 220-sequence subset (20 genomes x 11 variants) with the start-position fix applied.
The fix (step7_standardise.py rule) moved only some sequences; those were re-run through Phold (phold/out3) and replace their rows from the
first run (phold/out). Everything else is identical to step8_phold_analysis.py: share of genes labelled, AUROC natural vs damaged,
2,000 bootstrap resamples of the 20 genomes, seed 1337. Output: results/phold_auroc_rotated.csv"""
import numpy as np, pandas as pd
a = pd.read_csv("phold/out/phold_per_cds_predictions.tsv", sep="\t")
b = pd.read_csv("phold/out3/phold_per_cds_predictions.tsv", sep="\t")
moved = set(b.contig_id) | set(l[1:].split()[0] for l in open("phold/phold_in3.fna") if l.startswith(">"))
a["labelled"] = a.phrog != "No_PHROG"; b["labelled"] = b.phrog != "No_PHROG"
d = pd.concat([a[~a.contig_id.isin(moved)], b])
def table(df):
    s = df.groupby("contig_id").agg(g=("cds_id", "size"), l=("labelled", "sum")).reset_index(); s["share"] = s.l / s.g
    s["parent"] = s.contig_id.str.split("__").str[0]; s["var"] = s.contig_id.str.split("__").str[1]
    return s.pivot(index="parent", columns="var", values="share")
auc = lambda x, y: float((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean())
def run(piv, seed=1337):
    rng = np.random.default_rng(seed); out = {}
    for v in piv.columns:
        if v == "natural": continue
        x, y = piv["natural"].values, piv[v].values
        bs = [auc(x[i], y[i]) for i in (rng.integers(0, len(x), len(x)) for _ in range(2000))]
        out[v] = (round(auc(x, y), 2), round(np.percentile(bs, 2.5), 2), round(np.percentile(bs, 97.5), 2))
    return out
o = run(table(a)); n = run(table(d))
R = pd.DataFrame({"orig": {k: f"{v[0]:.2f} ({v[1]:.2f}-{v[2]:.2f})" for k, v in o.items()}, "start-fixed": {k: f"{v[0]:.2f} ({v[1]:.2f}-{v[2]:.2f})" for k, v in n.items()}})
R.to_csv("results/phold_auroc_rotated.csv"); print(R.to_string())
print("sequences re-run with the start fix:", len(moved), "| natural genes labelled median:", round(table(d)["natural"].median(), 3))
