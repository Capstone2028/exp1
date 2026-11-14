"""step12b_compare.py - original vs prompt-excluded Evo 2 protein-level numbers, G1 (1,000-letter prompt) vs G2 (300), side by side.
'orig' = every gene; 'ov' = genes overlapping generated letters (end after the prompt); 'wh' = wholly generated genes (start after the prompt).
G1-G2 difference: paired by genome (the same 40 genomes feed both), 95% interval from 2,000 bootstrap resamples of genomes, seed 1337.
A positive interval that excludes 0 means G1 is higher (more natural-looking) than G2."""
import numpy as np, pandas as pd
rng = np.random.default_rng(1337)
pe = pd.read_csv("evo2/prompt_excluded.csv")
pe["parent"] = pe.id.str.split("__").str[0]
tg = pd.read_csv("evo2/targets.tsv", sep="\t").set_index("id")
ph = pd.read_csv("phold/out2/phold_per_cds_predictions.tsv", sep="\t")
ph = ph[ph.contig_id.str.contains(r"__G[12]$")].copy()
ph["P"] = ph.contig_id.map(tg.prompt_len)
h = pd.read_csv("phold/out2/m12_hits.tsv", sep="\t", header=None, names=["q", "t", "e", "b"])
m12 = set(h.q.str.split(":", n=1).str[1])
ph["phold"] = ph.phrog != "No_PHROG"; ph["high"] = ph.phold & ph.annotation_confidence.eq("high"); ph["m12"] = ph.cds_id.isin(m12)
def share(df, col, tag):
    keep = {"orig": df, "ov": df[df.end > df.P], "wh": df[df.start > df.P]}[tag]
    return keep.groupby("contig_id")[col].mean()
tab = {}
for tag in ("orig", "ov", "wh"):
    t = pe.set_index("id")
    tab[("M6 core genes (of 5)", tag)] = t[f"m6_{tag}"]
    if tag != "orig":
        tab[("M11 Pfam share", tag)] = pd.to_numeric(t[f"pfam_{tag}"], errors="coerce")
    tab[("M12 simplified share (Phold proteins)", tag)] = share(ph, "m12", tag)
    tab[("M13 Phold share", tag)] = share(ph, "phold", tag)
    tab[("M13 Phold high-confidence share", tag)] = share(ph, "high", tag)
orig11 = pd.read_csv("data/metrics_pfam.csv").set_index("id").m11_pfam_share
tab[("M11 Pfam share", "orig")] = orig11[orig11.index.isin(pe.id)]
n_genes = {tag: pe.set_index("id")[f"genes_{tag}"] for tag in ("orig", "ov", "wh")}
print("genes per output (median): orig %s  ov %s  wh %s;  outputs with zero qualifying genes: ov %d, wh %d" % (
    n_genes["orig"].median(), n_genes["ov"].median(), n_genes["wh"].median(), (n_genes["ov"] == 0).sum(), (n_genes["wh"] == 0).sum()))
for g in ("G1", "G2"):
    x = pe[pe.set == g]; print(f"  {g}: genes orig {x.genes_orig.sum()}, overlapping generated {x.genes_ov.sum()}, wholly generated {x.genes_wh.sum()}; outputs with no wholly-generated gene: {(x.genes_wh == 0).sum()}")
def med(s, ids): v = s.reindex(ids).dropna(); return v.median() if len(v) else np.nan
def diff(s, tag):
    a = s.reindex([f"{p}__G1" for p in sorted(pe.parent.unique())]).values
    b = s.reindex([f"{p}__G2" for p in sorted(pe.parent.unique())]).values
    ok = ~np.isnan(a.astype(float)) & ~np.isnan(b.astype(float)); a, b = a[ok].astype(float), b[ok].astype(float)
    d = a - b
    bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)]
    return d.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(d)
g1 = pe[pe.set == "G1"].id; g2 = pe[pe.set == "G2"].id
print("\nmedian per output | G1 | G2 | G1-G2 mean difference [95% CI] (n pairs)")
for m in dict.fromkeys(k[0] for k in tab):
    for tag in ("orig", "ov", "wh"):
        if (m, tag) not in tab: continue
        s = tab[(m, tag)]; d, lo, hi, n = diff(s, tag)
        flag = "G1 higher" if lo > 0 else ("G2 higher" if hi < 0 else "no clear difference")
        print(f"{m:40s} {tag:5s} {med(s, g1):6.2f} {med(s, g2):6.2f}   {d:+.2f} [{lo:+.2f}, {hi:+.2f}] n={n}  {flag}")
    print()
pe.to_csv("results/evo2_prompt_excluded_per_sequence.csv", index=False)
