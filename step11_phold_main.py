"""step11_phold_main.py - items 1 and 2 of the user's list: Phold on the 60 natural test genomes, 80 Evo 2 outputs (G1, G2) and 30 calibration genomes.
Phold = ProstT5 + Foldseek (structure-based labels); M12 here = the simplified MMseqs2 search against PHROG family sequences, run on the SAME
proteins (Phold's own pyrodigal-gv gene calls), so the comparison is protein by protein. Higher share labelled = more natural-looking.
Natural pass rate (design section 7): share of the 60 natural TEST genomes whose share-labelled is inside the calibration (30 genomes) 5th-95th percentile."""
import numpy as np, pandas as pd
d = pd.read_csv("phold/out2/phold_per_cds_predictions.tsv", sep="\t")
h = pd.read_csv("phold/out2/m12_hits.tsv", sep="\t", header=None, names=["q", "t", "e", "b"])
m12 = set(h.q.str.split(":", n=1).str[1])
d["phold"] = d.phrog != "No_PHROG"
d["phold_high"] = d.phold & d.annotation_confidence.eq("high")
d["m12"] = d.cds_id.isin(m12)
d["known"] = d.phold & (d.function != "unknown function")
tg = pd.read_csv("evo2/targets.tsv", sep="\t").set_index("id")
def grp(c):
    if c.startswith("CAL_"): return "calibration"
    if c.endswith("__natural"): return "natural test"
    return "Evo2 " + c[-2:]
d["group"] = d.contig_id.map(grp)
d["inprompt"] = [(r.group.startswith("Evo2") and r.end <= tg.loc[r.contig_id, "prompt_len"]) for r in d.itertuples()]
per = d.groupby(["contig_id", "group"]).agg(genes=("cds_id", "size"), phold=("phold", "mean"), high=("phold_high", "mean"), m12=("m12", "mean"), known=("known", "mean")).reset_index()
print("sequences:", per.group.value_counts().to_dict(), "| genes:", d.group.value_counts().to_dict())
print("\nshare of genes labelled, per sequence (median [IQR]):")
for g in ["calibration", "natural test", "Evo2 G1", "Evo2 G2"]:
    x = per[per.group == g]
    f = lambda c: f"{x[c].median():.2f} [{x[c].quantile(.25):.2f}-{x[c].quantile(.75):.2f}]"
    print(f"  {g:13s} n={len(x):2d}  Phold {f('phold')}   Phold-high {f('high')}   M12-simplified {f('m12')}   genes/seq {x.genes.median():.0f}")
cal = per[per.group == "calibration"]
for col, nm in (("phold", "M13 Phold"), ("m12", "M12 simplified (same proteins)")):
    lo, hi = cal[col].quantile(.05), cal[col].quantile(.95)
    nat = per[per.group == "natural test"][col]
    print(f"\n{nm}: calibration 5th-95th = [{lo:.2f}, {hi:.2f}]; natural test pass (two-sided) {((nat >= lo) & (nat <= hi)).mean():.0%}, one-sided (>= 5th) {(nat >= lo).mean():.0%}")
    for g in ("Evo2 G1", "Evo2 G2"):
        e = per[per.group == g][col]
        print(f"   {g}: inside range {((e >= lo) & (e <= hi)).mean():.0%};  median {e.median():.2f} = {e.median() / nat.median():.2f}x the natural median {nat.median():.2f}")
print("\nprotein by protein, Phold vs simplified M12 (rows = share of that group's genes):")
for g in ["natural test", "Evo2 G1", "Evo2 G2"]:
    x = d[d.group == g]; n = len(x)
    print(f"  {g:13s} n={n:4d}  both {(x.phold & x.m12).sum() / n:.2f}   Phold only {(x.phold & ~x.m12).sum() / n:.2f}   M12 only {(~x.phold & x.m12).sum() / n:.2f}   neither {(~x.phold & ~x.m12).sum() / n:.2f}")
x = d[(d.group == "natural test") & d.phold & ~d.m12]
print("\nwhat Phold labels in natural genomes that M12 misses (n=%d): by function" % len(x)); print(x.function.value_counts().head(8).to_string())
print("   confidence:", x.annotation_confidence.value_counts().to_dict())
y = d[d.group.str.startswith("Evo2") & d.phold & ~d.m12]
print("Evo 2 genes Phold labels and M12 misses (n=%d): by function" % len(y)); print(y.function.value_counts().head(6).to_string())
print("\nEvo 2 genes wholly inside the copied prompt:", int(d.inprompt.sum()), "of", int(d.group.str.startswith("Evo2").sum()))
for g in ("Evo2 G1", "Evo2 G2"):
    x = d[d.group == g]
    print(f"  {g}: share labelled, genes inside prompt {x[x.inprompt].phold.mean():.2f} (n={int(x.inprompt.sum())}); genes reaching generated letters {x[~x.inprompt].phold.mean():.2f} (n={int((~x.inprompt).sum())}); high-confidence {x[~x.inprompt].phold_high.mean():.2f}")
nn = d[d.group == "natural test"]
print(f"  natural test genes: {nn.phold.mean():.2f} labelled, {nn.phold_high.mean():.2f} high confidence")
per.to_csv("results/phold_main_per_sequence.csv", index=False)
