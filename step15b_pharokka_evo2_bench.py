"""step15b_pharokka_evo2_bench.py - the proper M12 (real Pharokka) applied to the Evo 2 outputs and the friend's benchmark sequences.
1. Evo 2: share inside the natural (calibration 5th-95th percentile) range, medians, and the prompt-excluded recount (genes that reach / lie wholly in generated letters),
   G1 minus G2 paired by genome with 2,000-bootstrap 95% intervals (seed 1337).
2. Friend's phiX174 and E. coli sequences: per-sequence proper M12 next to the simplified value.
Input: data/metrics_pharokka.csv, pharokka/out_all/pharokka_cds_final_merged_output.tsv, evo2/targets.tsv."""
import numpy as np, pandas as pd
rng = np.random.default_rng(1337)
m = pd.read_csv("data/metrics_pharokka.csv"); old = pd.read_csv("data/metrics_phrog.csv")
cal = m[m.id.str.startswith("CAL|")].m12p_phrog_share
lo, hi = cal.quantile(.05), cal.quantile(.95)
nat = m[m.id.str.endswith("__natural")].m12p_phrog_share
print(f"proper M12 calibration 5th-95th = [{lo:.2f}, {hi:.2f}]; natural test genomes inside {((nat >= lo) & (nat <= hi)).mean():.0%} (one-sided >= 5th: {(nat >= lo).mean():.0%}); natural median {nat.median():.2f}")
oldcal = old[old.id.str.startswith("CAL|")].m12_phrog_share; olo, ohi = oldcal.quantile(.05), oldcal.quantile(.95)
print(f"simplified M12 range [{olo:.2f}, {ohi:.2f}]")
ev = m[m.id.str.endswith(("__G1", "__G2"))].copy(); ev["set"] = ev.id.str[-2:]
oev = old[old.id.isin(ev.id)].set_index("id").m12_phrog_share
print("\nEvo 2: share of outputs inside the natural range, and median share of genes with a PHROG")
for g in ("G1", "G2"):
    x = ev[ev.set == g]
    print(f"  {g}: proper inside {((x.m12p_phrog_share >= lo) & (x.m12p_phrog_share <= hi)).mean():.0%}, median {x.m12p_phrog_share.median():.2f} (= {x.m12p_phrog_share.median() / nat.median():.2f}x natural)   | simplified inside "
          f"{((oev[x.id] >= olo) & (oev[x.id] <= ohi)).mean():.0%}, median {oev[x.id].median():.2f}")
# prompt-excluded recount
d = pd.read_csv("pharokka/out_all/pharokka_cds_final_merged_output.tsv", sep="\t", low_memory=False)
d = d[d.contig.isin(ev.id)].copy()
tg = pd.read_csv("evo2/targets.tsv", sep="\t").set_index("id")
d["P"] = d.contig.map(tg.prompt_len); d["lo_"] = d[["start", "stop"]].min(axis=1); d["hi_"] = d[["start", "stop"]].max(axis=1)
d["has"] = d.phrog.astype(str) != "No_PHROG"
def share(sel): return d[sel].groupby("contig").has.mean()
ver = {"all genes": share(d.hi_ > 0), "reaching generated": share(d.hi_ > d.P), "wholly generated": share(d.lo_ > d.P)}
par = sorted({i.split("__")[0] for i in ev.id})
print("\nEvo 2, proper M12 share of genes with a PHROG, median G1 / G2 and mean paired difference G1-G2 [95% CI]:")
for k, s in ver.items():
    a = s.reindex([p + "__G1" for p in par]).values.astype(float); b = s.reindex([p + "__G2" for p in par]).values.astype(float)
    ok = ~np.isnan(a) & ~np.isnan(b); dd = (a - b)[ok]
    bs = [dd[rng.integers(0, len(dd), len(dd))].mean() for _ in range(2000)]
    print(f"  {k:20s} {np.nanmedian(a):.2f} / {np.nanmedian(b):.2f}   {dd.mean():+.2f} [{np.percentile(bs, 2.5):+.2f}, {np.percentile(bs, 97.5):+.2f}] (n={ok.sum()} pairs)")
# friend's sequences
f = m[m.id.str.startswith(("phiX174_", "ecoli_"))].copy()
print("\nFriend's benchmark sequences (proper M12; simplified in brackets from the earlier run):")
simp = {}
for k in ("full", "long"):
    t = pd.read_csv(f"bench/{k}/metrics_phrog.csv") if k == "full" else None
    if t is not None: simp.update(dict(zip(t.id, t.m12_phrog_share)))
try:
    s2 = pd.read_csv("results/bench_long_scores.csv"); simp.update(dict(zip(s2.id, s2.M12)))
except Exception: pass
for r in f.itertuples():
    print(f"  {r.id:18s} genes {r.m12p_genes:3d}  PHROG share {r.m12p_phrog_share:.2f}  known-function share {r.m12p_known_share:.2f}   [simplified {simp.get(r.id, float('nan')):.2f}]")
