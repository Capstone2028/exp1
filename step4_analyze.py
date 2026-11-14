"""step4_analyze.py - apply the design's section-7 decision rules to the metrics computed so far.
Run on Windows: python step4_analyze.py      Output: results/*.csv, results/grid.png, results/verdicts.md

FIXED BEFORE ANY AUROC WAS COMPUTED (medians by damage type had been looked at in earlier steps; no tests had):
  Score per metric (direction decided in advance; higher score = more natural):
    M1  -max(|z_length|, |z_GC|, |z_longest_run|)     z = (x - calibration median) / calibration IQR (or std if IQR = 0)
    M3  -end_overlap
    M4  -max(|z_genes|, |z_density|, |z_mean_len|)    same z rule, pyrodigal numbers
    M5  same with PHANOTATE numbers
    M6  core families found (0-5)           M7 valid-ORF share     M8 gene-order share
    M9  mean(gap share, codon share)        M10 RBS share          M11 Pfam share   M12 PHROG any-hit share
  M2 CheckV completeness % (not-determined counted as 0; added 30 Oct after CheckV finished, same rules).
  Not run (tools not installed): M13 Phold, M14 ESMFold, M15 Foldseek.
  AUROC = P(natural genome scores higher than its damaged version) + 0.5 P(tie), over the 60 test genomes;
          95% interval = 2,000 bootstrap resamples of genomes (natural and damaged resampled together), seed 1337.
  Middle severity: D1 30%, D2 2%, D3 6 (index n//2), D4 3, D5 3, D6 4 cuts, D7 50%, D8, D9; D10 revcomp and rotation separately.
  Sensitivity pass : for EVERY damage predicted 'catch' (H) in the design table: AUROC >= 0.85 and lower bound >= 0.75.
  False alarms pass: AUROC for D1(30%), D10 revcomp and D10 rotation all within 0.40-0.60.
  Independence pass: if predicted H at D7 or D8: at least one of them passes the sensitivity test; (n/a if no H there).
                     A metric with AUROC < 0.6 at D8 is labelled 'texture-level only' when it is meant to measure genes/proteins.
  Natural pass rate: share of the 60 natural test genomes whose score is inside the calibration 5th-95th percentile (two-sided, as written);
                     the one-sided share (>= 5th percentile) is also reported. Need >= 90%.
  Agreement: Spearman M4 vs M5 and among M11/M12, over all ladder sequences and over natural genomes only.
  Verdict: Keep = all applicable checks pass; Drop for validity = no predicted-H damage passes (or no H at all and no reaction anywhere),
           or both false-alarm and sensitivity fail; Keep with caveat = anything else (the failing checks are listed).
"""
import json, os
import numpy as np, pandas as pd

rng = np.random.default_rng(1337)
os.makedirs("results", exist_ok=True)
q = pd.read_csv("data/metrics_quick.csv").merge(pd.read_csv("data/metrics_phanotate.csv"), on="id") \
    .merge(pd.read_csv("data/metrics_pfam.csv"), on="id").merge(pd.read_csv("data/metrics_checkv.csv"), on="id").merge(pd.read_csv("data/metrics_phrog.csv"), on="id")
cal = q[q.id.str.startswith("CAL|")].copy()


def z(col, cal_col):
    med = cal_col.median(); iqr = cal_col.quantile(.75) - cal_col.quantile(.25)
    sd = iqr if iqr > 0 else (cal_col.std() or 1.0)
    return (col - med) / sd


def zmax(cols):
    return -pd.concat([z(q[c], cal[c]).abs() for c in cols], axis=1).max(axis=1)


S = pd.DataFrame({"id": q.id})
S["M1"] = zmax(["length", "m1_gc", "m1_maxrun"])
S["M2"] = q.m2_completeness
S["M3"] = -q.m3_end_overlap
S["M4"] = zmax(["m4_genes", "m4_density", "m4_meanlen"])
S["M5"] = zmax(["m5_genes", "m5_density", "m5_meanlen"])
S["M6"] = q.m6_core; S["M7"] = q.m7_valid; S["M8"] = q.m8_order
S["M9"] = (q.m9_gap + q.m9_codon) / 2; S["M10"] = q.m10_rbs
S["M11"] = q.m11_pfam_share; S["M12"] = q.m12_phrog_share
S.to_csv("results/scores_all.csv", index=False)
METRICS = [f"M{i}" for i in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
NAMES = {"M1": "Basic gates", "M2": "CheckV", "M3": "Closed-genome", "M4": "Pyrodigal", "M5": "PHANOTATE", "M6": "Core genes", "M7": "ORF integrity",
         "M8": "Gene order", "M9": "Start/stop", "M10": "RBS signal", "M11": "Pfam", "M12": "PHROG"}
PRED = json.load(open("data/predictions.json"))            # damage order D2..D9
DAM = ["D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9"]

man = pd.read_csv("data/ladder_manifest.tsv", sep="	", dtype={"severity": str}, quoting=3)
man = man.merge(S, on="id")
genomes = sorted(man.parent.unique()); gi = {g: i for i, g in enumerate(genomes)}
nat = man[man.damage == "natural"].set_index("parent").loc[genomes]
LEVELS = {"D1_30": ("D1", "30"), "D2_2": ("D2", "2"), "D3_6": ("D3", "6"), "D4_3": ("D4", "3"), "D5_3": ("D5", "3"),
          "D6_4": ("D6", "4"), "D7_50": ("D7", "50"), "D8": ("D8", "1"), "D9": ("D9", "1"),
          "D10_revcomp": ("D10", "revcomp"), "D10_rotation": ("D10", "rotation")}
ALL_LEVELS = {f"{d}_{s}": (d, s) for d, s in man[~man.damage.isin(["natural"])][["damage", "severity"]].drop_duplicates().values}


def auroc(a, b):                         # P(a > b) + 0.5 P(a == b), all pairs
    return float((a[:, None] > b[None, :]).mean() + 0.5 * (a[:, None] == b[None, :]).mean())


def boot(a, b, n=2000):
    k = len(a); out = np.empty(n)
    for i in range(n):
        ix = rng.integers(0, k, k); out[i] = auroc(a[ix], b[ix])
    return np.percentile(out, [2.5, 97.5])


rows, rows_all = [], []
for m in METRICS:
    ns = nat[m].values
    for lv, (d, s) in ALL_LEVELS.items():
        ds = man[(man.damage == d) & (man.severity == s)].set_index("parent").loc[genomes][m].values
        rows_all.append((m, lv, auroc(ns, ds)))
    for lv, (d, s) in LEVELS.items():
        ds = man[(man.damage == d) & (man.severity == s)].set_index("parent").loc[genomes][m].values
        lo, hi = boot(ns, ds)
        rows.append(dict(metric=m, level=lv, auroc=auroc(ns, ds), lo=lo, hi=hi))
G = pd.DataFrame(rows); G.to_csv("results/auroc_mid.csv", index=False)
pd.DataFrame(rows_all, columns=["metric", "level", "auroc"]).to_csv("results/auroc_all_severities.csv", index=False)
grid = G.pivot(index="metric", columns="level", values="auroc").loc[METRICS, list(LEVELS)]

# dose response
dose = {}
for m in METRICS:
    dose[m] = {}
    for d in ["D2", "D3", "D4", "D5", "D6", "D7"]:
        sub = man[man.damage == d]
        dose[m][d] = sub.severity.astype(float).corr(sub[m], method="spearman")
dose = pd.DataFrame(dose).T; dose.to_csv("results/dose_spearman.csv")

# natural pass rate
nat_pass = {}
for m in METRICS:
    c = S.loc[S.id.str.startswith("CAL|"), m]; lo_, hi_ = c.quantile(.05), c.quantile(.95)
    v = nat[m].values
    nat_pass[m] = (float(((v >= lo_) & (v <= hi_)).mean()), float((v >= lo_).mean()))

# verdicts
verd = []
for m in METRICS:
    pm = PRED.get(m, (None, "NNNNNNNN"))[1] if m != "M12" else PRED["M12"][1]
    H = [d for d, c in zip(DAM, pm) if c == "H"]
    g = {r.level: r for r in G[G.metric == m].itertuples()}
    sens = {d: (g[d if d in ("D8", "D9") else [k for k in LEVELS if k.startswith(d + "_")][0]]) for d in H}
    ok = {d: (r.auroc >= .85 and r.lo >= .75) for d, r in sens.items()}
    sens_pass = all(ok.values()) if H else None
    fa = {k: g[k].auroc for k in ("D1_30", "D10_revcomp", "D10_rotation")}
    fa_pass = all(.4 <= v <= .6 for v in fa.values())
    ind_ds = [d for d in H if d in ("D7", "D8")]
    ind_pass = any(ok[d] for d in ind_ds) if ind_ds else None
    texture = g["D8"].auroc < .6
    np2, np1 = nat_pass[m]
    npass = np2 >= .9
    fails = []
    if sens_pass is False: fails.append("sensitivity (" + ",".join(d for d in H if not ok[d]) + ")")
    if not fa_pass: fails.append("false alarms (" + ",".join(f"{k}={v:.2f}" for k, v in fa.items() if not .4 <= v <= .6) + ")")
    if ind_pass is False: fails.append("independence")
    if not npass: fails.append(f"natural pass rate {np2:.0%}")
    if (not H) and not any(abs(v - .5) > .1 for k, v in ((r.level, r.auroc) for r in g.values())):
        v_ = "Drop for validity (no reaction to any damage)"
    elif H and not any(ok.values()) and not fa_pass:
        v_ = "Drop for validity"
    elif H and not any(ok.values()):
        v_ = "Drop for validity (catches none of its predicted damages)"
    elif not fails:
        v_ = "Keep"
    else:
        v_ = "Keep with caveat: " + "; ".join(fails)
    verd.append(dict(metric=m, name=NAMES[m], predicted_catch=",".join(H) or "-", sens_pass=sens_pass, falsealarm_pass=fa_pass,
                     independence=ind_pass, nat_pass_two_sided=round(np2, 3), nat_pass_one_sided=round(np1, 3),
                     d8_texture_only=bool(texture), verdict=v_))
V = pd.DataFrame(verd); V.to_csv("results/verdicts.csv", index=False)

# agreement
allS = S[~S.id.str.startswith("CAL|") & ~S.id.str.contains(r"__G[12]$")]; natS = allS[allS.id.str.endswith("__natural")]   # ladder sequences only (Evo 2 outputs excluded)
agree = {"M4~M5 all": allS.M4.corr(allS.M5, method="spearman"), "M4~M5 natural": natS.M4.corr(natS.M5, method="spearman"),
         "M11~M12 all": allS.M11.corr(allS.M12, method="spearman"), "M11~M12 natural": natS.M11.corr(natS.M12, method="spearman")}

# outputs
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print("AUROC at middle severity (natural vs damaged; 0.5 = no difference; below 0.5 means damaged scored HIGHER):")
print(grid.round(2).to_string())
print("\nDose response (Spearman, severity vs score; negative = score falls with damage):")
print(dose.round(2).to_string())
print("\nNatural pass rate (two-sided, one-sided):", {k: tuple(round(x, 2) for x in v) for k, v in nat_pass.items()})
print("\nAgreement:", {k: round(v, 2) for k, v in agree.items()})
print("\nVerdicts:"); print(V[["metric", "name", "predicted_catch", "verdict"]].to_string(index=False))

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 5.2))
    im = ax.imshow(grid.values, cmap="RdBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(grid.shape[1]), grid.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(METRICS)), [f"{m} {NAMES[m]}" for m in METRICS])
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            ax.text(j, i, f"{grid.values[i, j]:.2f}", ha="center", va="center", fontsize=7, color="black")
    fig.colorbar(im, label="AUROC natural vs damaged (0.5 = cannot tell)"); fig.tight_layout(); fig.savefig("results/grid.png", dpi=150)
except Exception as e:
    print("plot skipped:", e)
with open("results/agreement.json", "w") as fh: json.dump(agree, fh, indent=1)
