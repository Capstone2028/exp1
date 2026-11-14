"""step3_quick_metrics.py - the metrics that need no big databases: M1, M3, M4, M6, M7, M8, M9, M10.
Run in WSL: ~/kmer-venv/bin/python step3_quick_metrics.py        (needs pyrodigal and blastp)
Input : data/ladder.fna (60 natural + 1,620 damaged test sequences), data/natural_calib.fna (30 calibration genomes)
Output: data/metrics_quick.csv  (one row per sequence)

Plain English for each metric (all scores are computed from the DNA alone, per sequence):
  M1  length, GC fraction, letters other than ACGT, longest run of one repeated letter
  M3  closed-genome check: longest stretch where the last letters equal the first letters (>= 8)
  M4  pyrodigal (meta mode) gene finder: gene count, coding density (share of letters inside genes), mean gene length
  M6  core-gene checklist: how many of the 5 shared gene families (F1-F5) have a called protein that matches a
      calibration-set protein of that family (blastp e<=1e-5, aligned part >= 50% of the shorter protein)
  M7  ORF integrity: share of called genes that are whole (not cut by the sequence end), start with ATG/GTG/TTG,
      have whole codons, end with a stop and contain no internal stop
  M8  gene order: the genome is read as a circle; neighbouring F1-F5 genes are paired; score = share of pairs
      that also occur (same order and direction) in at least half of the calibration genomes
  M9  start/stop layout (simplified, ignores where the genome starts): share of neighbouring-gene gaps inside the
      calibration 5th-95th percentile range (m9_gap) and share of genes whose start+stop codon pair is seen in >=1%
      of calibration genes (m9_codon)
  M10 ribosome-binding-site motif (GGAG or GAGG) 4-20 letters before the gene start: share of genes with one
"""
import collections, csv, re, subprocess, sys
import numpy as np, pyrodigal

COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
STARTS = ("ATG", "GTG", "TTG")
STOPS = ("TAA", "TAG", "TGA")
CODON = {}
for i, a in enumerate("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    CODON["".join("TCAG"[x] for x in (i // 16, (i // 4) % 4, i % 4))] = a


def read_fa(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        else:
            d[n].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


finder = pyrodigal.GeneFinder(meta=True)


def genes_of(seq):
    out = []
    for g in finder.find_genes(seq.encode()):
        a, e = g.begin - 1, g.end
        cs = seq[a:e] if g.strand == 1 else rc(seq[a:e])
        out.append(dict(a=a, e=e, st=g.strand, partial=bool(g.partial_begin or g.partial_end), cs=cs,
                        prot="".join(CODON.get(cs[i:i + 3], "X") for i in range(0, len(cs) - 3, 3))))
    return out


def valid_orf(g):
    cs = g["cs"]
    return (not g["partial"] and len(cs) % 3 == 0 and cs[:3] in STARTS and cs[-3:] in STOPS
            and all(CODON.get(cs[i:i + 3]) != "*" for i in range(0, len(cs) - 3, 3)))


# ---- calibration references ----------------------------------------------------------------------
calib = read_fa("data/natural_calib.fna")
ref_prot = {}
ref_pairs_count = collections.Counter()
gaps, codon_pairs, ncalg = [], collections.Counter(), 0
fam_by = collections.defaultdict(list)
for l in list(open("data/genes_oriented.tsv"))[1:]:
    g, st, i, a, b, s, f, pl = l.rstrip("\n").split("\t")
    if st == "calibration" and f:
        fam_by[g].append((int(a), int(b), int(s), f))
with open("data/ref_proteins.faa", "w") as fh:
    for g, lst in fam_by.items():
        for a, b, s, f in lst:
            cs = calib[g][a:b] if s == 1 else rc(calib[g][a:b])
            p = "".join(CODON.get(cs[i:i + 3], "X") for i in range(0, len(cs) - 3, 3))
            fh.write(f">{f}|{g}|{a}\n{p}\n")


def canon_pair(x, y):
    p1, p2 = (x, y), ((y[0], -y[1]), (x[0], -x[1]))
    return min(p1, p2)


def circ_pairs(lab):                   # lab: list of (family, strand) in genome order
    n = len(lab)
    return [canon_pair(lab[i], lab[(i + 1) % n]) for i in range(n)] if n >= 2 else []


for g, lst in fam_by.items():
    lab = [(f, s) for a, b, s, f in sorted(lst)]
    for p in set(circ_pairs(lab)):
        ref_pairs_count[p] += 1
REF_PAIRS = {p for p, c in ref_pairs_count.items() if c >= 0.5 * len(fam_by)}
print("calibration genomes:", len(fam_by), "reference gene-order pairs:", sorted(REF_PAIRS))

for g, seq in calib.items():
    gl = [x for x in genes_of(seq) if valid_orf(x)]
    gl.sort(key=lambda x: x["a"])
    for x in gl:
        codon_pairs[(x["cs"][:3], x["cs"][-3:])] += 1
    ncalg += len(gl)
    for u, v in zip(gl, gl[1:]):
        gaps.append(v["a"] - u["e"])
GAP_LO, GAP_HI = np.percentile(gaps, [5, 95])
OK_CODON = {k for k, c in codon_pairs.items() if c / ncalg >= 0.01}
print("calibration gap 5-95%:", GAP_LO, GAP_HI, "common start/stop pairs:", sorted(OK_CODON))

# ---- evaluate ---------------------------------------------------------------------------------------
seqs = {**read_fa("data/ladder.fna")}
import os
if os.path.exists("evo2/outputs.fna"):      # Evo 2 check (step 5): scored with exactly the same code
    seqs.update(read_fa("evo2/outputs.fna"))
for g, s in calib.items():
    seqs["CAL|" + g] = s
G = {n: genes_of(s) for n, s in seqs.items()}

# one blastp for all called proteins
with open("data/q_proteins.faa", "w") as fh:
    for n, gl in G.items():
        for i, x in enumerate(gl):
            if len(x["prot"]) >= 20:
                fh.write(f">{n}@{i}\n{x['prot']}\n")
subprocess.run("makeblastdb -in data/ref_proteins.faa -dbtype prot -out db/refprot >/dev/null", shell=True, check=True)
subprocess.run('blastp -query data/q_proteins.faa -db db/refprot -evalue 1e-5 -outfmt "6 qseqid sseqid evalue length qlen slen" '
               '-num_threads 8 -max_target_seqs 50 -out data/q_vs_ref.tsv', shell=True, check=True)
best = {}
for l in open("data/q_vs_ref.tsv"):
    q, s, ev, al, ql, sl = l.split()
    if int(al) >= 0.5 * min(int(ql), int(sl)):
        ev = float(ev)
        if q not in best or ev < best[q][1]:
            best[q] = (s.split("|")[0], ev)

rows = []
for n, seq in seqs.items():
    gl = G[n]; L = len(seq)
    run = max(len(m.group()) for m in re.finditer(r"(.)\1*", seq))
    ends = 0
    for k in range(8, min(2000, L // 2) + 1):
        if seq[:k] == seq[-k:]:
            ends = k
    cov = np.zeros(L, bool)
    for x in gl:
        cov[x["a"]:x["e"]] = True
    labs = []
    for i, x in enumerate(gl):
        x["fam"] = best.get(f"{n}@{i}", (None,))[0]
        if x["fam"]:
            labs.append((x["a"], x["fam"], x["st"]))
    labs.sort()
    pairs = circ_pairs([(f, s) for a, f, s in labs])
    srt = sorted(gl, key=lambda x: x["a"])
    gp = [v["a"] - u["e"] for u, v in zip(srt, srt[1:])]
    up = []
    for x in gl:
        u = seq[max(0, x["a"] - 20):max(0, x["a"] - 4)] if x["st"] == 1 else rc(seq[x["e"] + 4:x["e"] + 20])
        up.append(bool(re.search("GGAG|GAGG", u)))
    rows.append(dict(
        id=n, length=L, m1_gc=round((seq.count("G") + seq.count("C")) / L, 4), m1_nonacgt=sum(c not in "ACGT" for c in seq),
        m1_maxrun=run, m3_end_overlap=ends,
        m4_genes=len(gl), m4_density=round(cov.mean(), 4), m4_meanlen=round(np.mean([x["e"] - x["a"] for x in gl]), 1) if gl else 0,
        m6_core=len({x["fam"] for x in gl if x["fam"]}),
        m7_valid=round(np.mean([valid_orf(x) for x in gl]), 4) if gl else 0,
        m8_order=round(np.mean([p in REF_PAIRS for p in pairs]), 4) if pairs else 0,
        m9_gap=round(np.mean([GAP_LO <= g <= GAP_HI for g in gp]), 4) if gp else 0,
        m9_codon=round(np.mean([(x["cs"][:3], x["cs"][-3:]) in OK_CODON for x in gl]), 4) if gl else 0,
        m10_rbs=round(np.mean(up), 4) if up else 0))
with open("data/metrics_quick.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("rows:", len(rows))
