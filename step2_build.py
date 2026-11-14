"""step2_build.py - orient the chosen Microviridae genomes, split 30 calibration / 60 test, build the damage ladder.
Run in WSL: ~/kmer-venv/bin/python step2_build.py      (pyrodigal; everything is deterministic, seed 1337)

Plain English
  * Orientation: genes are read on one strand in these viruses. We flip each genome so the family-F2 gene (the
    306-letter protein found in nearly every genome; not yet named) reads forward, then rotate the circle so
    that gene starts at position 1. F2 is chosen for being universal, not because it is proven to be the
    replication gene (design doc says "replication gene"; see RUNLOG).
  * Damage D1-D10 follows design doc section 4. Genes = pyrodigal meta-mode genes on the oriented genome.
    Genes cut by the new circle origin are not seen by the caller; they are simply not edited.
"""
import collections, random, itertools
import numpy as np, pyrodigal

SEED = 1337
import sys
MING = int(sys.argv[1]) if len(sys.argv) > 1 else 6   # minimum genes (>=150 letters) per genome
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
BASES = "ACGT"
CODON = {}
for i, a in enumerate("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    CODON["".join("TCAG"[x] for x in (i // 16, (i // 4) % 4, i % 4))] = a
SYN = collections.defaultdict(list)
for c, a in CODON.items():
    SYN[a].append(c)


def read_fa(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        else:
            d[n].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


def write_fa(p, items):
    with open(p, "w") as fh:
        for n, s in items:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


# ---- families (F1..F5 from step1c, cov 0.7) -----------------------------------------------------
prot = read_fa("data/rep_proteins.faa")
fam_of_prot = {}
members = collections.defaultdict(set)
for l in list(open("data/families.tsv"))[1:]:
    f, p, _ = l.split()
    if f in ("F1", "F2", "F3", "F4", "F5"):
        fam_of_prot[p] = f
        members[p.split("|")[0]].add(f)
seq2fam = {}
for p, f in fam_of_prot.items():
    seq2fam.setdefault(prot[p], f)
group = sorted(g for g, s in members.items() if len(s) == 5)
print("clusters with all five families:", len(group))

genomes = read_fa("data/cluster_reps.fna")
finder = pyrodigal.GeneFinder(meta=True)


# 6-letter protein words -> family, from the F1-F5 proteins of step 1 (used to name each gene's family)
kidx = collections.defaultdict(collections.Counter)
for pid, f in fam_of_prot.items():
    q = prot[pid]
    for i in range(len(q) - 5):
        kidx[q[i:i + 6]][f] += 1


def family_of(q):
    votes = collections.Counter()
    for i in range(len(q) - 5):
        for f in kidx.get(q[i:i + 6], ()):
            votes[f] += 1
    if not votes:
        return None
    f, n = votes.most_common(1)[0]
    return f if n >= max(3, 0.2 * (len(q) - 5)) else None


def call_circular(s):
    """genomes are circular: call genes on the doubled sequence, keep genes that START in the first copy.
    Returns (start, end, strand, protein); end may exceed len(s) when a gene runs over the origin."""
    L, out = len(s), []
    for g in finder.find_genes((s + s).encode()):
        a, b = g.begin - 1, g.end
        if a >= L or g.partial_begin or g.partial_end or (b - a) % 3 or b - a > L:
            continue
        out.append((a, b, g.strand, g.translate().rstrip("*")))
    return out


oriented, dropped = {}, collections.Counter()
for gname in group:
    s = genomes[gname]; L = len(s)
    gl = [(a, b, st, family_of(q)) for a, b, st, q in call_circular(s)]
    f2 = [x for x in gl if x[3] == "F2"]
    if not f2:
        dropped["no F2 gene"] += 1; continue
    if f2[0][2] == -1:                                        # flip so F2 reads forward
        s = rc(s)
        gl = [((-b) % L, (-b) % L + (b - a), -st, f) for a, b, st, f in gl]
    f2 = [x for x in gl if x[3] == "F2"][0]
    r = f2[0] % L
    s = s[r:] + s[:r]                                         # rotate: F2 starts at position 1
    genes = []
    for a, b, st, f in gl:
        a2 = (a - r) % L; b2 = a2 + (b - a)
        if b2 > L:
            continue                                          # runs over the new origin: left unedited
        cs = s[a2:b2] if st == 1 else rc(s[a2:b2])
        if cs[:3] not in ("ATG", "GTG", "TTG") or CODON.get(cs[-3:]) != "*" or any(CODON[cs[i:i + 3]] == "*" for i in range(0, len(cs) - 3, 3)):
            dropped["gene failed start/stop check (gene dropped, genome kept)"] += 1; continue
        genes.append((a2, b2, st, f))
    if not all(any(g[3] == f for g in genes) for f in ("F1", "F2", "F3", "F4", "F5")):
        dropped["lost a core gene"] += 1; continue
    if not any(g[0] == 0 and g[3] == "F2" for g in genes):
        dropped["F2 not at position 1"] += 1; continue
    if sum(g[1] - g[0] >= 150 for g in genes) < MING:
        dropped["fewer than 8 callable genes"] += 1; continue
    oriented[gname] = (s, sorted(genes))
print("oriented genomes:", len(oriented), "dropped:", dict(dropped))

names = sorted(oriented)
random.Random(SEED).shuffle(names)
calib, test = sorted(names[:30]), sorted(names[30:90])
with open("data/split.tsv", "w") as fh:
    fh.write("genome\tset\tlength\tn_genes\n")
    for n in calib + test:
        fh.write(f"{n}\t{'calibration' if n in calib else 'test'}\t{len(oriented[n][0])}\t{len(oriented[n][1])}\n")
with open("data/genes_oriented.tsv", "w") as fh:
    fh.write("genome\tset\tgene\tstart0\tend\tstrand\tfamily\tprotein_len\n")
    for n in calib + test:
        for i, (a, b, st, f) in enumerate(oriented[n][1]):
            fh.write(f"{n}\t{'calibration' if n in calib else 'test'}\t{i}\t{a}\t{b}\t{st}\t{f or ''}\t{(b - a) // 3 - 1}\n")
write_fa("data/natural_calib.fna", [(n, oriented[n][0]) for n in calib])
write_fa("data/natural_test.fna", [(n, oriented[n][0]) for n in test])


# ---- damage functions: each takes (seq, genes, rng, severity) -> (new_seq, note) ----------------------
def cover(seq, genes):
    cov = np.zeros(len(seq), np.int16)
    for a, b, st, p in genes:
        cov[a:b] += 1
    return cov


def coding(arr, a, b, st):
    s = "".join(arr[a:b])
    return s if st == 1 else rc(s)


def put(arr, a, b, st, s):
    arr[a:b] = list(s if st == 1 else rc(s))


def d1(seq, genes, rng, p):
    arr, cov = list(seq), cover(seq, genes)
    cand = []                              # (gene, codon index): has a synonym and touches no overlapped letter
    for gi, (a, b, st, _) in enumerate(genes):
        cs = coding(arr, a, b, st)
        n = len(cs) // 3
        for ci in range(1, n - 1):
            c = cs[3 * ci:3 * ci + 3]
            lo = a + 3 * ci if st == 1 else b - 3 * ci - 3
            if len(SYN.get(CODON.get(c), [])) > 1 and cov[lo:lo + 3].max() == 1:
                cand.append((gi, ci))
    k = round(p * len(cand))
    by = collections.defaultdict(list)
    for gi, ci in rng.sample(cand, k):
        by[gi].append(ci)
    for gi, cis in by.items():
        a, b, st, _ = genes[gi]
        cs = list(coding(arr, a, b, st))
        for ci in cis:
            c = "".join(cs[3 * ci:3 * ci + 3])
            cs[3 * ci:3 * ci + 3] = list(rng.choice([x for x in SYN[CODON[c]] if x != c]))
        put(arr, a, b, st, "".join(cs))
    return "".join(arr), f"{k} codons swapped of {len(cand)} swappable"


def d2(seq, genes, rng, rate):
    n = round(rate / 100 * len(seq))
    arr = list(seq)
    for i in rng.sample(range(len(seq)), n):
        arr[i] = rng.choice([x for x in BASES if x != arr[i]])
    return "".join(arr), f"{n} substitutions"


def d3(seq, genes, rng, n):
    cov = cover(seq, genes)
    inside = [i for i in range(1, len(seq) - 1) if cov[i] > 0]
    arr, ops = list(seq), []
    for i in sorted(rng.sample(inside, n), reverse=True):
        if rng.random() < 0.5:
            ops.append(f"ins@{i}"); arr.insert(i, rng.choice(BASES))
        else:
            ops.append(f"del@{i}"); del arr[i]
    return "".join(arr), ",".join(ops)


def d4(seq, genes, rng, k):
    arr, ops = list(seq), []
    ok = [i for i, g in enumerate(genes) if g[1] - g[0] >= 150]
    k = min(k, len(ok) - 2)                                   # genomes have only 6-9 callable genes: always leave 2
    for gi in rng.sample(ok, k):
        a, b, st, _ = genes[gi]
        cs = list(coding(arr, a, b, st))
        n = len(cs) // 3
        if rng.random() < 0.5:
            ci = rng.randrange(1, n - 1)
            cs[3 * ci:3 * ci + 3] = list("TAA"); ops.append(f"g{gi}:stop@codon{ci}")
        else:
            cs[0:3] = list("GCT"); ops.append(f"g{gi}:start-removed")
        put(arr, a, b, st, "".join(cs))
    return "".join(arr), f"k={k} of nominal; " + ",".join(ops)


def d5(seq, genes, rng, k):
    k = min(k, len(genes) - 2)
    pick = rng.sample(range(len(genes)), k)
    merged = []
    for a, b in sorted((genes[i][0], genes[i][1]) for i in pick):
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    keep, last = [], 0
    for a, b in merged:
        keep.append(seq[last:a]); last = b
    keep.append(seq[last:])
    return "".join(keep), f"k={k}; genes {sorted(pick)} removed ({sum(b - a for a, b in merged)} letters)"


def d6(seq, genes, rng, n):
    L = len(seq)
    for _ in range(200):
        cuts = sorted(rng.sample(range(1, L), n))
        pieces = [seq[cuts[i]:cuts[i + 1]] for i in range(n - 1)] + [seq[cuts[-1]:] + seq[:cuts[0]]]
        order = list(range(n)); rng.shuffle(order)
        flips = [rng.random() < 0.5 for _ in range(n)]
        out = "".join(rc(pieces[o]) if f else pieces[o] for o, f in zip(order, flips))
        if out not in seq + seq and out not in rc(seq) + rc(seq):
            return out, f"{n} cuts at {cuts}, order {order}, flips {[int(f) for f in flips]}"
    raise RuntimeError("could not make a non-trivial shuffle")


def d7(seq, genes, rng, frac):
    arr, tot = list(seq), 0
    for a, b, st, _ in genes:                 # genes are processed one after another; overlaps see earlier edits
        cs = coding(arr, a, b, st)
        n = len(cs) // 3
        cod = [cs[3 * i:3 * i + 3] for i in range(n)]
        idx = rng.sample(range(1, n - 1), round(frac * (n - 2)))
        vals = [cod[i] for i in idx]; rng.shuffle(vals)
        for i, v in zip(idx, vals):
            cod[i] = v
        tot += len(idx)
        put(arr, a, b, st, "".join(cod))
    return "".join(arr), f"{tot} codons shuffled within genes (start and stop kept)"


def d10(seq, genes, rng, which):
    if which == "revcomp":
        return rc(seq), "reverse complement"
    off = rng.randrange(len(seq) // 10, len(seq) * 9 // 10)
    return seq[off:] + seq[:off], f"rotated by {off}"


# k-mer-only chain (D8): order-5 chain (6-letter words) learned from the CALIBRATION genomes only
cnt = np.full((1024, 4), 0.1)
for n in calib:
    v = ["ACGT".index(c) for c in oriented[n][0]]
    for i in range(5, len(v)):
        ctx = 0
        for x in v[i - 5:i]:
            ctx = ctx * 4 + x
        cnt[ctx, v[i]] += 1
prob = cnt / cnt.sum(1, keepdims=True)


def d8(seq, genes, rng, _):
    r = np.random.default_rng(rng.randrange(2 ** 32))
    v = [rng.randrange(4) for _i in range(5)]
    ctx = 0
    for x in v:
        ctx = ctx * 4 + x
    for _i in range(len(seq) - 5):
        x = int(r.choice(4, p=prob[ctx])); v.append(x); ctx = (ctx * 4 + x) % 1024
    return "".join(BASES[x] for x in v), "order-5 chain from calibration genomes"


# unrelated windows (D9): alternate 8 well-known phages and 3 bacteria (not E. coli)
phages = list(read_fa("data/phages/ncbi_dataset/data/genomic.fna").items())
bact = []
for f in ("family/GCF_000006945.2", "family/GCF_000240185.1", "wrong/GCF_000006765.1"):
    recs = read_fa(f"../kmer-evo2/data/genomes/{f}.fna")
    bact.append((f.split("/")[1], max(recs.values(), key=len)))


def d9_factory():
    counter = itertools.count()

    def d9(seq, genes, rng, _):
        pool = phages if next(counter) % 2 == 0 else bact
        for _t in range(100):
            nm, s = rng.choice(pool)
            L = rng.randint(5000, 6000); st = rng.randrange(len(s) - L)
            w = s[st:st + L]
            if set(w) <= set(BASES):
                return w, f"window of {nm.split()[0]} at {st}"
        raise RuntimeError("no clean window")
    return d9


LADDER = [("D1", d1, [10, 30]), ("D2", d2, [0.5, 1, 2, 5, 10]), ("D3", d3, [1, 3, 6, 12]), ("D4", d4, [1, 3, 6]),
          ("D5", d5, [1, 3, 6]), ("D6", d6, [2, 4, 8]), ("D7", d7, [25, 50, 100]), ("D8", d8, [1]),
          ("D9", d9_factory(), [1]), ("D10", d10, ["revcomp", "rotation"])]
out, rows = [], []
for n in test:
    seq, genes = oriented[n]
    out.append((n + "__natural", seq)); rows.append((n + "__natural", n, "natural", "-", len(seq), ""))
    for did, fn, sevs in LADDER:
        for sv in sevs:
            sev = sv / 100 if did in ("D1", "D7") else sv
            rng = random.Random(f"{SEED}|{n}|{did}|{sv}")
            s2, note = fn(seq, genes, rng, sev)
            assert set(s2) <= set(BASES), (n, did, sv)
            sid = f"{n}__{did}_{sv}"
            out.append((sid, s2)); rows.append((sid, n, did, sv, len(s2), note))
write_fa("data/ladder.fna", out)
with open("data/ladder_manifest.tsv", "w") as fh:
    fh.write("id\tparent\tdamage\tseverity\tlength\tnote\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print("ladder sequences:", len(out), dict(collections.Counter(r[2] for r in rows)))
