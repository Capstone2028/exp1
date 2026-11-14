"""step7_standardise.py - POST-HOC robustness check (not part of the frozen analysis): give every circular sequence the same start.

Why: M4 and M8 reacted to D10 rotation (AUROC 0.62 / 0.60) because a linear gene finder cuts genes at the sequence ends,
so where the circle is opened changes what it sees. Fix: before scoring, rotate each sequence so it starts at the
start codon of its family-F2 gene (left end if forward, right end if reverse strand; only when exactly one intact F2 gene, 0.8-1.25x the calibration median length, is found) (the anchor step 2 used to start every natural genome). Genomes that already
start there are unchanged; sequences with no recognisable F2 gene (e.g. k-mer chain, wrong family) are left as they are.
Applied identically to every damage type; the ORIGINAL results are kept and reported next to the standardised ones.

Run in WSL:  ~/kmer-venv/bin/python step7_standardise.py <exp1 dir> <out dir>
Writes <out>/data/ladder.fna (ONLY the sequences whose start moved) and <out>/data/std_changes.tsv (all ladder ids, offset, changed?)."""
import collections, os, sys
import pyrodigal

EXP, OUT = sys.argv[1], sys.argv[2]
os.makedirs(f"{OUT}/data", exist_ok=True)


def read_fa(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        else:
            d[n].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


prot = read_fa(f"{EXP}/data/rep_proteins.faa")
fam_of_prot = {}
for l in list(open(f"{EXP}/data/families.tsv"))[1:]:
    f, p, _ = l.split()
    if f in ("F1", "F2", "F3", "F4", "F5"):
        fam_of_prot[p] = f
kidx = collections.defaultdict(collections.Counter)
for pid, f in fam_of_prot.items():
    q = prot[pid]
    for i in range(len(q) - 5):
        kidx[q[i:i + 6]][f] += 1


def family_of(q):                       # identical rule to step2_build.py
    votes = collections.Counter()
    for i in range(len(q) - 5):
        for f in kidx.get(q[i:i + 6], ()):
            votes[f] += 1
    if not votes:
        return None
    f, n = votes.most_common(1)[0]
    return f if n >= max(3, 0.2 * (len(q) - 5)) else None


import statistics
F2LEN = statistics.median(3 * (int(l.split('	')[7]) + 1) for l in list(open(f"{EXP}/data/genes_oriented.tsv"))[1:] if l.split('	')[6] == 'F2' and l.split('	')[1] == 'calibration')
print('median F2 length (letters):', F2LEN)
finder = pyrodigal.GeneFinder(meta=True)
lad = read_fa(f"{EXP}/data/ladder.fna")
changed, rows = {}, []
for sid, s in lad.items():
    L = len(s)
    f2 = []
    for g in finder.find_genes((s + s + s).encode()):             # circular: call on three copies, keep genes starting in the middle one
        a, b = g.begin - 1 - L, g.end - L
        if not (0 <= a < L) or g.partial_begin or g.partial_end or b - a > L:
            continue
        if family_of(g.translate().rstrip("*")) == "F2":
            f2.append((a, b, g.strand))
    best = f2[0] if len(f2) == 1 and 0.8 * F2LEN <= f2[0][1] - f2[0][0] <= 1.25 * F2LEN else None   # exactly one intact F2 gene
    # forward gene: its start codon (left end) goes to position 1; reverse-strand gene: its start codon is the RIGHT end,
    # so the gene's end goes to the end of the sequence (a reverse complement of a standard genome stays as it is)
    off = ((best[0] if best[2] == 1 else best[1]) % L) if best else 0
    rows.append((sid, off if best else "noF2", int(off != 0)))
    if off != 0:
        changed[sid] = s[off:] + s[:off]
with open(f"{OUT}/data/ladder.fna", "w", newline="\n") as fh:
    for k, v in changed.items():
        fh.write(f">{k}\n{v}\n")
with open(f"{OUT}/data/std_changes.tsv", "w") as fh:
    fh.write("id\toffset\tchanged\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
cnt = collections.Counter(r[0].split("__")[1].split("_")[0] for r in rows if r[2])
tot = collections.Counter(r[0].split("__")[1].split("_")[0] for r in rows)
print("sequences whose start moved, by damage type:", {k: f"{cnt[k]}/{tot[k]}" for k in tot})
