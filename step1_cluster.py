"""step1_cluster.py - keep complete Microviridae genomes of 4-6.5 kb, remove near-copies (95% identity).

Plain English: two genomes that are >=95% identical count as one "cluster", and we keep one genome per
cluster so near-copies are not counted many times. Identity is estimated from shared 16-letter words
(a "sketch": only words whose hash is divisible by 8 are kept, which is a fair random sample). Words are
counted on both strands, so a reversed or rotated genome still matches. A sketch estimate is not BLAST,
so the threshold is approximate; see RUNLOG.md.
  shared fraction f of 16-mers ~ identity^16, Jaccard j = f/(2-f); identity 0.95 -> f 0.44 -> j 0.28
"""
import re, sys, json, collections
import numpy as np

FNA = "data/mv/ncbi_dataset/data/genomic.fna"
K, MOD, JMIN = 16, 8, 0.28
LUT = np.full(256, 4, np.uint8)
for i, c in enumerate("ACGT"):
    LUT[ord(c)] = i

def read_fasta(path):
    name, hdr, parts = None, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, hdr, "".join(parts).upper()
            hdr = line[1:].strip(); name = hdr.split()[0]; parts = []
        else:
            parts.append(line.strip())
    if name:
        yield name, hdr, "".join(parts).upper()

def sketch(seq):
    a = LUT[np.frombuffer(seq.encode(), np.uint8)].astype(np.uint64)
    n = len(a) - K + 1
    f = np.zeros(n, np.uint64); r = np.zeros(n, np.uint64)
    for j in range(K):
        seg = a[j:j + n]
        f = f * np.uint64(4) + seg
        r = r + (np.uint64(3) - seg) * (np.uint64(4) ** np.uint64(j))
    can = np.minimum(f, r)
    h = (can * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(7)
    return set(h[(h % np.uint64(MOD)) == 0].tolist())

recs, drop = [], collections.Counter()
for name, hdr, seq in read_fasta(FNA):
    if not (4000 <= len(seq) <= 6500):
        drop["length"] += 1; continue
    if re.search(r"[^ACGT]", seq):
        drop["non-ACGT"] += 1; continue
    recs.append((name, hdr, seq))
print("kept", len(recs), "dropped", dict(drop))

recs.sort(key=lambda r: (-len(r[2]), r[0]))        # deterministic: longest first, then accession
index = collections.defaultdict(list)             # hash -> representative ids
reps, sizes, members = [], [], []
for name, hdr, seq in recs:
    s = sketch(seq)
    hits = collections.Counter()
    for h in s:
        for rid in index.get(h, ()):
            hits[rid] += 1
    best, bj = None, 0.0
    for rid, sh in hits.items():
        j = sh / (len(s) + len(reps[rid][3]) - sh)
        if j > bj:
            best, bj = rid, j
    if best is not None and bj >= JMIN:
        members[best].append(name)
    else:
        rid = len(reps)
        reps.append((name, hdr, seq, s)); members.append([name])
        for h in s:
            index[h].append(rid)

with open("data/clusters.tsv", "w") as fh:
    fh.write("rep\tn_members\tlength\theader\n")
    for (name, hdr, seq, _), m in zip(reps, members):
        fh.write(f"{name}\t{len(m)}\t{len(seq)}\t{hdr}\n")
with open("data/cluster_reps.fna", "w") as fh:
    for name, hdr, seq, _ in reps:
        fh.write(f">{name} {hdr.split(' ', 1)[-1]}\n{seq}\n")
print("clusters", len(reps), "largest", sorted(map(len, members))[-5:])
