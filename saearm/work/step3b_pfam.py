"""step3b_pfam.py - M11: share of called genes whose protein has a Pfam family hit (pyhmmer, Pfam gathering cutoffs).
Run in WSL: ~/kmer-venv/bin/python step3b_pfam.py     Input: data/q_proteins.faa (from step3), db/pfam/Pfam-A.hmm
Output: data/metrics_pfam.csv (id, m11_genes, m11_pfam_hits, m11_pfam_share)
Plain English: Pfam is a library of protein-family fingerprints. A real gene's protein usually matches one;
a scrambled one usually does not. The "gathering cutoff" is the score each family's curators say is safe."""
import collections, csv
import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile

alpha = Alphabet.amino()
with SequenceFile("data/q_proteins.faa", digital=True, alphabet=alpha) as sf:
    seqs = sf.read_block()
print("proteins:", len(seqs))
hit_names = set()
with HMMFile("db/pfam/Pfam-A.hmm") as hf:
    for hits in pyhmmer.hmmsearch(hf, seqs, cpus=8, bit_cutoffs="gathering"):
        for h in hits:
            if h.included:
                hit_names.add(h.name.decode() if isinstance(h.name, bytes) else h.name)
print("proteins with a Pfam hit:", len(hit_names))
tot, hit = collections.Counter(), collections.Counter()
for s in seqs:
    name = s.name.decode() if isinstance(s.name, bytes) else s.name
    sid = name.rsplit("@", 1)[0]
    tot[sid] += 1
    hit[sid] += name in hit_names
ids = [l.split(",")[0] for l in list(open("data/metrics_quick.csv"))[1:]]
with open("data/metrics_pfam.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["id", "m11_genes", "m11_pfam_hits", "m11_pfam_share"])
    for i in ids:
        w.writerow([i, tot[i], hit[i], round(hit[i] / tot[i], 4) if tot[i] else 0])
