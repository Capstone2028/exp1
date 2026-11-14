"""step12_prompt_excluded.py - recompute the Evo 2 protein-level numbers counting only genes that reach the GENERATED letters
(user correction: G1 prompts hand over intact natural sequence, so detections inside the copied prompt are not generation).
Run in WSL: ~/kmer-venv/bin/python step12_prompt_excluded.py        Output: evo2/prompt_excluded.csv
Three versions per Evo 2 output (prompt length P from evo2/targets.tsv):
  orig : all called genes (what was reported before)
  ov   : genes whose END is after letter P (they overlap generated letters; includes genes that start in the prompt)
  wh   : genes that START at or after letter P (wholly generated)
Metrics here (gene calls = pyrodigal meta, the same calls as M4/M6/M11):
  m6_core : number of the 5 core families (F1-F5) matched by a gene's protein (same blastp hits and rule as M6: data/q_vs_ref.tsv)
  pfam    : share of those genes with a Pfam hit (pyhmmer, gathering cutoffs), re-searched for the ov genes only
Phold (M13) and the simplified M12 on Phold's proteins are handled in step12b (Windows) from phold/out2."""
import collections, csv
import pyrodigal, pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile, TextSequence
from pyhmmer.plan7 import HMMFile

COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
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


ev = read_fa("evo2/outputs.fna")
PL = {l.split("\t")[0]: int(l.split("\t")[2]) for l in list(open("evo2/targets.tsv"))[1:]}
finder = pyrodigal.GeneFinder(meta=True)
best = {}                                              # same rule as step3_quick_metrics.py
for l in open("data/q_vs_ref.tsv"):
    q, s, ev_, al, ql, sl = l.split()
    if int(al) >= 0.5 * min(int(ql), int(sl)):
        e = float(ev_)
        if q not in best or e < best[q][1]:
            best[q] = (s.split("|")[0], e)
genes = {}
for sid, seq in ev.items():
    gl = []
    for i, g in enumerate(finder.find_genes(seq.encode())):
        a, e = g.begin - 1, g.end
        cs = seq[a:e] if g.strand == 1 else rc(seq[a:e])
        prot = "".join(CODON.get(cs[k:k + 3], "X") for k in range(0, len(cs) - 3, 3))
        gl.append(dict(i=i, a=a, e=e, prot=prot, fam=best.get(f"{sid}@{i}", (None,))[0]))
    genes[sid] = gl
# Pfam for genes that overlap generated letters (>= 20 aa, as in the main pipeline)
todo = {f"{sid}@{x['i']}": x["prot"] for sid, gl in genes.items() for x in gl if x["e"] > PL[sid] and len(x["prot"]) >= 20}
alpha = Alphabet.amino()
seqs = [TextSequence(name=k.encode(), sequence=v.replace("X", "A")).digitize(alpha) for k, v in todo.items()]
hits = set()
with HMMFile("db/pfam/Pfam-A.hmm") as hf:
    for th in pyhmmer.hmmsearch(hf, seqs, cpus=8, bit_cutoffs="gathering"):
        for h in th:
            if h.included:
                hits.add(h.name.decode() if isinstance(h.name, bytes) else h.name)
print("generated-overlapping proteins:", len(todo), "with Pfam hit:", len(hits))
rows = []
for sid, gl in genes.items():
    P = PL[sid]
    r = {"id": sid, "set": sid[-2:], "prompt_len": P, "genes_orig": len(gl)}
    for tag, keep in (("orig", gl), ("ov", [x for x in gl if x["e"] > P]), ("wh", [x for x in gl if x["a"] >= P])):
        r[f"genes_{tag}"] = len(keep)
        r[f"m6_{tag}"] = len({x["fam"] for x in keep if x["fam"]})
        if tag != "orig":
            kk = [x for x in keep if len(x["prot"]) >= 20]
            r[f"pfam_{tag}"] = round(sum(f"{sid}@{x['i']}" in hits for x in kk) / len(kk), 4) if kk else ""
    rows.append(r)
with open("evo2/prompt_excluded.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("rows", len(rows))
