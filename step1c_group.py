"""step1c_group.py - rebuild protein families from data/blastp_all.tsv, save them, and describe the biggest ones
(how many genomes, typical protein length, which co-occur) so we can choose one coherent Microviridae group."""
import collections, statistics, sys
COV = float(sys.argv[1])
prot = {}
name = None
for l in open("data/rep_proteins.faa"):
    if l.startswith(">"): name = l[1:].strip()
    else: prot[name] = len(l.strip())
parent = {p: p for p in prot}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for l in open("data/blastp_all.tsv"):
    q, s, al, ql, sl = l.split(); al, ql, sl = int(al), int(ql), int(sl)
    if q != s and al >= COV * max(ql, sl) and min(ql, sl) >= 60: parent[find(q)] = find(s)
fam = collections.defaultdict(list)
for p in prot: fam[find(p)].append(p)
order = sorted(fam, key=lambda r: len({p.split("|")[0] for p in fam[r]}), reverse=True)
with open("data/families.tsv", "w") as fh:
    fh.write("family\tprotein\tlength\n")
    for i, r in enumerate(order[:40]):
        for p in fam[r]: fh.write(f"F{i+1}\t{p}\t{prot[p]}\n")
gsets = [{p.split("|")[0] for p in fam[r]} for r in order[:12]]
print("fam  genomes  proteins  median_len  p10-p90_len")
for i, r in enumerate(order[:12]):
    L = sorted(prot[p] for p in fam[r])
    print(f"F{i+1:<3} {len(gsets[i]):7d} {len(fam[r]):8d} {statistics.median(L):10.0f}  {L[len(L)//10]}-{L[9*len(L)//10]}")
print("genomes containing F1 and Fj:", [(f"F{j+1}", len(gsets[0] & gsets[j])) for j in range(1, 12)])
print("genomes with >=3 of F1..F6:", sum(1 for g in set().union(*gsets[:6]) if sum(g in s for s in gsets[:6]) >= 3))
print("genomes with >=4 of F1..F8:", sum(1 for g in set().union(*gsets[:8]) if sum(g in s for s in gsets[:8]) >= 4))
