"""step1b_families.py - call genes in every cluster representative, group similar proteins into families.
Run in WSL (pyrodigal): ~/kmer-venv/bin/python step1b_families.py
A 'protein family' = proteins linked by blastp hits (e<=1e-5, >=50% of the shorter protein aligned)."""
import collections, subprocess, pyrodigal
fin = pyrodigal.GeneFinder(meta=True)
genomes = {}
name = None
for l in open("data/cluster_reps.fna"):
    if l.startswith(">"): name = l[1:].split()[0]; genomes[name] = []
    else: genomes[name].append(l.strip())
prot, ngenes = {}, {}
with open("data/rep_proteins.faa", "w") as fh:
    for g, parts in genomes.items():
        seq = "".join(parts)
        genes = fin.find_genes(seq.encode())
        ngenes[g] = len(genes)
        for i, gene in enumerate(genes):
            pid = f"{g}|{i}"; p = gene.translate().rstrip("*"); prot[pid] = p
            fh.write(f">{pid}\n{p}\n")
print("genomes", len(genomes), "proteins", len(prot))
subprocess.run("makeblastdb -in data/rep_proteins.faa -dbtype prot -out db/repprot >/dev/null", shell=True, check=True)
subprocess.run('blastp -query data/rep_proteins.faa -db db/repprot -evalue 1e-5 -outfmt "6 qseqid sseqid length qlen slen" '
               '-num_threads 8 -max_target_seqs 2000 -out data/blastp_all.tsv', shell=True, check=True)
parent = {p: p for p in prot}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for l in open("data/blastp_all.tsv"):
    q, s, al, ql, sl = l.split(); al, ql, sl = int(al), int(ql), int(sl)
    if q != s and al >= 0.5 * min(ql, sl): parent[find(q)] = find(s)
fam = collections.defaultdict(set)
for p in prot: fam[find(p)].add(p.split("|")[0])
top = sorted(fam.values(), key=len, reverse=True)[:15]
print("top protein families by number of genomes containing them:", [len(t) for t in top])
import statistics
print("genes per genome: median", statistics.median(ngenes.values()))
