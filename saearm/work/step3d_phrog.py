"""step3d_phrog.py - M12: share of called genes labelled by a PHROG family (phage protein families).
Run in WSL: ~/kmer-venv/bin/python step3d_phrog.py      Needs mmseqs (sudo apt install mmseqs2).
This is a simplified version of the PHROG search Pharokka does. Pharokka searches protein-vs-profile; that ran out of memory
(7 GB WSL) / took >15 min for 200 proteins here, so we search the PHROG profile sequence set instead (less sensitive).
Same database (MMseqs2 protein search against the PHROG family sequences, -s 7.5, e<=1e-5, Pharokka's PHROG v4
database from Zenodo 7563578) on the proteins pyrodigal called (data/q_proteins.faa). It is NOT the full Pharokka
pipeline (no PHANOTATE-based calls, no CRISPR/tRNA steps); genes are the same ones used for M4/M11.
A hit to a PHROG whose annotation is 'unknown function' still counts as a hit (it is a phage-protein family).
Output: data/metrics_phrog.csv (id, m12_genes, m12_phrog_hits, m12_phrog_share, m12_known_share)"""
import collections, csv, subprocess, os
DB = "db/pharokka/pharokka_v1.2.0_databases"
os.makedirs("data/mm", exist_ok=True)
sh = lambda c: subprocess.run(c, shell=True, check=True)
sh("mmseqs createdb data/q_proteins.faa data/mm/q >/dev/null")
sh(f"cp data/mm/q_h.dbtype {DB}/phrogs_profile_db_seq_h.dbtype")      # the download lacks this tiny header-type file
sh(f"mmseqs search data/mm/q {DB}/phrogs_profile_db_seq data/mm/res data/mm/tmp -e 1e-5 -s 7.5 --threads 8 >/dev/null")
sh(f'mmseqs convertalis data/mm/q {DB}/phrogs_profile_db_seq data/mm/res data/phrog_hits.tsv --format-output "query,target,evalue,bits" --threads 8 >/dev/null')
annot = {}
for l in list(open(f"{DB}/phrog_annot_v4.tsv"))[1:]:
    p = l.rstrip("\n").split("\t")
    annot[p[0]] = p[2] if len(p) > 2 else ""
best = {}
for l in open("data/phrog_hits.tsv"):
    q, t, ev, bits = l.split("\t")
    if q not in best or float(ev) < best[q][1]:
        best[q] = (t, float(ev))
tot, hit, known = collections.Counter(), collections.Counter(), collections.Counter()
for l in open("data/q_proteins.faa"):
    if l.startswith(">"):
        q = l[1:].strip(); sid = q.rsplit("@", 1)[0]
        tot[sid] += 1
        if q in best:
            hit[sid] += 1
            num = best[q][0].split("_")[-1]
            known[sid] += annot.get(num, "unknown function") not in ("unknown function", "", "NA", "hypothetical protein")
ids = [l.split(",")[0] for l in list(open("data/metrics_quick.csv"))[1:]]
with open("data/metrics_phrog.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["id", "m12_genes", "m12_phrog_hits", "m12_phrog_share", "m12_known_share"])
    for i in ids:
        w.writerow([i, tot[i], hit[i], round(hit[i] / tot[i], 4) if tot[i] else 0, round(known[i] / tot[i], 4) if tot[i] else 0])
print("proteins", sum(tot.values()), "with PHROG hit", sum(hit.values()))
