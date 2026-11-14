"""step3c_phanotate.py - M5: PHANOTATE gene calls (the gene finder inside Pharokka), same three numbers as M4:
gene count, coding density, mean gene length. Run in WSL: ~/kmer-venv/bin/python step3c_phanotate.py
Each sequence is run on its own (PHANOTATE treats it as one linear contig). Output: data/metrics_phanotate.csv"""
import csv, os, subprocess, tempfile
from multiprocessing import Pool
import numpy as np

EXE = os.path.expanduser("~/kmer-venv/bin/phanotate.py")

def read_fa(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        else: d[n].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}

def run(item):
    n, seq = item
    with tempfile.NamedTemporaryFile("w", suffix=".fna", delete=False) as fh:
        fh.write(f">s\n{seq}\n"); path = fh.name
    out = subprocess.run([EXE, "-f", "tabular", path], capture_output=True, text=True).stdout
    os.unlink(path)
    cov = np.zeros(len(seq), bool); lens = []
    for l in out.splitlines():
        if l.startswith("#") or not l.strip(): continue
        a, b = sorted(int(x) for x in l.split("\t")[:2])
        cov[a - 1:b] = True; lens.append(b - a + 1)
    return dict(id=n, m5_genes=len(lens), m5_density=round(cov.mean(), 4), m5_meanlen=round(float(np.mean(lens)), 1) if lens else 0)

if __name__ == "__main__":
    seqs = read_fa("data/ladder.fna")
    for g, s in read_fa("data/natural_calib.fna").items(): seqs["CAL|" + g] = s
    if os.path.exists("evo2/outputs.fna"): seqs.update(read_fa("evo2/outputs.fna"))
    with Pool(8) as p: rows = p.map(run, list(seqs.items()), chunksize=8)
    with open("data/metrics_phanotate.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("rows", len(rows))
