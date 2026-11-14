"""gen_exp1.py - Evo 2 7B through EvSpark, continuing each prompt to its natural genome's length (design section 6).
Reuses the tested plumbing of ../kmer-evo2/generate_evspark.py (baseline condition = 1 candidate per round, no k-mers).
GPU run (Linux, ~40 GB VRAM; NOT YET RUN on a GPU):
  PYTHONPATH=/path/EvSpark:/path/EvSpark/scripts:/path/kmer-evo2 python gen_exp1.py --out evo2/out.fa \
      --checkpoint /path/L27_g12_150M_s1.pt --target-weights /path/evo2_7b.pt
CPU check of the plumbing: add --mock --max-n 2.   Resumable: finished records are skipped."""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "kmer-evo2"))
from kmerlib import KmerTable, read_fasta, write_fasta
import generate_evspark as ge

ap = argparse.ArgumentParser()
ap.add_argument("--prompts", default="evo2/prompts.fa"); ap.add_argument("--targets", default="evo2/targets.tsv")
ap.add_argument("--out", required=True); ap.add_argument("--chunk-len", type=int, default=512)
ap.add_argument("--checkpoint", default="L27_g12_150M_s1"); ap.add_argument("--target-weights", default="models/evo2_7b/evo2_7b.pt")
ap.add_argument("--gamma", type=int, default=12); ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--max-n", type=int, default=None); ap.add_argument("--mock", action="store_true")
a = ap.parse_args(); a.max_length = None
todo = {l.split("\t")[0]: int(l.split("\t")[4]) for l in list(open(a.targets))[1:]}
prompts = list(read_fasta(a.prompts))[: a.max_n]
records = list(read_fasta(a.out)) if os.path.exists(a.out) else []
done = {n for n, _ in records}
gen = None
if not a.mock:
    a.max_length = max(len(p) for _, p in prompts) + max(todo.values()) + 2 * a.chunk_len + 4 * a.gamma + 256
    gen = ge.load_real(a, 1)
for i, (sid, prompt) in enumerate(prompts):
    if sid in done: continue
    g = gen if gen is not None else ge.MockGenerator(KmerTable([6]).finalize())
    t0 = time.time()
    seq, rounds = ge.grow(g, None, prompt, todo[sid], a.chunk_len, 1, "paper", a.seed * 1_000_003 + i)
    dt = time.time() - t0
    records.append((sid, seq)); write_fasta(a.out, records)
    with open(a.out + ".log.jsonl", "a") as fh:
        fh.write(json.dumps({"id": sid, "letters": todo[sid], "seconds": dt, "letters_per_sec": todo[sid] / dt, "engine": "mock" if a.mock else "evspark"}) + "\n")
    print(f"{sid}: {todo[sid]} letters in {dt:.1f}s ({todo[sid] / dt:.0f}/s)", flush=True)
