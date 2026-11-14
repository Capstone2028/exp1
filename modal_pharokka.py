"""modal_pharokka.py - run the real Pharokka (PHANOTATE genes + PHROG profile search) on Modal, so M12 (and M5) can be done "properly".
UNTESTED until first run. From the exp1 folder, with the kmer-evo2 venv:
  ../kmer-evo2/.venv/Scripts/modal run modal_pharokka.py::fetch_db                                  # CPU: Pharokka database into a Volume
  ../kmer-evo2/.venv/Scripts/modal run modal_pharokka.py::run_pharokka --inp pharokka_test.fna --tag test     # 3-sequence smoke test: time + peak memory
  ../kmer-evo2/.venv/Scripts/modal run modal_pharokka.py::run_pharokka --inp pharokka_in.fna --tag all        # all 1,804 sequences
  ../kmer-evo2/.venv/Scripts/modal volume get pharokka-vol out/all ./pharokka/out_all
Plain English: Pharokka needs a lot of memory for the PHROG profile search (it ran out at 7 GB on the local WSL). Here it gets a big
CPU machine. Pharokka is installed from conda (bioconda) INSIDE the Modal image, so nothing is installed locally.
Pharokka is run with --skip_extra_annotations --skip_mash: tRNAscan-SE, MinCED, ARAGORN and the mash genome match do not feed M5 or M12.
Cost guard: hard timeouts; each function prints seconds, peak memory and approximate dollars."""
import os
import subprocess
import time

import modal

HERE = os.path.dirname(os.path.abspath(__file__))
CPUS, MEM_GB = 16, 48
CPU_PER_H = CPUS * 0.047 + MEM_GB * 0.008        # $/h at Modal's CPU and memory rates as read 25 Oct 2026

app = modal.App("exp1-pharokka")
vol = modal.Volume.from_name("pharokka-vol", create_if_missing=True)
V = "/vol"

image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install("pharokka", channels=["conda-forge", "bioconda"])
    .add_local_file(os.path.join(HERE, "pharokka", "pharokka_test.fna"), "/app/pharokka_test.fna")
    .add_local_file(os.path.join(HERE, "pharokka", "pharokka_in.fna"), "/app/pharokka_in.fna")
    .add_local_file(os.path.join(HERE, "pharokka", "pharokka_core.fna"), "/app/pharokka_core.fna")
    .add_local_file(os.path.join(HERE, "pharokka", "pharokka_test12.fna"), "/app/pharokka_test12.fna")
)


def _sh(cmd, timeout):
    print("$", cmd, flush=True)
    r = subprocess.run(cmd, shell=True, timeout=timeout)
    if r.returncode:
        raise RuntimeError(f"failed ({r.returncode}): {cmd}")


@app.function(image=image, volumes={V: vol}, cpu=4, memory=16384, timeout=3600)
def fetch_db():
    t0 = time.time()
    _sh("pharokka.py --version || true", 120)
    _sh(f"(pharokka install -o {V}/db || install_databases.py -o {V}/db) && du -sh {V}/db && ls {V}/db | head -30", 3000)
    vol.commit()
    print(f"[cost] {time.time() - t0:.0f} s CPU container (about $0.05)")


@app.function(image=image, volumes={V: vol}, cpu=CPUS, memory=MEM_GB * 1024, timeout=7200)
def run_pharokka(inp: str = "pharokka_test.fna", tag: str = "test"):
    import resource                                # Linux-only: imported inside the cloud function
    t0 = time.time()
    out = f"{V}/out/{tag}"
    _sh("pharokka.py -h | grep -iE 'skip_extra|skip_mash|--meta|--split|gene_predictor|reverse_mmseqs2|sensitivity' | head", 60)
    _sh(f"rm -rf {out}; mkdir -p {V}/out && pharokka.py -i /app/{inp} -o {out} -d {V}/db -t {CPUS} -m -s -g phanotate "
        f"--skip_extra_annotations --skip_mash -f", 6800)
    peak_gb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1e6   # kB -> GB (Linux)
    vol.commit()
    _sh(f"ls {out} | head -30; echo; du -sh {out}", 120)
    s = time.time() - t0
    print(f"[cost] {s:.0f} s on {CPUS} cores / {MEM_GB} GiB, peak memory of the largest child process about {peak_gb:.1f} GB = "
          f"about ${s / 3600 * CPU_PER_H:.2f} (check the Modal dashboard)")


@app.local_entrypoint()
def main():
    print(__doc__)
