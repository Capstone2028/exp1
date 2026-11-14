"""modal_phold.py - run Phold (M13) on Modal for the 220 ladder sequences (20 genomes x 11 variants).
UNTESTED on Modal until first run. Run from the exp1 folder with the kmer-evo2 venv:
  ../kmer-evo2/.venv/Scripts/modal run modal_phold.py::install_db     # CPU only: downloads ProstT5 + 15 GB Phold DB into a Volume
  ../kmer-evo2/.venv/Scripts/modal run modal_phold.py::run_phold      # GPU (T4) for ProstT5 + CPU cores for Foldseek
  ../kmer-evo2/.venv/Scripts/modal volume get phold-vol out ./phold/out
Plain English: Phold turns each protein into a 3D-shape "alphabet" with the ProstT5 model (GPU helps here), then searches a
library of known protein shapes with Foldseek (CPU) to name what the protein does.
Cost guard: hard timeouts; each function prints elapsed seconds and approximate dollars."""
import os
import subprocess
import time

import modal

HERE = os.path.dirname(os.path.abspath(__file__))
GPU = os.environ.get("PHOLD_GPU", "T4")
PRICE = {"T4": 0.59, "L4": 0.80}[GPU]            # $/h, modal.com/pricing as read 25 Oct 2026 (T4); L4 approximate
CPU_PER_H = 8 * 0.047 + 32 * 0.008               # 8 cores + 32 GiB memory, approx $/h

app = modal.App("exp1-phold")
vol = modal.Volume.from_name("phold-vol", create_if_missing=True)
V = "/vol"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("curl", "ca-certificates")
    .pip_install("torch", "phold")
    .run_commands(
        "mkdir -p /opt/foldseek && curl -fsSL https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz | tar xz -C /opt/foldseek --strip-components=1",
        "/opt/foldseek/bin/foldseek version",
    )
    .env({"PATH": "/opt/foldseek/bin:/usr/local/bin:/usr/bin:/bin"})
    .add_local_file(os.path.join(HERE, "phold", "phold_in.fna"), "/app/phold_in.fna")
    .add_local_file(os.path.join(HERE, "phold", "phold_in2.fna"), "/app/phold_in2.fna")
    .add_local_file(os.path.join(HERE, "phold", "phold_in3.fna"), "/app/phold_in3.fna")
    .add_local_file(os.path.join(HERE, "phold", "phold_in4.fna"), "/app/phold_in4.fna")
    .add_local_file(os.path.join(HERE, "phold", "phold_in5.fna"), "/app/phold_in5.fna")
)


def _sh(cmd, timeout):
    print("$", cmd, flush=True)
    r = subprocess.run(cmd, shell=True, timeout=timeout)
    if r.returncode:
        raise RuntimeError(f"failed ({r.returncode}): {cmd}")


@app.function(image=image, volumes={V: vol}, cpu=4, memory=16384, timeout=5400)
def install_db():
    t0 = time.time()
    _sh(f"phold install -d {V}/db", 5000)
    vol.commit()
    _sh(f"du -sh {V}/db; ls {V}/db | head -30", 120)
    print(f"[cost] {time.time() - t0:.0f} s CPU container (about $0.05-0.20)")


@app.function(image=image, volumes={V: vol}, gpu=GPU, cpu=8, memory=32768, timeout=7200)
def run_phold(batch_size: int = 8, inp: str = "phold_in.fna", tag: str = "phold_run"):
    t0 = time.time()
    out = f"{V}/out/{tag}"
    _sh(f"rm -rf {out}; mkdir -p {V}/out && phold run -i /app/{inp} -o {out} -d {V}/db -t 8 --batch_size {batch_size} -f", 7000)
    vol.commit()
    _sh(f"ls -la {out} | head -30", 120)
    s = time.time() - t0
    print(f"[cost] {s:.0f} s on {GPU} = about ${s / 3600 * (PRICE + CPU_PER_H):.2f} (GPU + 8 cores + 32 GiB; check the Modal dashboard)")


@app.local_entrypoint()
def main():
    print(__doc__)
