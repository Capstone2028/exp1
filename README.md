# exp1 - which predicted-plausibility metrics can be trusted?

Damage-ladder calibration on Microviridae genomes plus a small Evo 2 check.
Results and interpretation: **`REPORT.md`**. Chronological record of every step, decision and failed attempt: **`RUNLOG.md`**.
Design document: `../exp1_design.docx` (`design.txt` is a plain-text extract of it).

## What is in the repository, and what is not

Tracked: all scripts, the report and log, the small data the analysis needs (the ladder and its manifest, the split, gene coordinates,
metric tables, results, Evo 2 prompts/outputs, Phold per-gene tables, CheckV summaries). About 25 MB.

Not tracked (see `.gitignore`; all rebuildable or re-downloadable):

| Item | Size | How to get it back |
|---|---|---|
| `db/pfam/` Pfam-A.hmm | 2.2 GB | `ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`, unpack |
| `db/checkv/` CheckV db v1.5 | about 6 GB | `checkv download_database ./db/checkv` |
| `db/pharokka/` PHROG database (Pharokka DB v1.2.0, Zenodo 7563578) | 0.3 GB | download the tarball from Zenodo, unpack |
| `db/phold/` Phold database + ProstT5 | 15 GB | `phold install -d db/phold` (the Modal run keeps its own copy in a Volume) |
| `data/mv/` NCBI Microviridae download | 13 MB zip | `datasets download virus genome taxon Microviridae --complete-only` |
| `data/blastp_all.tsv` all-vs-all protein search | 1.4 GB | rerun `step1b_families.py` (about 1 h on 8 cores) |
| `data/cluster_reps.fna`, scratch and CheckV `tmp/` folders, raw logs | - | produced by the scripts |

## Environments used

- **Windows (Python 3.13, pandas, numpy, matplotlib):** `step1_cluster.py`, `step4_analyze.py`, `step5_prep_evo2.py`, `step6`, `step8`, `step8b`, `step9b`, `step11`, `step12b`, `step13`, `gen_exp1.py` (CPU test only).
- **WSL Ubuntu (a venv with pyrodigal, pyhmmer, phanotate, checkv; mmseqs2, blast+, diamond, hmmer from apt):** `step1b` to `step3d`, `step2_build.py`, `step7`, `step9_f2_findable.py`, `step10`, `step12_prompt_excluded.py`, `test_ladder.py`. PHANOTATE needs `setuptools<81`.
- **Phold (WSL, separate venv with CPU torch and the Foldseek binary on PATH)** was only used for tests; the real runs were on Modal.
- **Modal (T4 for Phold, L40S for Evo 2):** `modal_phold.py` here, and `../kmer-evo2/modal_app.py` (function `exp1_gen`) for Evo 2 + EvSpark. Both are untested outside this experiment's runs.

Windows Application Control blocked pyrodigal's compiled module on Windows, which is why gene calling runs in WSL.

## Order of the scripts

1. `step1_cluster.py` filter and cluster; `step1b_families.py` (WSL) call genes and group proteins into families; `step1c_group.py` stricter family rule and summary.
2. `step2_build.py` (WSL) orient, split 30/60, build the 27-version damage ladder; `test_ladder.py` checks each damage did what it should.
3. Metrics: `step3_quick_metrics.py` (M1, M3, M4, M6 to M10), `step3c_phanotate.py` (M5), `step3b_pfam.py` (M11), `step3d_phrog.py` (M12), CheckV (M2, commands in `RUNLOG.md`), `modal_phold.py` + `step8_phold_analysis.py` (M13).
4. `step4_analyze.py` AUROCs, dose response, natural pass rate, agreement, verdicts, heat map (writes `results/`).
5. Evo 2: `step5_prep_evo2.py`, `gen_exp1.py` via Modal, `step6_evo2_analysis.py`.
6. Follow-ups: `step7_standardise.py` (start-position fix, then rescore the moved sequences and rerun `step4_analyze.py` in `std_analysis/`), `step8b_phold_rotated.py`, `step9_f2_findable.py`, `step9b_evo2_restricted.py`, `step10_m10_d1.py`, `step11_phold_main.py`, `step12_prompt_excluded.py` + `step12b_compare.py`, `step13_checkv_controls.py`.

Everything is seeded (1337); rerunning gives the same ladder and bootstrap intervals.
