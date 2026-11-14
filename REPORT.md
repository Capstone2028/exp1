# Experiment 1: which predicted-plausibility metrics can be trusted?

A damage-ladder calibration on Microviridae (phiX174-like) phage genomes, with a small Evo 2 check.
Design: `../exp1_design.docx` (29 Oct 2026). Day-by-day record, every decision and every failed attempt: `RUNLOG.md`.
How to rerun: `README.md`.

**No claim is made anywhere that a generated sequence is viable.** Every score here is predicted from DNA alone.

---

## 1. What was done, in plain English

Programs that judge "does this DNA look like a real genome?" are often used to screen DNA written by a model such as Evo 2.
Before trusting any of them, we need to know what each one actually notices. So we took real phage genomes, damaged them in
ten known ways (the "ladder"), and measured how well each metric tells a damaged genome from its natural original.

- **AUROC** is the score used throughout. For a metric it is the chance that a natural genome scores better than its damaged
  copy. 0.5 means the metric cannot tell them apart, 1.0 means it always can, and below 0.5 means the damaged copy scored *better*.
- A good metric catches the damage it should (AUROC at least 0.85 with the lower end of the 95% interval at least 0.75),
  stays near 0.5 on harmless changes (0.40 to 0.60), and passes at least 90% of natural genomes.
- The rules and the predicted catches were written down in the design before any data (`design.txt` is a text copy of its tables).
  The way I turned the loosely-worded "Drop for validity" rule into code, and the scoring direction of each metric, are written
  at the top of `step4_analyze.py` and were fixed before any AUROC was computed.

### The ten damages (D1 to D10)

| ID | Damage | Levels | Role |
|---|---|---|---|
| D1 | Synonymous recoding (same protein, different codons) | 10%, 30% of swappable codons | harmless control |
| D2 | Random point mutations | 0.5, 1, 2, 5, 10% | gradual damage |
| D3 | Frameshifts (insert/delete 1 letter inside genes) | 1, 3, 6, 12 | broken reading frames |
| D4 | Premature stops or removed start codons | 1, 3, 6 genes (see deviations) | truncated proteins |
| D5 | Whole genes deleted | 1, 3, 6 genes (see deviations) | missing parts |
| D6 | Segment shuffle (cut, reorder, flip) | 2, 4, 8 cuts | gene order broken, genes kept |
| D7 | Codon shuffle inside genes | 25, 50, 100% | proteins destroyed, gene-like signal kept |
| D8 | k-mer-only chain (6-letter word statistics, no neural network) | 1 | texture without structure |
| D9 | Wrong family (windows of unrelated phages and bacteria) | 1 | organism specificity |
| D10 | Reverse complement; rotation | 2 | harmless control |

## 2. Data

| Step | Count |
|---|---|
| Complete Microviridae genomes downloaded (NCBI) | 7,581 |
| 4,000 to 6,500 letters, only A/C/G/T | 7,194 |
| Clusters at about 95% identity (16-mer sketch, an approximation of BLAST) | 5,856 |
| Clusters that are strictly phiX174-like (6 or more of its 11 proteins) | 17 (design needs 90; its switch rule is under 60) |
| **Broader group used instead:** clusters carrying all five shared protein families F1 to F5 | 2,723 |
| Usable: all five as full genes, F2 first, at least 6 genes of 150+ letters | 229 |
| Split (seed 1337) | 30 calibration, 60 test |

Because only 17 clusters were phiX174-like, I followed the design's own fallback and used a broader Microviridae group with a
data-driven core-gene list ("option 1", chosen by the user). **M6, M8 and M9 therefore mean "the shared genes of this
Microviridae group", not "phiX174's 11 genes".**

The five core families were found by protein matching, not by name. Phold later annotated them on one genome: F2 is a
*replication initiation protein* (used as the anchor: every genome was flipped and rotated so F2 starts at position 1),
F1 a major head (capsid) protein, F3 a pilot protein for DNA ejection, F4 a DNA-binding protein, F5 a virion structural protein.

Each of the 60 test genomes gets 27 damaged versions: **1,620 damaged plus 60 natural = 1,680 sequences**
(`data/ladder.fna`, `data/ladder_manifest.tsv`).

## 3. The metric panel

| ID | Metric | Status |
|---|---|---|
| M1 | Basic gates (length, GC, longest run) | done |
| M2 | CheckV completeness (v1.1.1, db 1.5) | done |
| M3 | Closed-genome check (repeated ends) | done |
| M4 | Pyrodigal gene count, coding density, mean gene length | done |
| M5 | PHANOTATE (the gene finder inside Pharokka) | done; the real Pharokka gave the same gene calls (99.8% identical counts) |
| M6 | Core-gene checklist (of 5 families) | done |
| M7 | ORF integrity | done |
| M8 | Gene order against the reference | done |
| M9 | Start/stop layout (simplified, ignores where the genome starts) | done |
| M10 | Ribosome-binding-site motif upstream of genes | done |
| M11 | Pfam families (pyhmmer, gathering cutoffs) | done |
| M12 | PHROG families | done with the **real Pharokka** (run on Modal); an earlier simplified search is kept for comparison |
| M13 | Phold structure labels | done on **20 of the 60 genomes only** (220 sequences) |
| M14 | ESMFold confidence | not run |
| M15 | Foldseek fold match on our own proteins | not run |

## 4. Results: sensitivity grid

AUROC, natural versus damaged, at the middle severity (60 test genomes). The last column is the same test for rotation after the
post-hoc start-position fix (section 7).

| Metric | D1 30% | D2 2% | D3 6 | D4 3 | D5 3 | D6 4 | D7 50% | D8 | D9 | D10 rev | D10 rot | D10 rot, start-fixed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M1 Basic gates | 0.49 | 0.51 | 0.50 | 0.51 | 0.95 | 0.49 | 0.49 | 0.55 | 0.81 | 0.50 | 0.50 | 0.50 |
| M2 CheckV | 0.49 | 0.59 | 0.60 | 0.50 | 0.97 | 0.57 | 0.81 | 1.00 | 1.00 | 0.50 | 0.49 | 0.49 |
| M3 Closed-genome | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 |
| M4 Pyrodigal | 0.53 | 0.85 | 0.93 | 0.76 | 0.92 | 0.77 | 0.59 | 1.00 | 0.87 | 0.50 | 0.62 | 0.53 |
| M5 PHANOTATE | 0.50 | 0.80 | 0.89 | 0.68 | 0.79 | 0.71 | 0.53 | 1.00 | 0.78 | 0.48 | 0.60 | 0.51 |
| M6 Core genes | 0.50 | 0.56 | 0.58 | 0.56 | 0.98 | 0.61 | 0.92 | 1.00 | 1.00 | 0.50 | 0.50 | 0.50 |
| M7 ORF integrity | 0.54 | 0.32 | 0.21 | 0.45 | 0.76 | 0.65 | 0.57 | 0.67 | 0.85 | 0.50 | 0.67 | 0.51 |
| M8 Gene order | 0.50 | 0.67 | 0.71 | 0.59 | 0.75 | 0.71 | 0.82 | 1.00 | 1.00 | 0.50 | 0.60 | 0.50 |
| M9 Start/stop | 0.49 | 0.67 | 0.74 | 0.66 | 0.67 | 0.85 | 0.57 | 0.99 | 0.86 | 0.50 | 0.80 | 0.52 |
| M10 RBS signal | 0.62 | 0.70 | 0.75 | 0.71 | 0.64 | 0.65 | 0.71 | 0.94 | 0.68 | 0.50 | 0.51 | 0.50 |
| M11 Pfam | 0.43 | 0.60 | 0.65 | 0.54 | 0.53 | 0.55 | 1.00 | 1.00 | 0.05 | 0.50 | 0.51 | 0.49 |
| M12 PHROG (Pharokka) | 0.51 | 0.76 | 0.73 | 0.66 | 0.56 | 0.69 | 0.94 | 0.97 | 0.27 | 0.51 | 0.57 | not rerun |
| M12 PHROG (first, simplified) | 0.48 | 0.69 | 0.70 | 0.53 | 0.47 | 0.64 | 0.95 | 1.00 | 0.36 | 0.50 | 0.57 | 0.51 |
| M13 Phold (20 genomes) | 0.51 | 0.67 | 0.85 | 0.64 | 0.54 | 0.73 | 1.00 | 1.00 | 0.12 | n/a | 0.64 | 0.50 |

Intervals are in `results/auroc_mid.csv` (2,000 bootstrap resamples of genomes). M13's intervals are wide (20 genomes):
D3 0.85 (0.76 to 0.95), D4 0.64 (0.53 to 0.76). All severities: `results/auroc_all_severities.csv`; dose response (Spearman
correlation of severity with score): `results/dose_spearman.csv`; heat map: `results/grid.png`.

## 5. Verdicts

| Metric | Verdict | Why |
|---|---|---|
| M2 CheckV | **Keep** | Passes every check: D5 0.97, D8 1.00, D9 1.00; controls 0.49 to 0.50; 97% of natural genomes pass. Weak on D2/D3 (about 0.6), which were only partly predicted. |
| M1 Basic gates | **Keep, nearly empty** | Nothing was predicted for it. It only catches D5 (0.95) and D9 (0.81). |
| M9 Start/stop | **Keep, nearly empty** (reclassified by the user from Drop) | Its rotation false alarm disappears after the start fix (0.80 to 0.52). Only weak sensitivity remains (D4 predicted, 0.66). |
| M4 Pyrodigal | Keep with caveat | Misses D4 (0.76); only 82% of natural genomes pass; rotation false alarm 0.62 goes to 0.53 after the start fix. |
| M5 PHANOTATE | Keep with caveat | Misses D4 (0.68) and D5 (0.79); 87% natural pass. Its coding density barely drops on the k-mer chain. |
| M6 Core genes | Keep with caveat | Misses D3 (0.58) and D4 (0.56): protein matching tolerates a few broken genes; 87% natural pass. |
| M8 Gene order | Keep with caveat | Misses D5 (0.75) and D6 (0.71). Rotation 0.60 goes to 0.50 after the start fix. It reads circular neighbour pairs, so the fix cannot undo a shuffle for it (D6 0.71 to 0.72). |
| M10 RBS signal | Keep with caveat | Reacts to D1 (0.62). This is a leak in the control, not a bug (section 6). 88% natural pass. |
| M11 Pfam | Keep with caveat | Misses D3 (0.65) and D4 (0.54). Blind to D9 as predicted (0.05). |
| M12 PHROG (Pharokka) | Keep with caveat | Misses D3 (0.73) and D4 (0.66); weak on D5 (0.56). Only **70%** of natural genomes pass (below the 90% rule; the first simplified search gave 88%). Blind to wrong family (0.27). The start fix was not rerun for it. |
| M13 Phold | Keep with caveat (indicative) | Misses D4 (0.62 after the start fix). Natural pass rate 90% (exactly at the limit). Only 20 genomes. |
| M7 ORF integrity | **Drop** | Damaged genomes score *higher* than natural (D3 0.21, D2 0.32): the gene finder splits a broken gene into more whole-looking genes. Rotation does not explain it. |
| M3 Closed-genome | **Drop** | Never reacts to anything (records are circular with no repeated ends, as the design guessed). |
| M14, M15 | not run | |

Natural pass rates (share of the 60 natural test genomes inside the calibration set's 5th to 95th percentile): M1 97%, M2 97%,
M3 100%, M4 82%, M5 87%, M6 87%, M7 92%, M8 92%, M9 90%, M10 88%, M11 92%, M12 88%, M13 90%.

### Findings that matter

1. **Gene finders (M4, M5, M7) are blind to codon shuffling (D7)**, exactly as predicted (0.59, 0.53, 0.57): the reading frame and
   codon usage survive. **Protein-level metrics see it:** M6 0.92, Pfam 1.00, PHROG 0.95, Phold 1.00.
2. **Generic family metrics are blind to the wrong family (D9).** Pfam 0.05, PHROG 0.36 and Phold 0.12: foreign phage proteins are
   labelled *more* readily than the natural ones. Only the family-specific metrics catch it (M2 1.00, M6 1.00, M8 1.00). This is
   why Phold is a generous labeller, and why a low label rate is stronger evidence than a high one.
3. **The k-mer chain (D8) is caught by almost everything that reads genes or proteins** (0.94 to 1.00 for M2, M4, M5, M6, M8, M9, M10, M11, M12, M13) but not by the basic gates (M1 0.55) or ORF integrity (M7 0.67). M4 and M5 catch it through gene statistics, not composition: PHANOTATE still finds 76% coding density in k-mer noise where pyrodigal finds 25%.
4. **Agreement.** M4 and M5 correlate 0.74 over all ladder sequences but only 0.24 across natural genomes alone; Pfam and PHROG
   0.46 overall and 0.14 on natural genomes. They agree on whether something is damaged, not on which natural genome is better.
5. **Frameshifts and stops are the weak spot.** D3 and D4 are missed by most protein metrics because these genomes have only 6 to
   9 callable genes and a few broken ones leave most of the signal intact.

## 6. Why does M10 react to D1?

M10 reads 16 letters upstream of each gene for a GGAG/GAGG motif, so any reaction comes from changes inside those windows by construction.
54% of those window letters lie inside the upstream gene's coding sequence, so a synonymous swap can change them. But D1 does **not**
cluster there: 1.3% of the changed letters fall in windows, about the background rate (1.07 times the coding-letter baseline).
Natural windows are rich in the motif, so incidental hits lose it more often than they create it: D1 at 30% removes it from 16% of
genes that have it and adds it to 5% of genes that lack it (matched random point mutations remove it from 31%).
Conclusion: D1 is not a clean control for DNA-motif metrics. Whether it harms translation initiation cannot be told from these data.
(`step10_m10_d1.py`)

## 7. Post-hoc check: standardising where the circle starts

M4 (0.62) and M8 (0.60) failed the D10 rotation control because a linear gene finder cuts genes at the ends of the sequence.
**Not part of the frozen plan; report it next to the original.** Each sequence is rotated so the start codon of its single intact
F2 gene is at position 1 (reverse-strand genes are anchored at the other end so reverse complements are untouched). Sequences with
no recognisable F2 (the k-mer chain, wrong-family windows, heavily broken genes) are left alone. 224 of 1,680 sequences moved;
no natural genome moved. Result: rotation false alarms fall to 0.50 to 0.53 for M4, M5, M7, M8, M9, M12 and M13. Cost: D6
detection by M7 falls 0.65 to 0.50 and by M9 0.85 to 0.73. It works only where F2 is recognisable, and rotation back to the original
start gives Phold an identical sequence (so its 0.50 is by construction). Results: `std_analysis/results/`. (`step7_standardise.py`)

## 8. The Evo 2 check

Evo 2 7B (no fine-tuning), run through the EvSpark speed-up (temperature 1.0, top-k 4, one candidate at a time), continued the first
1,000 letters (G1) or 300 letters (G2) of 40 test genomes to each genome's own length: **80 outputs, 353,504 letters generated**.
EvSpark's own check passed (256 of 256 tokens identical to plain Evo 2 when run greedily; 116 against 42 tokens/s).
Four outputs each contain one NUL character that the model produced; they were kept as generated and replaced by `N` only so the
tools could read them (`evo2/outputs_raw.fna` is untouched).

### 8.1 Against the natural range (share of outputs inside the calibration 5th to 95th percentile)

| Metric | Natural genomes | G1 | G2 |
|---|---|---|---|
| M1 Basic gates | 97% | 100% | 95% |
| M2 CheckV completeness (median) | 100% | 82% (55% inside) | 17% (25% inside) |
| M4 Pyrodigal | 82% | 42% | 40% |
| M5 PHANOTATE | 87% | 48% | 32% |
| M6 Core genes (median of 5) | 5 (87% inside) | 2 (0%) | 1 (0%) |
| M7 ORF integrity | 92% | 85% | 65% |
| M8 Gene order | 92% | 32% | 8% |
| M10 RBS signal | 88% | 72% | 60% |
| M11 Pfam | 92% | 20% | 10% |
| M12 PHROG (Pharokka) | 70% | 15% | 5% |

The "equivalent damage" (the damage level whose average metric profile is closest to the output) is **codon shuffle (D7)
for 68 of 80 outputs** (33 of 40 G1, 35 of 40 G2; 66 of 80 with the first, simplified M12); the rest scatter (D2 3, D5 3, D8 3, D6 1, D4 1, natural 1). In words: gene-like structure and
natural-looking composition survive, while protein families and gene order are largely lost. This is a crude summary, not a
statement that Evo 2 "applies" that damage.

### 8.2 Phold on the 80 outputs (and the 60 natural test and 30 calibration genomes)

Share of genes labelled by Phold, per sequence (median): natural 0.65, calibration 0.62, **Evo 2 G1 0.25, G2 0.16**
(high-confidence only: 0.47, 0.47, 0.12, 0.04). Evo 2 outputs fall inside the natural range 12% (G1) and 10% (G2) of the time,
against 90% for natural test genomes. Because Phold labels foreign proteins readily (D9 = 0.12), a rate this far below natural is
strong evidence the proteins are broken, not merely novel. The roughly 19% of generated genes that are labelled (mostly head and
packaging, and replication) show that something phage-like is present, not that the right genes are there.

Compared protein by protein with the simplified PHROG search (M12) on the same proteins: in natural genomes both label 32% of
genes, Phold alone 30%, M12 alone 7%. What Phold adds is mostly head and packaging proteins (86 of 132) and replication-type
proteins (39), the divergent proteins that structure recognises and sequence search misses. Part of this gap is the simplified M12,
not structure-versus-sequence.

This comparison used the first, simplified M12 on Phold's own gene calls. The real Pharokka M12 uses a different (PHANOTATE) gene set, so it was compared only sequence by sequence: it labels a median 0.33 of natural genes (calibration range 0.15 to 0.56), against 0.08 for Evo 2 G1 and 0.00 for G2; inside the natural range 15% (G1) and 5% (G2).

### 8.3 Correction: part of every G1 result is the copied prompt

The G1 prompt (1,000 letters) is verbatim natural sequence. In 28 of 41 cases F2 ends before letter 1,000 (median F2 length 978),
so a detection can be Phold or a metric reading back the prompt. Recomputed counting only genes that reach generated letters
("overlap") or start after the prompt ("wholly generated"). Median per output, G1 / G2, and the mean paired G1 minus G2 difference
with its 95% interval over 40 genomes:

| Metric | Version | G1 / G2 | G1 minus G2 |
|---|---|---|---|
| Core genes (of 5) | original | 2 / 1 | +0.75 [+0.40, +1.12] |
| | overlap | 1 / 1 | +0.17 [-0.20, +0.55] |
| | wholly generated | 1 / 0 | +0.45 [+0.10, +0.80] |
| Pfam share | original | 0.25 / 0.11 | +0.14 [+0.07, +0.21] |
| | overlap | 0.18 / 0.11 | +0.07 [0.00, +0.14] |
| | wholly generated | 0.12 / 0.00 | +0.08 [+0.01, +0.15] |
| PHROG share (simplified) | original | 0.13 / 0.12 | +0.06 [+0.01, +0.10] |
| | overlap | 0.10 / 0.12 | -0.03 [-0.07, +0.02] |
| | wholly generated | 0.00 / 0.00 | +0.04 [0.00, +0.08] |
| Phold share | original | 0.25 / 0.16 | +0.09 [+0.02, +0.16] |
| | overlap | 0.18 / 0.16 | +0.01 [-0.06, +0.09] |
| | wholly generated | 0.13 / 0.00 | +0.04 [-0.03, +0.11] |

With the real Pharokka M12 the picture is the same. G1 minus G2 (paired, 95% interval): all genes +0.04 [0.00, +0.08]; genes reaching generated letters +0.01 [-0.02, +0.05] (median 0.00 / 0.00); wholly generated +0.01 [-0.03, +0.05]. No prompt-length effect once copied sequence is excluded.

For genes that reach generated letters, Phold labels about 0.19 (G1) and 0.18 (G2), against 0.63 for natural genes; 32 of the 34
genes lying wholly inside a G1 prompt are labelled (97%, natural copies). Restricting to outputs with a findable F2 (38 of 80 have
exactly one intact F2; 30 reach generated letters) does not change the headline: core genes stay at a median of 1 and closest
damage stays D7. The earlier statement that the shorter prompt does worse on every protein-level metric was largely "the longer
prompt hands over more intact sequence". A small residual G1 advantage on Pfam and core genes under the strict rule cannot be
separated from a region difference (G1's generated region starts at letter 1,000, G2's at 300) or from real conditioning.

### 8.4 CheckV controls (does the 82% versus 17% gap survive?)

CheckV scores whole sequences by comparison with references, so genes cannot be excluded; three controls instead. Median completeness:

| | G1 | G2 |
|---|---|---|
| Original output | 82.2 | 17.4 |
| Prompt only (also the length-matched natural start) | 20.3 | 6.1 |
| Generated letters only | 64.4 | 14.4 |
| Natural genome, same region and length as the generated part | 81.4 | 95.0 |

(11 of 40 G1 and 15 of 40 G2 generated parts were "not determined" and count as 0.) G1 minus G2 gap, percentage points,
mean of 40 paired genomes [95% interval]: original +29.8 [+13.7, +44.5]; prompt floor +14.5 [+13.3, +15.6]; original minus prompt-only
+15.3 [-0.2, +30.4]; generated only +20.6 [+4.8, +35.5]; length-matched natural -12.1 [-13.4, -10.4].
So about half the gap is the prompt's own length, the rest is borderline after subtracting it, and on generated letters alone G1
is still above G2 (and more so once length is allowed for). **This does not support "no prompt-length effect" for CheckV, and
it disagrees with Phold.** I cannot say which is right: CheckV rewards similarity to known genomes, the generated-only scores are
bimodal, and about 30% of G1 prompts end inside F2 so the output starts by finishing a gene begun from natural context.

### 8.5 What the Evo 2 check supports

- The outputs are phage-like in composition and gene structure, and far below natural genomes on every protein-level and
  family-level measure (core genes 1 to 2 of 5, Phold 0.16 to 0.25 against 0.65, Pfam 0.11 to 0.25 against 0.57).
- A longer prompt gave little or no clear advantage once copied sequence is excluded, except possibly on Pfam, core genes and CheckV.
- Nothing here says the outputs are viable or non-viable.

## 9. Deviations and limitations (read before citing any number)

1. **Family changed.** Only 17 clusters were strictly phiX174-like; a broader Microviridae group was used (design fallback). M6/M8/M9 are about this group's five shared families.
2. **Few genes.** These genomes have only 6 to 9 callable genes (the design assumed about 11). D4 and D5 "6 genes" would wreck the whole genome, so they were capped to leave at least 2 genes (real counts in `ladder_manifest.tsv`). D3/D4 therefore hit few genes.
3. **Only 229 of 2,723 candidate genomes were usable**, and the 60 test genomes are drawn from them; clustering is a 16-mer sketch, not BLAST.
4. **M12 and Pharokka.** M12 was first a simplified search (the profile search ran out of memory at 7 GB locally) and M5 was run as raw PHANOTATE. Both were later redone with the real Pharokka 1.10.1 on a Modal cloud machine (run with `--skip_extra_annotations --skip_mash`, since the tRNA, CRISPR and genome-match steps do not feed M5 or M12). M5 is unchanged (99.8% identical gene counts). The proper M12 labels slightly fewer genes (more, shorter genes are called) and its natural pass rate falls to 70%. The start-position fix and the Phold protein-by-protein comparison were not rerun with it.
5. **M13 (Phold) is on 20 genomes only**; M14 and M15 were not run.
6. **Post-hoc and user-requested changes:** the start-position standardisation (section 7), the reclassification of M9, adding M2 after CheckV finished, and the prompt-exclusion and CheckV controls (section 8.3 and 8.4).
7. **D1 is not neutral for DNA-motif metrics** (section 6).
8. **The Evo 2 check is one run per prompt, unfine-tuned, on 40 genomes**; 4 of 80 outputs contain a NUL character.
9. M8's reference order comes from the calibration set (calibration genomes score 1.0 by construction); the natural test genomes are the honest baseline.

## 10. Cost

About **$2.3 to $2.6** of the $30 Modal credit (summed from the per-job cost lines each Modal run printed; the Modal dashboard is the authority): Evo 2 setup, checks and trial about $0.25, the full 80-sequence run $1.45 (2,673 s on an L40S), and Phold about $0.15 to $0.30 in total (database setup, the 220-sequence run, the 170-sequence run, the rotated re-run). Pharokka on Modal (database, two timing tests and the full run) about $0.25. Everything else ran on the local CPU (WSL).

## 11. Where the numbers are

| What | File |
|---|---|
| Ladder sequences and what each damage did | `data/ladder.fna`, `data/ladder_manifest.tsv` |
| Split and gene coordinates | `data/split.tsv`, `data/genes_oriented.tsv` |
| Metric values per sequence | `data/metrics_quick.csv` (M1, M3, M4, M6 to M10), `metrics_phanotate.csv` (M5), `metrics_pfam.csv` (M11), `metrics_phrog.csv` (M12), `metrics_checkv.csv` (M2), `metrics_phold.csv` (M13) |
| AUROC grid, dose response, verdicts | `results/` |
| Start-position fix | `std/`, `std_analysis/` |
| Evo 2 prompts, outputs, F2 findability, prompt-excluded numbers | `evo2/` |
| Phold per-gene predictions | `phold/out/`, `phold/out2/`, `phold/out3/` |
