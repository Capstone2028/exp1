"""step10_m10_d1.py - why does M10 (ribosome-binding-site motif GGAG/GAGG, 4-20 letters before each gene start) react to D1?
Run in WSL: ~/kmer-venv/bin/python step10_m10_d1.py
For the 60 natural test genomes and their D1 versions (10% and 30% of swappable codons), using the natural genome's own gene calls
(pyrodigal meta, same as M4/M10), whose coordinates D1 does not move:
  A. How much of the RBS window (16 letters/gene) sits where D1 is ALLOWED to change letters (single-coverage coding letters)?
  B. Of the letters D1 actually changed, what share lie inside RBS windows, versus the share of all coding single-coverage letters
     that are in windows (enrichment 1.0 = D1 does not target windows, it just hits them at the background rate)?
  C. Per gene: motif present in natural and lost in D1 / absent in natural and gained in D1.
  D. Same C for D2 point mutations at the same number of changed letters (is D1 special, or does any change do this?).
Reading: M10 only looks at the window, so every M10 reaction comes from changes inside windows BY CONSTRUCTION. What can differ is
whether the windows are hit more than chance (B), and whether motifs are lost more than gained (C): natural windows are enriched
for the motif, random changes destroy more than they create."""
import collections, random, re
import numpy as np, pyrodigal

COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
finder = pyrodigal.GeneFinder(meta=True)


def read_fa(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        else:
            d[n].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


lad = read_fa("data/ladder.fna")
tests = sorted({k.split("__")[0] for k in lad})
MOT = re.compile("GGAG|GAGG")


def window(seq, a, e, st):                   # exactly the M10 definition
    if st == 1:
        return seq[max(0, a - 20):max(0, a - 4)], list(range(max(0, a - 20), max(0, a - 4)))
    return rc(seq[e + 4:e + 20]), list(range(e + 4, min(len(seq), e + 20)))


res = {lv: collections.Counter() for lv in ("D1_10", "D1_30", "D2_matched")}
rng = random.Random(1337)
for g in tests:
    nat = lad[g + "__natural"]
    genes = [(x.begin - 1, x.end, x.strand) for x in finder.find_genes(nat.encode())]
    cov = np.zeros(len(nat), int)
    for a, e, st in genes: cov[a:e] += 1
    single = cov == 1
    win_pos = set()
    for a, e, st in genes: win_pos.update(window(nat, a, e, st)[1])
    for lv in ("D1_10", "D1_30"):
        alt = lad[f"{g}__{lv}"]
        changed = [i for i in range(len(nat)) if nat[i] != alt[i]]
        c = res[lv]
        c["n_changed"] += len(changed)
        c["changed_in_window"] += sum(i in win_pos for i in changed)
        c["window_letters"] += len(win_pos)
        c["window_letters_single_coding"] += sum(single[i] for i in win_pos)
        c["single_coding_letters"] += int(single.sum())
        c["genome_letters"] += len(nat)
        for a, e, st in genes:
            w0, _ = window(nat, a, e, st); w1, _ = window(alt, a, e, st)
            m0, m1 = bool(MOT.search(w0)), bool(MOT.search(w1))
            c["genes"] += 1; c["motif_nat"] += m0; c["lost"] += (m0 and not m1); c["gained"] += (m1 and not m0)
        if lv == "D1_30":                    # D2 with the same number of changed letters, anywhere
            arr = list(nat)
            for i in rng.sample(range(len(nat)), len(changed)):
                arr[i] = rng.choice([x for x in "ACGT" if x != arr[i]])
            d2 = "".join(arr); c2 = res["D2_matched"]
            for a, e, st in genes:
                w0, _ = window(nat, a, e, st); w1, _ = window(d2, a, e, st)
                m0, m1 = bool(MOT.search(w0)), bool(MOT.search(w1))
                c2["genes"] += 1; c2["motif_nat"] += m0; c2["lost"] += (m0 and not m1); c2["gained"] += (m1 and not m0)
            c2["n_changed"] += len(changed)
for lv, c in res.items():
    print(f"\n== {lv}  ({len(tests)} genomes, {c['genes']} genes)")
    if "changed_in_window" in c:
        wfrac = c["window_letters_single_coding"] / c["window_letters"]
        print(f"A. window letters that D1 is allowed to change (single-coverage coding): {wfrac:.0%} of {c['window_letters']} window letters")
        print(f"B. changed letters inside RBS windows: {c['changed_in_window']} of {c['n_changed']} = {c['changed_in_window'] / c['n_changed']:.1%}; "
              f"windows are {c['window_letters'] / c['genome_letters']:.1%} of all letters and {c['window_letters_single_coding'] / c['single_coding_letters']:.1%} of single-coverage coding letters")
        print(f"   enrichment vs all letters {c['changed_in_window'] / c['n_changed'] / (c['window_letters'] / c['genome_letters']):.2f}x; "
              f"vs coding letters share of changes in windows (changes only occur in coding): "
              f"{(c['changed_in_window'] / c['n_changed']) / (c['window_letters_single_coding'] / c['single_coding_letters']):.2f}x")
    print(f"C. genes with the motif in natural: {c['motif_nat']} ({c['motif_nat'] / c['genes']:.0%}); lost {c['lost']} ({c['lost'] / max(c['motif_nat'], 1):.0%} of those); "
          f"gained {c['gained']} ({c['gained'] / (c['genes'] - c['motif_nat']):.0%} of genes without it)")
