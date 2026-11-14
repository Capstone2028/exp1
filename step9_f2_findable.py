"""step9_f2_findable.py - for the 80 Evo 2 outputs: is a family-F2 gene findable (same rule as step7_standardise.py)?
Run in WSL: ~/kmer-venv/bin/python step9_f2_findable.py      Output: evo2/f2_findable.tsv
  any_F2    : at least one gene labelled F2 (6-letter protein word vote, same as step2/step7)
  intact_F2 : exactly one F2 gene whose length is 0.8-1.25x the calibration median (948 letters) = the condition the rotation fix needs
  core_found: how many of F1-F5 have any gene labelled with that family"""
import sys, statistics
import os
sys.argv = ["x", os.getcwd(), "/tmp/zz9"]          # run from the exp1 folder
src = open("step7_standardise.py").read().split("lad = read_fa")[0]
exec(src)
ev = read_fa("evo2/outputs.fna")
PL = {l.split("	")[0]: int(l.split("	")[2]) for l in list(open("evo2/targets.tsv"))[1:]}
rows = []
for sid, s in ev.items():
    L = len(s); fams = []
    for g in finder.find_genes((s + s + s).encode()):
        a, b = g.begin - 1 - L, g.end - L
        if not (0 <= a < L) or g.partial_begin or g.partial_end or b - a > L:
            continue
        f = family_of(g.translate().rstrip("*"))
        if f: fams.append((f, b - a, a, b))
    f2 = [x for x in fams if x[0] == "F2"]
    pl = PL[sid]
    intact = len(f2) == 1 and 0.8 * F2LEN <= f2[0][1] <= 1.25 * F2LEN
    in_prompt = int(any(x[3] <= pl for x in f2))               # an F2 gene lying wholly inside the copied prompt
    gen_part = int(any(x[3] > pl for x in f2))                 # an F2 gene reaching into generated letters
    rows.append((sid, len(f2), int(intact), len({x[0] for x in fams}), in_prompt, gen_part, int(intact and f2[0][3] > pl)))
with open("evo2/f2_findable.tsv", "w") as fh:
    fh.write("id\tn_F2\tintact_F2\tcore_found\tF2_inside_prompt\tF2_reaches_generated\tintact_F2_reaching_generated\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
for tag in ("G1", "G2"):
    rr = [r for r in rows if r[0].endswith(tag)]
    print(tag, "n", len(rr), "any F2:", sum(r[1] > 0 for r in rr), "intact single F2:", sum(r[2] for r in rr), "| F2 wholly inside prompt:", sum(r[4] for r in rr), "| F2 reaching generated letters:", sum(r[5] for r in rr), "| intact F2 reaching generated:", sum(r[6] for r in rr), "median core families:", statistics.median(r[3] for r in rr))
