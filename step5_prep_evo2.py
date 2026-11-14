"""step5_prep_evo2.py - prepare the small Evo 2 check (design section 6): NO GPU needed, costs nothing.
40 test genomes (seeded random pick) x 2 prompt lengths: G1 = first 1,000 letters, G2 = first 300 letters.
Evo 2 continues each prompt to the natural genome's own length. Writes evo2/prompts.fa and evo2/targets.tsv and prints the cost table."""
import random
def read_fa(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        else: d[n].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}
nat = read_fa("data/natural_test.fna")
pick = sorted(random.Random(1337).sample(sorted(nat), 40))
rows = []
with open("evo2/prompts.fa", "w") as fa, open("evo2/targets.tsv", "w") as tg:
    tg.write("id\tparent\tprompt_len\ttarget_len\tto_generate\n")
    for g in pick:
        for tag, pl in (("G1", 1000), ("G2", 300)):
            sid = f"{g}__{tag}"
            fa.write(f">{sid}\n{nat[g][:pl]}\n")
            tg.write(f"{sid}\t{g}\t{pl}\t{len(nat[g])}\t{len(nat[g]) - pl}\n"); rows.append(len(nat[g]) - pl)
tot = sum(rows)
print("sequences", len(rows), "letters to generate", tot, "(range per sequence", min(rows), "-", max(rows), ")")
print("\nspeed (letters/s)  hours   L40S $1.95/h   A100-80GB $2.50/h   H100 $3.95/h   [generation only, c=1]")
for sp in (60, 93, 109, 120):
    h = tot / sp / 3600
    print(f"{sp:>10}      {h:5.2f}   ${h*1.95:5.2f}        ${h*2.5:5.2f}           ${h*3.95:5.2f}")
print("\nplus setup/model load/lossless check: about 0.3-0.5 h more on the same GPU (assumed, untested)")
