"""test_ladder.py - check that each damage type did what the design says. Run in WSL after step2_build.py."""
import collections, sys
sys.argv = ["x"]
exec(open("step2_build.py").read().split("# ---- damage functions")[0].split("oriented, dropped")[0])   # helpers only
nat = read_fa("data/natural_test.fna")
lad = read_fa("data/ladder.fna")
man = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in list(open("data/ladder_manifest.tsv"))[1:]}
genes = collections.defaultdict(list)
for l in list(open("data/genes_oriented.tsv"))[1:]:
    g, st, i, a, b, s, f, pl = l.rstrip("\n").split("\t")
    if st == "test":
        genes[g].append((int(a), int(b), int(s)))
tr = lambda s: "".join(CODON[s[i:i + 3]] for i in range(0, len(s) - 2, 3))


def prots(seq, gl):
    return [tr(seq[a:b] if st == 1 else rc(seq[a:b])) for a, b, st in gl]


def solo(gl):
    cov = collections.Counter(i for a, b, _ in gl for i in range(a, b))
    return [k for k, (a, b, _) in enumerate(gl) if all(cov[i] == 1 for i in range(a, b))]


bad, c = [], collections.Counter()
for sid, seq in lad.items():
    _, parent, did, sev, L, note = man[sid]
    ref, gl = nat[parent], genes[parent]
    P0 = prots(ref, gl)
    ok = set(seq) <= set("ACGT")
    if did == "natural":
        ok &= seq == ref
    elif did == "D1":
        P = prots(seq, gl); ok &= len(seq) == len(ref) and seq != ref and all(P[k] == P0[k] for k in solo(gl))
    elif did == "D2":
        ok &= len(seq) == len(ref) and sum(x != y for x, y in zip(seq, ref)) == round(float(sev) / 100 * len(ref))
    elif did == "D3":
        ok &= abs(len(seq) - len(ref)) <= int(sev) and seq != ref
    elif did == "D4":
        ok &= len(seq) == len(ref) and seq != ref and sum(x != y for x, y in zip(seq, ref)) <= 3 * int(sev)
    elif did == "D5":
        ok &= len(seq) < len(ref)
    elif did == "D6":
        ok &= len(seq) == len(ref) and sorted(seq) == sorted(ref) or sorted(seq) == sorted(rc(ref)) or len(seq) == len(ref)
        ok &= seq not in ref + ref and seq not in rc(ref) + rc(ref)
    elif did == "D7":
        P = prots(seq, gl)
        ok &= len(seq) == len(ref) and all(sorted(P[k][1:-1]) == sorted(P0[k][1:-1]) and P[k][0] == P0[k][0] for k in solo(gl))
        ok &= any(P[k] != P0[k] for k in solo(gl))
    elif did == "D8":
        ok &= len(seq) == len(ref) and seq != ref
    elif did == "D9":
        ok &= 5000 <= len(seq) <= 6000
    elif did == "D10":
        ok &= (seq == rc(ref)) if sev == "revcomp" else (len(seq) == len(ref) and seq in ref + ref and seq != ref)
    c[(did, bool(ok))] += 1
    if not ok: bad.append(sid)
print(sorted(c.items()))
print("FAILED:", bad[:10])
# D7 100% on the first test genome: how different are the proteins from the natural ones?
g = sorted(nat)[0]
print("D7 100% protein identity (should be low) :", [round(sum(x == y for x, y in zip(a, b)) / len(a), 2) for a, b in zip(prots(nat[g], genes[g]), prots(lad[g + "__D7_100"], genes[g]))])
print("D1 30% nucleotide identity to natural:", round(sum(x == y for x, y in zip(nat[g], lad[g + "__D1_30"])) / len(nat[g]), 3))
