"""step15_pharokka_metrics.py - proper M12 (and M5) from the real Pharokka run on Modal (pharokka/out_all, 1,804 contigs).
M12 (proper): per sequence, share of Pharokka's PHANOTATE genes that carry a PHROG assignment (mmseqs profile search OR pyhmmer HMM search,
              Pharokka's own combined 'phrog' column), and the share whose PHROG has a known function (category other than 'unknown function').
M5 (through Pharokka): gene count, coding density (0-1) and mean gene length from Pharokka's PHANOTATE calls.
Output: data/metrics_pharokka.csv (ids use the same naming as the other metric tables: CAL| for calibration genomes).
Also compares Pharokka's gene calls with the earlier raw-PHANOTATE run (data/metrics_phanotate.csv)."""
import numpy as np, pandas as pd
d = pd.read_csv("pharokka/out_all/pharokka_cds_final_merged_output.tsv", sep="\t", low_memory=False)
L = pd.read_csv("pharokka/out_all/pharokka_length_gc_cds_density.tsv", sep="\t")
d["has"] = d.phrog.astype(str) != "No_PHROG"
d["known"] = d.has & (d.category.astype(str) != "unknown function")
d["glen"] = (d.stop - d.start).abs() + 1
g = d.groupby("contig").agg(m12p_genes=("gene", "size"), m12p_hits=("has", "sum"), m12p_known=("known", "sum"), m5p_meanlen=("glen", "mean")).reset_index()
g["m12p_phrog_share"] = (g.m12p_hits / g.m12p_genes).round(4)
g["m12p_known_share"] = (g.m12p_known / g.m12p_genes).round(4)
g["m5p_genes"] = g.m12p_genes
g["m5p_meanlen"] = g.m5p_meanlen.round(1)
g = g.merge(L[["contig", "cds_coding_density"]], on="contig"); g["m5p_density"] = (g.cds_coding_density / 100).round(4)
g["id"] = g.contig.str.replace(r"^CAL_", "CAL|", regex=True)
out = g[["id", "m12p_genes", "m12p_hits", "m12p_phrog_share", "m12p_known_share", "m5p_genes", "m5p_density", "m5p_meanlen"]]
out.to_csv("data/metrics_pharokka.csv", index=False)
print("rows", len(out), "| genes per sequence median", out.m12p_genes.median())
# compare with the raw PHANOTATE run and the simplified M12
old5 = pd.read_csv("data/metrics_phanotate.csv"); old12 = pd.read_csv("data/metrics_phrog.csv")
m = out.merge(old5, on="id").merge(old12, on="id")
print("sequences compared:", len(m))
print("PHANOTATE gene count identical (Pharokka vs raw):", round((m.m5p_genes == m.m5_genes).mean(), 3),
      "| median abs difference in coding density:", round((m.m5p_density - m.m5_density).abs().median(), 4))
m["grp"] = m.id.map(lambda s: "CAL" if s.startswith("CAL|") else ("EVO2_" + s[-2:] if s.endswith(("__G1", "__G2")) else s.split("__")[1].split("_")[0]))
print("\nmedian share of genes with a PHROG: simplified M12 vs proper M12 (Pharokka)")
t = m.groupby("grp").agg(simplified=("m12_phrog_share", "median"), proper=("m12p_phrog_share", "median"),
                         simpl_known=("m12_known_share", "median"), proper_known=("m12p_known_share", "median")).round(2)
print(t.loc[["CAL", "natural", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "EVO2_G1", "EVO2_G2"]].to_string())
