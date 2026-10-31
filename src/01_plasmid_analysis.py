# BEN 2026-11-01
# Change aim: Analyse GC content, ORFs, codon usage and proteins of plasmid pPCP1.
# File: 01_plasmid_analysis.py
"""Genome-level analysis of the Yersinia pestis plasmid pPCP1 (NCBI NC_005816, 9,609 bp).
GC content, GC skew, ORF scan vs annotated genes, codon usage, protein properties."""
import os
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

os.makedirs("results", exist_ok=True)
rec = SeqIO.read("data/NC_005816.gb", "genbank")
seq = rec.seq
out = [f"{rec.id}  {rec.description}", f"length={len(seq)} bp  GC={gc_fraction(seq)*100:.2f}%"]

# --- sliding-window GC and GC skew (G-C)/(G+C)
W, S = 300, 50
pos, gc, skew = [], [], []
for i in range(0, len(seq) - W, S):
    w = str(seq[i:i + W]); g, c = w.count("G"), w.count("C")
    pos.append(i + W // 2); gc.append((g + c) / W * 100); skew.append((g - c) / (g + c))
fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
ax[0].plot(pos, gc); ax[0].set_ylabel("GC %")
ax[1].plot(pos, np.cumsum(skew) * S / W); ax[1].set_ylabel("cumulative GC skew"); ax[1].set_xlabel("position (bp)")
cds_feats = [f for f in rec.features if f.type == "CDS"]
for f in cds_feats:
    for a in ax:
        a.axvspan(int(f.location.start), int(f.location.end), color="grey", alpha=0.12)
ax[0].set_title("pPCP1: GC content and GC skew (grey = annotated CDS)")
plt.tight_layout(); plt.savefig("results/gc_skew.png", dpi=150); plt.close()

# --- ORF scan (6 frames, >=100 aa) vs annotation
# BEN 2026-11-01
# Change aim: Scan six reading frames for ORFs of at least minlen amino acids.
# File: 01_plasmid_analysis.py
# Function: orfs
def orfs(s, minlen=100):
    found = set()
    for strand, ss in ((1, s), (-1, s.reverse_complement())):
        for fr in range(3):
            prot = str(ss[fr:len(ss) - ((len(ss) - fr) % 3)].translate())
            p = 0
            for chunk in prot.split("*"):
                m = chunk.find("M")
                if m != -1 and len(chunk) - m >= minlen:
                    st = fr + 3 * (p + m); en = fr + 3 * (p + len(chunk)) + 3
                    a, b = (st, en) if strand == 1 else (len(s) - en, len(s) - st)
                    found.add((a, b, strand))
                p += len(chunk) + 1
    return sorted(found)
found = orfs(seq)
ann = {(int(f.location.start), int(f.location.end), f.location.strand) for f in cds_feats}
# stop codon position: end coordinate on the + strand, start coordinate on the - strand
def stop_key(a, b, strand):
    return (b if strand == 1 else a, strand)
ann_stops = {stop_key(a, b, st) for a, b, st in ann}
hits = sum(1 for a, b, st in found if stop_key(a, b, st) in ann_stops)
ann100 = sum(1 for f in cds_feats if len(f.qualifiers["translation"][0]) >= 100)
out.append(f"annotated CDS={len(ann)} (>=100 aa: {ann100})  ORFs>=100aa (ATG start) found={len(found)}  sharing a stop codon with an annotated CDS={hits}")

# --- per-CDS table: protein length, MW, pI, GC
rows = []
for f in cds_feats:
    prot = f.qualifiers["translation"][0]
    pa = ProteinAnalysis(prot)
    nt = f.extract(seq)
    rows.append((f.qualifiers.get("gene", ["-"])[0], f.qualifiers.get("product", ["-"])[0],
                 len(prot), pa.molecular_weight() / 1000, pa.isoelectric_point(), gc_fraction(nt) * 100))
with open("results/cds_table.tsv", "w") as fh:
    fh.write("gene\tproduct\taa_len\tMW_kDa\tpI\tGC_pct\n")
    for r in rows:
        fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.1f}\t{r[4]:.2f}\t{r[5]:.1f}\n")

# --- codon usage across all CDS
cod = Counter()
for f in cds_feats:
    n = str(f.extract(seq))
    cod.update(n[i:i + 3] for i in range(0, len(n) - 2, 3))
top = cod.most_common(15)
plt.figure(figsize=(7, 3.5))
plt.bar([c for c, _ in top], [n for _, n in top])
plt.title("15 most used codons in pPCP1 CDS"); plt.ylabel("count"); plt.tight_layout()
plt.savefig("results/codon_usage.png", dpi=150); plt.close()

open("results/plasmid_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
for r in rows: print(r)
