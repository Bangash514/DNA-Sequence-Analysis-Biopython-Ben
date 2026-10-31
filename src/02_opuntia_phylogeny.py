# BEN 2026-11-01
# Change aim: Build UPGMA and Neighbour-Joining trees from an Opuntia sequence alignment.
# File: 02_opuntia_phylogeny.py
"""Mini phylogeny of 7 Opuntia (prickly pear, plant) sequences using a pre-made
Clustal W alignment (Biopython test data, GenBank AF191658-AF191665)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

os.makedirs("results", exist_ok=True)
aln = AlignIO.read("data/opuntia.aln", "clustal")
for r in aln:
    r.id = r.id.split("|")[3]            # short accession, e.g. AF191659.1
print(f"{len(aln)} sequences, {aln.get_alignment_length()} alignment columns")
var = sum(len(set(aln[:, i])) > 1 for i in range(aln.get_alignment_length()))
print("variable columns:", var)

dm = DistanceCalculator("identity").get_distance(aln)
open("results/distance_matrix.txt", "w").write(str(dm) + f"\nvariable columns: {var}\n")
con = DistanceTreeConstructor()
for name, tree in (("upgma", con.upgma(dm)), ("nj", con.nj(dm))):
    fig = plt.figure(figsize=(6, 3.5)); ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title(f"Opuntia - {name.upper()} tree (identity distance)")
    plt.tight_layout(); plt.savefig(f"results/tree_{name}.png", dpi=150); plt.close()
    Phylo.write(tree, f"results/tree_{name}.nwk", "newick")
print(dm)
