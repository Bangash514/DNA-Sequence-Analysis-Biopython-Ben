# DNA-Sequence-Analysis-Biopython-Ben

**DNA sequence analysis with Biopython: a bacterial plasmid genome and a plant (Opuntia) phylogeny**

Two small, reproducible bioinformatics analyses on real public sequence data:

1. **Plasmid pPCP1 genome analysis** (*Yersinia pestis*, NCBI NC_005816, 9,609 bp): GC content, GC skew, six-frame ORF scan compared with the official annotation, codon usage and protein properties.
2. **Plant mini phylogeny** (*Opuntia*, prickly pear cactus): distance matrix and UPGMA / Neighbour-Joining trees from a 7-sequence alignment.

> Learning and portfolio project. It demonstrates a workflow and is not a research result.

## Quick start
```bash
git clone https://github.com/Bangash514/DNA-Sequence-Analysis-Biopython-Ben.git
cd DNA-Sequence-Analysis-Biopython-Ben
pip install -r requirements.txt
python src/01_plasmid_analysis.py
python src/02_opuntia_phylogeny.py
```
Run the commands from the project root. Outputs are written to `results/`.

## Methods
- **GC content / GC skew:** sliding window of 300 bp, step 50 bp. GC skew is (G - C) / (G + C); the cumulative curve is plotted.
- **ORF scan:** all six reading frames, ATG to stop, at least 100 amino acids. Results are compared with the annotated CDS by their stop codon position and strand.
- **Protein properties:** length, molecular weight and isoelectric point (pI) of each annotated protein, using Biopython `ProteinAnalysis`.
- **Codon usage:** counted over all annotated CDS.
- **Phylogeny:** identity distance matrix from a Clustal W alignment, then UPGMA and Neighbour-Joining trees (Biopython `DistanceTreeConstructor`).

## Results
**Plasmid pPCP1** (files: `results/plasmid_summary.txt`, `results/cds_table.tsv`)
- Length 9,609 bp, overall GC 45.26%.
- 10 annotated CDS, 7 of them at least 100 aa. The ORF scan (ATG start, at least 100 aa) found 7 ORFs, and all 7 end at the same stop codon as an annotated gene.
- The scan does not recover the annotated 123-aa hypothetical protein, which starts with GTG rather than ATG. One found ORF (reverse strand, ending where the annotated 90-aa hypothetical protein ends) is longer than the annotation because its in-frame ATG lies further upstream. This shows why start-codon choice matters for ORF calling.
- Window GC content varies from about 31% to 60% along the plasmid. Such variation can indicate regions with different origins, but this analysis does not test that.

Figures: `results/gc_skew.png`, `results/codon_usage.png`.

| Gene | Product | aa | MW (kDa) | pI | GC % |
|---|---|---|---|---|---|
| - | putative transposase | 340 | 39.9 | 9.63 | 52.8 |
| - | transposase/IS protein | 260 | 29.4 | 9.22 | 50.4 |
| rop | putative replication regulatory protein | 64 | 7.3 | 4.51 | 50.3 |
| - | hypothetical protein | 123 | 14.5 | 11.93 | 58.6 |
| pim | pesticin immunity protein | 145 | 16.3 | 9.52 | 32.9 |
| pst | pesticin | 357 | 40.0 | 5.78 | 39.2 |
| - | hypothetical protein | 138 | 15.8 | 10.11 | 49.2 |
| pla | outer membrane protease | 312 | 34.6 | 5.88 | 41.9 |
| - | putative transcriptional regulator | 99 | 11.1 | 5.85 | 45.7 |
| - | hypothetical protein | 90 | 10.5 | 9.63 | 44.0 |

**Opuntia alignment** (files: `results/distance_matrix.txt`, `results/tree_upgma.png`, `results/tree_nj.png`, `results/tree_*.nwk`)
- 7 sequences, 156 alignment columns, 18 variable columns.
- Identity distances range from 0.6% to 9.6%.
- Caveat: 156 bp is far too short for reliable phylogenetic conclusions, and no branch support (bootstrap) is computed.

## Limitations
- The input data are small test records, not a full genome or a large dataset.
- The ORF scan uses an ATG-only start and a fixed length cutoff, so alternative starts (GTG, TTG) and genes shorter than 100 aa are not detected.
- Distance trees on identity distance are a teaching baseline; model-based methods would be used for research.

## Possible extensions
Bootstrap support, longer sequences such as chloroplast genomes, BLAST searches, and RNA-seq differential expression.

## Data and credits
- `data/NC_005816.gb` (GenBank record) and `data/opuntia.aln` (Clustal W alignment, GenBank AF191658-AF191665) come from the test data of the Biopython repository, https://github.com/biopython/biopython (Biopython licence). The original sequences are from NCBI GenBank.
- Cock P.J.A. et al., *Biopython: freely available Python tools for computational molecular biology and bioinformatics*, Bioinformatics, 2009.

## Author
BEN, Software Engineer. GitHub: [Bangash514](https://github.com/Bangash514)
