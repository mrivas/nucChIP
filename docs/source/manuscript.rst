Manuscript
==========

Introduction
------------

We have used MNase digested ChIP-seq libraries to produce high definition histone maps to discover the function of H3K4me3 and H3K9me3 on alterntive splicing
Unlike sonicated ChIP-seq libraries our libraries avoid the overshadow that TSS produces on internal exonic regions.
This allow us to see for the first time histone modifications at the nucleosome level.

There is evidence for the involvement of H3K4me3 and H3K9me3 as regulators of alternative splicing. However, until now there isn't a genome-wide study of this effect. Here, we are doing so.

There is evidence of a sequence code for alternative splicing. However, unlike sequence, alternative splicing is tissue specific. On the other hand, the epigenome is reversible and tissue specific, being a natural candidate for splicing code.

We found that:
1. H3K4me3 is enriched on included exons, being its enrichment decrecing from constitutively expressed, to spliced-in, and to spliced-out.
2. H3K9me3 differentiate aternatively spliced exons from constitutive exons. Using constitytive exons as reference we found that on alternative spliced-exons, enrichment and depletion of H3K9me3 produces exon inclussion and exclusion, respectively.
3. MNase data shows that there is extra-spacing before and after exons, must likely to give splice for the splicing machinery to bind the DNA sequence
4. Compared to constitutive exons, nucleosomes on alternatively spliced exons are weakly bound to the genome. This is coherent with the kinetic theory of alternative splicing.
5. There isn't evidence of proteins binding intron-exons junctions. The splicesosome must likely binds the mRNA, not the DNA.


Results
-------

Libraries: MNase (public), H3K4me3, H3K27Ac, H3K9me3, RNA-seq(public)

Fragment sizes

Nucleosome locations

MNase digested ChIP-seq improves resolution over sonicated ChIP-seq
MNase digested ChIP-seq are coherent with sonicated ChIP-seq

Empty spaces are not bound by proteins (wide range MNase show so)

H3K4me3 is proportional to exon inclusion

MNase and H3K9me3 are slightly enriched on spliced-in exons

Procedure
1. Map data to mm9 with bowtie2, default parameters
2. Remove duplicates with picards tools
3. Count reads per nucleosome, getCounts
4. Compute expected values, with R script
5. Plot coverage per nucleosome,
6. Plot fragment size distribution, vPlot2

Discussion
----------

Materials and Methods
---------------------

Discovery of nucleosomes (iNPs, MNase)
Normalization of histone enrichment signals

