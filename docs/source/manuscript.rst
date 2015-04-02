Manuscript
==========

Introduction
------------

Alternative splicing is a mechanism that produces diverse transcripts of a gene, being exon skipping its most common form. In eukaryotes, the structure of a typical protein-coding gene is a chain of exons, each flanked by large intronic sequences. Whereas some of these exons are constitutively included in the mature messenger RNAs, other are removed along intronic regions in a temporal and spatial-specific manner by the spliceosome: an arrange of five snRNAs and ~300 proteins :cite:`Hoskins2012`. Alternative splicing affects more than 90% of human genes :cite:`Pan2008`, and its misregulation is linked to several diseases :cite:`Wang2007` :cite:`Xiong2014`.

The spliceosome delivers specificity by trans-acting RNA-binding proteins that recognize cis-acting elements on the messenger RNA :cite:`Chen2009` :cite:`Barash2010`. However, in recent years researchers have started to paint a picture where local histone modifications constitute an additional regulatory layer :cite:`Luco2011` :cite:`Zhou2014`. One mechanism involves modulation of mRNA elongation rate by local chromatin state. Since most eukaryotes genes are spliced cotranscriptionally :cite:`CarrilloOesterreich2010`, faster elongation rates hurdle the formation of the spliceosome by inducing secondary structures on the nascent RNA :cite:`Buratti2004` and reducing the time-window for the recognition of splicing signals :cite:`DelaMata2003`. Interestingly,  nucleosomes positioning correlates well with exon boundaries :cite:`Schwartz2009` :cite:`Tilgner2009`, and acetylation of their histones tails (which ease Poll II passage by destabilized DNA-histone interactions) reduces exon inclusion :cite:`Li2007`. A more direct involvement of histone marks, on the other hand, have come from gene-specific studies showing that the recruitment of core spliceosome components and splicing regulators is mediated by chromatin remodeling proteins that bind to specific histone marks :cite:`Sims2007` :cite:`Batsche2006` :cite:`Saint-Andre2011` :cite:`Luco2010`.

Using ChIP-seq data to correlate chromatin states to splicing outputs, several studies (:cite:`Ernst2010` :cite:`Ernst2011` :cite:`Dhami2010` :cite:`Enroth2012` :cite:`Pan2008` :cite:`Shindo2013` :cite:`Zhou2012` :cite:`Podlaha2014`) have reflected the role of histone marks on alternative splicing at genome-wide level. However, a main limiting factor of these studies has been the resolution of the ChIP-seq libraries. On average, exons not only overlap nucleosomes :cite:`Schwartz2009` :cite:`Tilgner2009`, but their typical length is similar to the number of DNA base pairs (147) needed to wrap a single nucleosome :cite:`Luco2011`. However, typical ChIP-seq protocols produces DNA fragments sizes on the range ~200-700 bp :cite:`Zhang2011` that allows the characterization histone marks over stretches of DNA spanning several nucleosomes, making difficult to resolve histone signals at exon-intron junctions. To overcome this barrier, here we used a high resolution paired-end ChIP-seq protocol to produced histone postranscriptional modification maps at single nucleosome resolution. Along with gene expression data and newly developed computational methods, we found clear correlations between enrichment of H3K4me3 and H3K27Ac on spliced-in over spliced-out exons, and weaker correlations between H3K9me3 enrichment and exon inclusion of alternatively skipped exons. Using a publicly available MNase-seq library we found no important difference in nucleosome enrichment between spliced-in, spliced-out and constitutively expressed exons.

We found that:
1. Nucleosomes are not differentially distributed between spliced-in and spliced-out exons.
1. H3K4me3 is enriched on included exons, being its enrichment decrecing from constitutively expressed, to spliced-in, and to spliced-out.
2. H3K9me3 differentiate aternatively spliced exons from constitutive exons. Using constitytive exons as reference we found that on alternative spliced-exons, enrichment and depletion of H3K9me3 produces exon inclussion and exclusion, respectively.
3. MNase data shows that there is extra-spacing before and after exons, must likely to give splice for the splicing machinery to bind the DNA sequence
4. Compared to constitutive exons, nucleosomes on alternatively spliced exons are weakly bound to the genome. This is coherent with the kinetic theory of alternative splicing.
5. There isn't evidence of proteins binding intron-exons junctions. The splicesosome must likely binds the mRNA, not the DNA.

Results
-------

Nucleosomes genome-wide positions
*********************************

First, we set out to determine nucleosomes' positions genome-wide using Carone et al :cite:`Carone2014` MNase-seq data. After mapped against the mouse genome, MNase-seq reads exhibited strong protection of mononucleosome size fragments, with 144 bp as the main fragment value (supplementary Figure 1). The MNase-seq data also showed protection of sub-nucleosomal size fragments as Carone's et al protocol :cite:`Carone2014` was intended to capture the footprints of molecules smaller than nucleosomes. Consequently, we only used MNase-seq reads with fragments on the range 135-155 bp as nucleosomes' footprints.

We used iNPs (a nucleosome-positioning algorithm; :cite:`Chen2014` ) to estimate nucleosomes' starting and ending coordinates based on the nucleosomes' footprints. We found 10,468,598 nucleosome locations genome-wide. 
As shown in Figure :num:`#nuc-widths`, the nucleosome's widths peaks at ~75 nt, which is coherent the length used by iNPS to represent the enrichment signals (to improve the signal over background ratio, iNPS reduces each fragment length to 75 (nt) around their midpoint). The sharp peaks is signal that most nucleosomes are well positioned and isolated --not overlapping flanking nucleosomes. On the other hand, the distance between adjacent nucleosomes (Figure :num:`#nuc-dists`) peaks at ~ 180 (nt), being this coherent with the typical combined length of nucleosomal (~147 nt) and linker DNA segments (~38 nt; :cite:`Jiang2009`).
This amount is concordant with the size of mouse genome (Figure 2): the total number of nucleosomes times the combined length of each nucleosomal DNA (147 nt) and its linker sequence (38 nt as the typical distance between neighbors nucleosomes; :cite:`Jiang2009` ) covered approximately 77% of the mouse genome length (2.5 Gb; :cite:`Waterston2002`). 

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/manuscript/fragDist.svg
   :width: 45%
Figure 1: Fragment size distribution of MNase-seq and MNChIP-seq libraries. All libraries strong protection of mononucleosomal size fragments.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/manuscript/nucLocations.svg
   :width: 45%
Figure 2: Nucleosome properties.

Histone mapping at single nucleosome resolution
***********************************************

Using MNChIP-seq, we generated single-nucleosome resolution maps of three histone marks: H3K3me3 (two replicates), H3K27Ac (two replicates), and H3K9me3. Once mapped to the mouse genome, all libraries produced fragment sizes distribution on the mononucleosomal range, with their mean values ranging from 144 to 167 bp (supplementary Figure 1). 

Since position-specific properties of the genome, such as chromatin compactness or GC content, can produce biases on the sampling of DNA fragments, we develop a computational method (see Methods) to normalized the counts of histone footprints per nucleosome (Figure 2). In our method, the histone (:math:`x_j`) and nuclesome (:math:`n_j`) footprints are counted over each nucleosome interval (:math:`j`) they happen to overlap. Then, :math:`x_j` is normalized as the ratio (:math:`r_j`):

.. math::
   
   r_j = \frac{ x_j } { E(X|n_j) } 

where :math:`E(X|n_j)` is the expected value of histone footprints (:math:`X`) given :math:`n_j`.

Thus normalized, we benchmarked our MNChIP-seq results by computing signal enrichment at the transcription start site (TSS) of silent and active genes. Our results (Figure 3) recapitulated prior studies on the role of H3K3me3 and H3K27Ac as activator, and H3K9me3 as repressor of gene expression. Both H3K4me3 and H3K27Ac show strong enrichment on active genes compared to silent genes. On the contrary, H3K9me3 is prevalent on silent genes, with lower enrichment signals on expressed genes.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/manuscript/cartoon.svg
   :width: 45%
Figure 3: Mapping of chromatin states at single nucleosome resolution. (A) Chromatin was digested by MNase to produce footprints in the nucleosomal range (~147 bp). The digested chromatin was then follow one of two paths. On the first (B), immunoprecipitation is used to filter histones having the marks: H3K4me3, H3K27Ac, or H3K9me3. The filtered fragments, are then (D) mapped against mouse genome to reveal the histone footprints. On the other path, (C) the digested chromatin fragments are directly mapped against the genome to detect the nucleosomes' footprints, which are eventually used (E) to detect the nucleosome positions (iNPs). At the final step (F), we counted the number of histone and nucleosome footprints over each nucleosome location, defined as +/-75 bp area around a nucleosome position. The histone counts are normalized (see Methods) over the nucleosome counts to correct for possible position-specific biases in the sampling process.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/manuscript/tss.svg
   :width: 45%
Figure 4: Normalized signal enrichment at TSS. Genes were clustered according to their gene expression. Silent genes correspond to genes with no detectable expression. Low, middle, and high expression genes were obtained from the 1st, 5th, and 9th quantiles of genes with expression above zero.  

Histone tails codes alternative splicing 
****************************************

Having showed the consistency of our data with previous studies, we next surveyed the chromatin structure of constitutively and alternatively spliced exon on mouse embryonic stem cells. We use RNA-seq data to classify Katz et al :cite:`Katz2010` database  (mouse mm9 version 1.0) of alternative spliced exons into spliced-in and spliced-out. We also created a database of internal (first and last exons of each transcript were filtered out) consititutively expreesed exons by removing all know alternatively spliced exons from the know pool of exons. To avoid the confounding effect of gene expression on enrichment of the histone signals, we selected only exons sitting on genes with similar gene expression (Figure).

The mean coverage of histone enrichment among all three types of exons (Figure) clearly showed that exons are marked by histones. What's more, both H3K4me3, and H3K27Ac were strongly enrichment on spliced-in compared to spliced-out exons. Constitutively expressed genes, on the other hand, were always more enriched that spliced-exons, suggesting that their presence promotes the inclusion of the exons in the final transcripts. As for H3K9me3, the differences in enrichment between spliced-in and spliced-out was weaker but still observable.

To test the statistical significance of these trends we called canonical nucleosome positions around exons using the MNase data. The MNase profile was smoothed and consequetives peaks were called the center position of each nucleosomes. Then, for each canonical nucleosomes we asigned the mean enrichment among spliced-in, spliced-out, and consitituvely expressed exons. We estimated the distribution of this test statistics by a bootstrap methods (1500 resamplings)

We also used the Kolmogorov-Smirnov test to test if the distribution of these statistic was significantly differnt betwee the three types of exons. Our results show that H3K4me3 was significantly differnt between the nucleosomes sitting directly on top of spliced-in and spliced-out exons, but not around them. What's more, all nuclesomes sitting on top of consititutively expressed exons were different compared to spliced-in or spliced.out exons. A similar trend was found for H3K27Ac.

Conversely, we found statistically differences on H3K9me3 enrichment only between spliced-in and spliced-out exons. In both cases, constitutively expressed genes were in a middle point. As it has been previouly reported, enrichment of H3K9me3 enrichment correlates with exon inclusion :cite:`Saint-Andre2011`.



Discussion
----------

We developed several open-source tool for analysis, and visualization of MNChIP-seq data.

H3K4me3 and H3K9me3 both enriched at nuc 1, but whereas H3K4me3 is enriched at nuc 2 on high expression genes, H3K9me3 is the opossite.

MNase digested ChIP-seq improves resolution over sonicated ChIP-seq
MNase digested ChIP-seq are coherent with sonicated ChIP-seq

Empty spaces are not bound by proteins (wide range MNase show so)

H3K4me3 is proportional to exon inclusion

MNase and H3K9me3 are slightly enriched on spliced-in exons



Materials and Methods
---------------------

Mapping (bowtie2 default parameters)
Removed duplicates (picard tools)
Gene expression (cufflinks)
Discovery of nucleosomes (iNPs, MNase)

Normalization of histone enrichment signals
-------------------------------------------

Position-specific properties of the genome, such as chromatin compactness or GC content, affect the sampling efficiency of MNase ChIP-seq fragments. This may introduced spurious signals of enrichment as, for instance, some region may appear enriched compared with --otherwise equal-- under-sampled regions. A widely used correction method is to normalized the counts of histone reads by the count of a control library (for instance IgG or MNase). That's, on any given nucleosome, :math:`j`, the number of histone reads, :math:`x_j`, is normalized as the ratio:

.. math::

   r_j = \frac{ x_j } { f(n_j|x,n) } 

where the function :math:`f(n_j|x,n)` computes the sampling efficiency given the number of control reads on the current nucleosome, :math:`n_j`. :math:`f(n_j|x,n)` is parametrized by the vectors :math:`x` and :math:`n`, which contain the genome-wide count per nucleosome of signal and control libraries, respectively. 

Typically, :math:`f(n_j|x,n)` is assumed to be a linear function of :math:`n_j`:

.. math::

   f(n_j|x,n) = n_j \frac{\sum_j x_j}{\sum_j n_j}

In other words, the rate of change (sensitivity) of :math:`f(n_j|x,n)` with respect to :math:`n_j` is assumed to be constant. However, this may not be an accurate model. In particular, it's plausible to assume that :math:`n_j` is a reflexion of the functional importance of a nucleosome. Whereas nucleosomes with low :math:`n_j` may not be functionally relevant but the result of baseline coverage, nucleosomes with large values of :math:`n_j` are more likely to play position-specific functions (such as activator/repressors at TSS, enhancers, etc). This is important as the sensitivity of :math:`f(n_j|x,n)` may not be the same among baseline and and function-specific nucleosomes. Here, we hypothesized that the sensitivity of :math:`f(n_j|x,n)` with respect to :math:`n_j` is not constant.

To test our hypothesis, we took into account stochastic variations on the read counts of signal and control libraries by estimating the relation between :math:`f(n_j|x,n)` and :math:`n_j` as the expected number of histone reads per nucleosome, :math:`X`, given :math:`n_j`. 

.. math::

   f(n_j|x,n) = E(X|n_j)

For the observed range of :math:`n_j` (using all genomic nucleosomes defined by iNPs) we computed :math:`E(X|n_j)` as: 

.. math::

   E(X|n_j) = \frac{1}{||J(n_j)||} \sum_{j \in J(n_j)} x_j

Here, :math:`J(n_j)` is the subset of nucleosomes with :math:`n_j` control reads.

By compromising all genomic nucleosomes, :math:`E(X|n_j)` is not only un-bias towards any position-specific biological function, but also tailors 1 as a reference point for :math:`r_j`; values of :math:`r_j` above and below 1 can be interpreted as enriched or depleted, respectively, of histone marks.


Additionally, the distribution of :math:`E(X|n_j)` is dependant on the total number of counts per nucleosome of both signal and control libraries. Using it as the denominator on :math:`r_j` produces a metric already normalized by library sizes.

We tested our hypothesis in all our MNase ChIP-seq libraries, and using as control MNase. For each library, after removing outliers (read counts per nucleosome over 99% quantiles), :math:`E(X|n_j)` resulted (Figures :num:`#m1-h3k4me3-exp`, :num:`#n1-h3k4me3-exp`, :num:`#m1-h3k27ac-exp`, and :num:`#f9-h3k9me3-exp`) in monotonic transformations of the number of MNase reads per nucleosome. Interestingly, the rate of change of :math:`E(X|n_j)` with regard to :math:`n_j` decreases along the :math:`x`-axis. This is direct support for our hypothesis.

Variable sensitivity may be the result of differences in functional specificity. To understand this idea, first we have to realize that when using MNase as control, the sensitivity of :math:`E(X|n_j)` with respect to :math:`n_j` can be interpreted as the average proportion of nucleosomes on a particular position having the corresponding histone tail mark as a result of background coverage. Seen from this perspective, the question is: why at low :math:`n_j` nucleosomes show a higher proportion of background histone tails modifications than at higher :math:`n_j`. According to the underlying assumption of our hypothesis, as nucleosomes with low :math:`n_j` values are unlikely to play position-specific biological functions, their histone tails may be indiscriminately modified, resulting in histone baseline coverage, on average, similar to nucleosome coverage. As a result, on these nucleosomes :math:`E(X|n_j)` closely follows changes on :math:`n_j`. Conversely, this effect is dampened at larger :math:`n_j` values, where baseline coverage of the signal became less prevalent as proportion of nucleosome coverage.

By taking into account this changes in sensitivity, :math:`E(X|n_j)` improves the measurement of :math:`r_j` when compared to linear transformation of :math:`n_j` as denominator for :math:`r_j`. This difference is specially important among nucleosomes with large values of :math:`n_j`, where using the later method would over-estimate the sensitivity of :math:`f(n_j|x,n)` with respect to :math:`n_j`, resulting in artificially larger differences in enrichment.  

The linear relationship between :math:`f(n_j|x,n)` and :math:`n_j` would only holds if the proportion of position-specific nucleosomes remains constant with respect to :math:`n_j`. If this may be the case for a particular library, :math:`E(X|n_j)` will be simply reduce to a linear trend. Thus, :math:`E(X|n_j)` can be interpreted as general formulation of :math:`f(n_j|x,n)`.


Procedure
1. Map data to mm9 with bowtie2, default parameters
2. Remove duplicates with picards tools
3. Count reads per nucleosome, getCounts
4. Compute expected values, with R script
5. Plot coverage per nucleosome,
6. Plot fragment size distribution, vPlot2

We compute histone enrichment at single nucleosomes
1. Map MNase-seq and MNChIP-seq against mm9
1. We call nucleosomes positions from nucleosome footprints
2. We counted histone and nucleosome footprints per nucleosome
3. We corrected MNase cofounding effect
4. We used RNA-seq to classify exons as spliced-in, spliced-out and constitutive. On each class we computed histone enrichment.
4. We found differential enrichment between spliced-in, spliced, out and constitutive among H3K4me3 and H3K27Ac. H3K9me3 didn't show strong differences

Figures:
Cartoon
Fragment size distributions
Nucleosomes locations
Fragment sizes on exons (good distribution of sizes)
Expected values
Enrichment on exons
Genome browser example to validate our results on TSS and exons
Bootstrap distributions of test-statistics
P-values among distributions

Supplemental figures:
Overlapp analysis among H3K4me3 and H3K27Ac (Venn diagrams)

Bibliography
============

.. bibliography:: Mendeley.bib
   :cited:
   :style: plain
