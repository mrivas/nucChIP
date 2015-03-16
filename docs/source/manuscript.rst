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

Nucleosome positions
********************

Unlike our libraries, the MNase data produced by Carone group (:cite:`Carone2014`, 8_mnase ) was not purified against sub-nucleosomal size fragments. After mapping the reads against mouse genome, the MNase data exhibited strong protection of mononucleosome and subnucleosomal size fragments. Its fragments' distribution had several peaks (Figure :num:`#mnase-fragments`), being the highest on the mononucleosomal range (144 bp). 

As we were interested on histone tail marks, we only used reads to the mono-nucleosome range (135-155 bp) to call nucleosome locations. Filtering out reads with mapping qualities below 20, we found 10,468,598 nucleosome locations genome-wide. This amount is coherent with the expected value for the mouse genome. The total number of nucleosomes times the combined length of each nucleosomal DNA (147 nt) and its linker sequence (38 nt as the typical distance between neighbors nucleosomes; :cite:`Jiang2009` ) covered approximately 77% of the mouse genome length (2.5 Gb; :cite:`Waterston2002`).

As shown in Figure :num:`#nuc-widths`, the nucleosome's widths peaks at ~75 nt, which is coherent the length used by iNPS to represent the enrichment signals (to improve the signal over background ratio, iNPS reduces each fragment length to 75 (nt) around their midpoint). The sharp peaks is signal that most nucleosomes are well positioned and isolated --not overlapping flanking nucleosomes. On the other hand, the distance between adjacent nucleosomes (Figure :num:`#nuc-dists`) peaks at ~ 180 (nt), being this coherent with the typical combined length of nucleosomal (~147 nt) and linker DNA segments (~38 nt; :cite:`Jiang2009`).

Histone mapping at nucleosome resolution
****************************************

We generated high resolution genome-wide histone maps. In our protocol, we used micrococcal nuclease (MNase) digestion to produce ChIP-seq (MNChIP-seq) insert fragments at mono-nucleosome size. To benchmark our results, we focused our analyzes on histone marks with biological functions well characterized. We choose  H3K4me3 (2-replicates), H3K27Ac, and H3K9me3 as their role as activator and represor of gene expression has being well established.

All our libraries resulted in strong protection of mono-nucleosome size fragments (Figures :num:`#n1-h3k4me3-frag`, :num:`#m1-h3k4me3-frag`, :num:`#m1-h3k27ac-frag`, :num:`9-h3k9me3-frags`). Whereas H3K9me3 fragment sized peaked at 173 nt, the other libraies peaked at around 147 nt. 

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

We tested our hypothesis in all our MNase ChIP-seq libraries, and using as control MNase. For each library, after removing outliers (read counts per nucleosome over 99% quantiles), :math:`E(X|n_j)` resulted (Figures :num:`#m1-h3k4me3-exp`, :num:`#n1-h3k4me3-exp`, :num:`#m1-h3k27ac-exp`, and :num:`#9-h3k9me3-exp`) in monotonic transformations of the number of MNase reads per nucleosome. Interestingly, the rate of change of :math:`E(X|n_j)` with regard to :math:`n_j` decreases along the :math:`x`-axis. This is direct support for our hypothesis.

Variable sensitivity may be the result of differences in functional specificity. To understand this idea, first we have to realize that when using MNase as control, the sensitivity of :math:`E(X|n_j)` with respect to :math:`n_j` can be interpreted as the average proportion of nucleosomes on a particular position having the corresponding histone tail mark as a result of background coverage. Seen from this perspective, the question is: why at low :math:`n_j` nucleosomes show a higher proportion of background histone tails modifications than at higher :math:`n_j`. According to the underlying assumption of our hypothesis, as nucleosomes with low :math:`n_j` values are unlikely to play position-specific biological functions, their histone tails may be indiscriminately modified, resulting in histone baseline coverage, on average, similar to nucleosome coverage. As a result, on these nucleosomes :math:`E(X|n_j)` closely follows changes on :math:`n_j`. Conversely, this effect is dampened at larger :math:`n_j` values, where baseline coverage of the signal became less prevalent as proportion of nucleosome coverage.

By taking into account this changes in sensitivity, :math:`E(X|n_j)` improves the measurement of :math:`r_j` when compared to linear transformation of :math:`n_j` as denominator for :math:`r_j`. This difference is specially important among nucleosomes with large values of :math:`n_j`, where using the later method would over-estimate the sensitivity of :math:`f(n_j|x,n)` with respect to :math:`n_j`, resulting in artificially larger differences in enrichment.  

The linear relationship between :math:`f(n_j|x,n)` and :math:`n_j` would only holds if the proportion of position-specific nucleosomes remains constant with respect to :math:`n_j`. If this may be the case for a particular library, :math:`E(X|n_j)` will be simply reduce to a linear trend. Thus, :math:`E(X|n_j)` can be interpreted as general formulation of :math:`f(n_j|x,n)`.

Validation of MNChIP-seq libraries
**********************************

To check whether our MNChIP-seq libraries are a truly reflection of the epigenome, we benchmarked our results against regular ChIP-seq. We computed the normalized enrichment of our libraries around the transcription start site of high, midium, low expressed genes, as well as silent genes.

Our results show that all libraries recapitulates their expected profile. The activation marks H3K4me3 and H3K27me3 are enriched on active genes compared to silent genes. However, unlike regular ChIP-seq, the higher resolution of our data shows that is only nucleosome +1 what really makes a difference. Conversely, the represor mark H3K9me3 is depleted in active genes, but depleted on silent genes (Figure :num:`#`).

Histone tails codes alternative splicing 
****************************************

:ref:`coverage_exons`

Statisticall validation
***********************

:ref:`countsPerNuc`

Discussion
----------

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

Procedure
1. Map data to mm9 with bowtie2, default parameters
2. Remove duplicates with picards tools
3. Count reads per nucleosome, getCounts
4. Compute expected values, with R script
5. Plot coverage per nucleosome,
6. Plot fragment size distribution, vPlot2

Bibliography
============

.. bibliography:: Mendeley.bib
   :style: plain
