Manuscript
==========

Introduction
------------

Explain alternative splicing, why is important

How it's control -> we need to learn more about the genome

What are we doing here, high resolution epigenome

What have we found?

Alternative splicing is a mechanism that produces diverse transcripts of a gene. This happens during the translation of the genetics information into messenger RNA, where portions of the RNA sequence are removed, and the remaining sub-sequences are spliced. Variations in this process result in alternative readouts of the gene. More than 90% of human genes are subject to alternative splicing, and its misregulation is linked to several diseases, including cystic fibrosis, Parkinsonism, and cancer. 

The most common form of alternative splicing is differential inclusion of exons. In eukaryotes, the structure of a typical protein-coding gene is a chain of exons, each flanked by large intronic sequences. Whereas some of these exons are constitutively included in the mature messenger RNA, other are skipped along intronic regions in a tissue or developmental stage-specific manner. The removal of intronic regions and alternatively skipped exons is controlled by the spliceosome, an arrange of five snRNAs and ~300 proteins. The formation of the spliceosome is seeded by trans-acting RNA-binding proteins that recognize cis-acting elements on the messenger RNA. From there, the spliceosome outgrow into a catalytic assembly through a process that is influenced by the local epigenetic landscape in two ways. On one hand, during transcription post-trancriptional histone modifications, such as acetylation, can ease Poll II passage over nucleosomes by destabilizing DNA-histone interactions. Since most eukaryotes genes are spliced cotranscriptionally, faster elongation rates hurdle the formation of the spliceosome by inducing secondary structures on the nascent RNA (Luco) and providing less time for the splicing machinery to recognize splicing signals. On the other hand, core spliceosome components and splicing regulators are recruited by chromatin remodeling proteins that bind specific histone marks cotranscriptionally.

Several studies have consistently showed how alternative splicing is affected by sequences features, such as cis-acting RNA sequence motifs, exon/intron length, propensity to form secondary, and mutations on the DNA sequence. Efforts to bring into the picture the role of chromatin and histone modifications include 

equences to atract the splicing machinery has been well studied. As for the epigenome there is clear evidence relation histone landscape and exon definition, but not of the role of each histone mark. The main gun to explore the effect of epignome on alternative splicing is ChIP-seq. However, association studies are dampened by low resoultion of the data. For instance, the high coverage of H3K4me3 at TSS of expressed genes overshadows it study on downstream exons, even when discarding exon 1. Here, we used a protocol to produce high resoulution ChIP-seq maps of the epigenome. 


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

As we were interested on histone tail marks, we only used reads on the mono-nucleosome range (135-155 bp) to call nucleosome locations. Filtering out reads with mapping qualities below 20, we found 10,468,598 nucleosome locations genome-wide. This amount is coherent with the expected value for the mouse genome: the total number of nucleosomes times the combined length of each nucleosomal DNA (147 nt) and its linker sequence (38 nt as the typical distance between neighbors nucleosomes; :cite:`Jiang2009` ) covered approximately 77% of the mouse genome length (2.5 Gb; :cite:`Waterston2002`).

As shown in Figure :num:`#nuc-widths`, the nucleosome's widths peaks at ~75 nt, which is coherent the length used by iNPS to represent the enrichment signals (to improve the signal over background ratio, iNPS reduces each fragment length to 75 (nt) around their midpoint). The sharp peaks is signal that most nucleosomes are well positioned and isolated --not overlapping flanking nucleosomes. On the other hand, the distance between adjacent nucleosomes (Figure :num:`#nuc-dists`) peaks at ~ 180 (nt), being this coherent with the typical combined length of nucleosomal (~147 nt) and linker DNA segments (~38 nt; :cite:`Jiang2009`).


Histone mapping at nucleosome resolution
****************************************

We generated high resolution genome-wide histone maps. In our protocol, we used micrococcal nuclease (MNase) digestion to produce ChIP-seq (MNChIP-seq) insert fragments at mono-nucleosome size. To benchmark our results, we focused our analyzes on histone marks with biological functions well characterized. We choose  H3K4me3 (2-replicates), H3K27Ac, and H3K9me3 as their role as activator and represor of gene expression has being well established.

All our libraries resulted in strong protection of mono-nucleosome size fragments (Figures :num:`#n1-h3k4me3-frag`, :num:`#m1-h3k4me3-frag`, :num:`#m1-h3k27ac-frag`, :num:`f9-h3k9me3-frag`). Whereas H3K9me3 fragment sized peaked at 173 nt, the other libraies peaked at around 147 nt. 

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

Validation of MNChIP-seq libraries
**********************************

To check whether our MNChIP-seq libraries are a truly reflection of the epigenome, we benchmarked our results against sonicated ChIP-seq. We computed the normalized enrichment of our libraries around the transcription start site of high, medium, low expressed genes, as well as silent genes.

Our results show that all libraries recapitulates their expected profile. The activation marks H3K4me3 and H3K27me3 are enriched on active genes compared to silent genes. However, unlike regular ChIP-seq, the higher resolution of our data shows that is only nucleosome +1 what really makes a difference. Conversely, the represor mark H3K9me3 is depleted in active genes, but depleted on silent genes (Figure :num:`#`).

:num:`m1-h3k4me3-tss`
:num:`n1-h3k4me3-tss`
:num:`m1-h3k27ac-tss`
:num:`f9-h3k9me3-tss`

Histone tails codes alternative splicing 
****************************************

We use RNA-seq data to compute gene expression genome-wide. Using the database Katz database of alternative spliced exons (:cite:`Katz2010`, database mouse mm9 version 1.0) we cluster internal exons (filtered out first and last exon of each gene) into spliced-in and spliced-out if their phi value was greater or lower than 0.7. To filter noisy values we requested at least 10 covering each exon and confidense intervals not wider than 0.2. We also created a database of consititutively expreesed exons by removing all know altervnatively spliced exons from the pool.

To avoid the cofounding effect of gene expression on enrichment of the epigenome, we selected only exons sitting on genes with similar gene expression.

We found that compared to regular ChIP-seq, nucChIP-seq was able to show the structure of exons. Both, H3K4me3, and H3K27Ac clearly show enrichment of spliced-in compared to spliced-out exons. What's more, constitutively expressed genes were alwasy more enriched that spliced-exons, meaning that the presence of both marks promotes the inclusion of the exons in the final transcripts. As for H3K9me3 we found not difference in enrichment between spliced (in and out) and consitituvely expressed genes. 



:num:`m1-h3k4me3-exon-5p`
:num:`m1-h3k4me3-exon-3p`
:num:`n1-h3k4me3-exon-5p`
:num:`n1-h3k4me3-exon-3p`
:num:`m1-h3k27ac-exon-5p`
:num:`m1-h3k27ac-exon-3p`
:num:`f9-h3k9me3-exon-5p`
:num:`f9-h3k9me3-exon-3p`

To test the statistical significance of these trends we call canonical nucleosome positions around exons using the MNase data. The MNase profile was smoothed and consequetives peaks were called the center position of each nucleosomes. Then, for each canonical nucleosomes we asigned the mean enrichment among spliced-in, spliced-out, and consitituvely expressed exons. We estimated the distribution of this test statistics by a bootstrap methods (1500 resamplings)

:num:`m1-h3k4me3-boxplot-5p`
:num:`m1-h3k4me3-boxplot-3p`
:num:`n1-h3k4me3-boxplot-5p`
:num:`n1-h3k4me3-boxplot-3p`
:num:`m1-h3k27ac-boxplot-5p`
:num:`m1-h3k27ac-boxplot-3p`
:num:`f9-h3k9me3-boxplot-5p`
:num:`f9-h3k9me3-boxplot-3p`

We also used the Kolmogorov-Smirnov test to test if the distribution of these statistic was significantly differnt betwee the three types of exons. Our results show that H3K4me3 was significantly differnt between the nucleosomes sitting directly on top of spliced-in and spliced-out exons, but not around them. What's more, all nuclesomes sitting on top of consititutively expressed exons were different compared to spliced-in or spliced.out exons. A similar trend was found for H3K27Ac.

Conversely, we found statistically differences on H3K9me3 enrichment only between spliced-in and spliced-out exons. In both cases, constitutively expressed genes were in a middle point. As it has been previouly reported, enrichment of H3K9me3 enrichment correlates with exon inclusion :cite:`Saint-Andre2011`.

:num:`m1-h3k4me3-pvalues-5p`
:num:`m1-h3k4me3-pvalues-3p`
:num:`n1-h3k4me3-pvalues-5p`
:num:`n1-h3k4me3-pvalues-3p`
:num:`m1-h3k27ac-pvalues-5p`
:num:`m1-h3k27ac-pvalues-3p`
:num:`f9-h3k9me3-pvalues-5p`
:num:`f9-h3k9me3-pvalues-3p`


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
