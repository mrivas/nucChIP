.. methods::

Methods
=======

Enrichment of histone marks
---------------------------

Position-specific properties of the genome, such as chromatin compactness or GC content, results in local under-sampling of MNase ChIP-seq fragments. This process may introduced spurious signals of enrichment as some region may appear highly enriched compared with --otherwise equal-- under-sampled regions. A widely used correction method is to normalized the number of histone reads by a control that is expected to reflect the level of local under-sampling (typically IgG or MNase). That's, on any given nucleosome, :math:`j`, the number of histone reads, :math:`x_j`, is normalized as the ratio:

.. math::

   r_j = \frac{ x_j } { n_j }

where :math:`n_j` is the concomitant number of control reads. 

This method assumes that :math:`x_j` is linearly dependent on :math:`n_j`. However, this may not be a accurate representation of their relation. In particular, is plausible to assume that the proportion of under-sampled genomic regions is more prevalent among nucleosomes with low values of :math:`n_j` as these sets gather the under-sampled nucleosomes from other sets with higher levels of :math:`n_j`. According to this, the values of :math:`x_j` at nucleosomes with low values of :math:`n_j` are a mixture of truly low converage nucleosomes and under-sampled nucleosomes, whereas at sets nucleosomes with larger values of :math:`n_j` they are increasingly representative of non-under-sampled nucleosomes. In consequence, we hypothesized that changes on :math:`x_j` should became less sensitive to changes on :math:`n_j` as :math:`n_j` increases. In other words, :math:`x_j` should follow a saturation curve with respect to :math:`n_j`.

As the sequencing step introduces random noise in both  signal and control libraries (once ready for sequencing, the fragments of a library are chosen randomly according to a binomial distribution) we estimated their relation statistically. Here, we propose to estimate under-sampling level at the nucleosome level as the expected value of histone reads given the value of the control at the same nucleosome. That's, to normalized the :math:`x_j` values as:

.. math::

   r_j = \frac{ x_j } { E(X|n_j)}

where, :math:`E(X|n_j)` is the expected number of histone reads per nucleosome, :math:`X`, given :math:`n_j` concomitant control reads. 

:math:`E(X|n_j)`, in turn, was computed for the observed range of :math:`n_j` using all genomic nucleosomes (iNPs): 

.. math::

   E(X|n_j) = \frac{1}{||J(n_j)||} \sum_{j \in J(n_j)} x_j

Here, :math:`J(n_j)` is the subset of nucleosomes with :math:`n_j` supporting MNase reads. 

By compromising all genomic nucleosomes, :math:`E(X|n_j))` is not only unbias towards any position-specific biological function, but also tailors 1 as a reference point for :math:`r_j`: values of :math:`r_j` above and below 1 can be interpreted as enriched or depleted, respectively, of histone marks.

Additionally, since :math:`E(X|n_j)` is directly proportional to the total number of reads per histone library, the ratios :math:`r_j` are a metric already normalized by library size.

For each MNase ChIP-seq library, after removing outliers (read counts per nucleosome over 99% quantiles) :math:`E(X|n_j)` resulted (Figures 1-4) in monotonic transformations of the number of MNase reads per nucleosome. Interestingly, :math:`E(X|n_j)` saturates in the direction of the :math:`x`-axis, meaning that as the number of supporting of MNase reads gets larger its confounding effect on histone reads gets attenuated. By taking into account this effect, :math:`E(X|n_j)` improves the sensitivity of :math:`r_j` when compared to :math:`n_j` as denominator for :math:`r_j`. This difference is specially important among nucleosomes with large values of :math:`n_j`, where using :math:`n_j` as denominator would over-amplifies the difference between their histone enrichment signals. 

The linear relationship between :math:`n_j` and :math:`x_j` would only holds when the proportion of under-sampled nucleosomes remains constant with respect to :math:`n_j`. If this may be the case for a particular library, :math:`E(X|n_j)` will be simply reduce to :math:`n_j`. Thus, :math:`E(X|n_j)` can be interpreted as general alternative to :math:`n_j` for cases when the proportion of under-sampled nucleosomes doesn't remain constant with :math:`n_j`.

As these results support our hypothesis, we concluded that the ability of :math:`E(X|n_j)` to capture the sensitivity changes of :math:`x_j` with respect to :math:`n_j` can improve the inference of :math:`r_j`.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/17_H3K4me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/H3K4me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/n1_H3K4me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/n2_H3K4me3.expectedCounts.svg
   :width: 45%
Figure 1: Expected number of H3K4me3 reads given supporting MNase reads per nucleosome.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/14_H3K27Ac.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/6_H3K27Ac.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/H3K27Ac.expectedCounts.svg
   :width: 45%
Figure 2: Expected number of H3K27Ac reads given supporting MNase reads per nucleosome.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/4_H3K9me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/9_H3K9me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/H3K9me3.expectedCounts.svg
   :width: 45%
Figure 3: Expected number of H3K9me3 reads given supporting MNase reads per nucleosome.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/12_H3K27me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/5_H3K27me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/H3K27me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/all_nucs/n3_H3K27me3.expectedCounts.svg
   :width: 45%
Figure 4: Expected number of H3K27me3 reads given supporting MNase reads per nucleosome.

Results using only high confidense nucleosomes can be seen here:

.. toctree::

   highQual_nuc.rst
