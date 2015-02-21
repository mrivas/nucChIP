.. methods::

Methods
=======

Enrichment of histone marks
---------------------------

Position-specific properties of the genome, such as chromatin compactness or GC content, results in local under-sampling of MNase ChIP-seq fragments. This process may introduced spurious signals of enrichment as some region may appear highly enriched compared with --otherwise equal-- under-sampled regions. A widely used correction method is to normalized the number of histone reads by an estimate of the level of local under-sampling (typically MNase or IgG). That's, on any given nucleosome, :math:`j`, the number of histone reads, :math:`x_j`, is normalized as the ratio:

.. math::

   r_j = \frac{ x_j } { n_j }

where :math:`n_j` is the concomitant number of control reads. 

However, this method doesn't take into account the random noise introduced during the sequencing step. Once ready for sequencing, the fragments of a library are chosen randomly according to a binomial distribution. Thus, any estimate of local under-sampling should be treated statistically. 

Here, we propose to estimate under-sampling at the nucleosome level as the expected value of histone reads given the value of the control at the same nucleosome. That's, the normalized the :math:`x_j` values as:

.. math::

   r_j = \frac{ x_j } { E(X|n_j)}

where, :math:`E(X|n_j)` is the expected number of histone reads per nucleosome, :math:`X`, given :math:`n_j` supporting MNase reads. 

:math:`E(X|n_j)` was computed for the observed range of :math:`n_j` using all genomic nucleosomes (iNPs): 

.. math::

   E(X|n_j) = \frac{1}{||J(n_j)||} \sum_{j \in J(n_j)} x_j

Here, :math:`J(n_j)` is the subset of nucleosomes with :math:`n_j` supporting MNase reads. 

By compromising all genomic nucleosomes, :math:`E(X|n_j))` is not only unbias towards any position-specific biological function, but also tailors 1 as a reference point for :math:`r_j`: values of :math:`r_j` above and below 1 can be interpreted as enriched or depleted, respectively, of histone marks.

Additionally, since :math:`E(X|n_j)` is directly proportional to the total number of reads per histone library, the ratios :math:`r_j` are a metric already normalized by library size.

For each MNase ChIP-seq library, after removing outliers (read counts per nucleosome over 99% quantiles) :math:`E(X|n_j)` resulted (Figures 1-4) in monotonic transformations of the number of MNase reads per nucleosome. Interestingly, :math:`E(X|n_j)` saturates in the direction of the :math:`x`-axis, meaning that as the number of supporting of MNase reads gets larger its confounding effect on histone reads gets attenuated. By taking into account this effect, :math:`E(X|n_j)` improves the sensitivity of :math:`r_j` compared to the widely used method of using :math:`n_j` as denominator for :math:`r_j`. This difference is specially important among nucleosomes with large values of :math:`n_j`, where using :math:`n_j` as denominator for :math:`r_j` would otherwise over-amplifies the difference between their histone enrichment signals. 

The linear relationship between :math:`n_j` and :math:`x_j` would only holds when the proportion of under-sampled nucleosomes remains constant with respect to :math:`n_j`. If this may be the case for a particular library, :math:`E(X|n_j)` will be simply reduce to :math:`n_j`. Thus, :math:`E(X|n_j)` can be interpreted as general alternative to :math:`n_j` for cases when the proportion of under-sampled nucleosomes doesn't remain constant with :math:`n_j`.

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
