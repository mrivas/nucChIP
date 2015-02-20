.. methods::

Methods
=======

Enrichment of histone marks
---------------------------

The number of MNase reads supporting the position of a nucleosome is a confounding variable for the number of histone reads overlapping the same nucleosome. On average, the larger the number of supporting MNase reads the larger the enrichment of histone reads. Thus, to filter out spurious signals of enrichment, for any given nucleosome, :math:`j`, we normalized the number of histone reads, :math:`x_j`, as the ratio:

.. math::

   r_j = \frac{ x_j } { E(X|n_j)}

where, :math:`E(X|n_j)` is the expected number of histone reads per nucleosome, :math:`X`, given :math:`n_j` supporting MNase reads. For the observed range of :math:`n_j`, we computed :math:`E(X|n_j)` using all genomic nucleosomes (iNPs, p-values :math:`<=` 0.05) according to: 

.. math::

   E(X|n_j) = \frac{1}{||J(n_j)||} \sum_{j \in J(n_j)} x_j

Here, :math:`J(n_j)` is the subset of nucleosomes with :math:`n_j` supporting MNase reads. 

By compromising all genomic nucleosomes, :math:`E(X|n_j))` is not only unbias towards any position-specific biological function, but also tailors 1 as a reference point for :math:`r_j`: values of :math:`r_j` above and below 1 can be interpreted as enriched or depleted, respectively, of histone marks.

Additionally, since :math:`E(X|n_j)` is directly proportional to the total number of reads per histone library, the ratios :math:`r_j` are a metric already normalized by library size.

For each MNase ChIP-seq library, after removing outliers (read counts per nucleosome over 99% quantiles) :math:`E(X|n_j)` resulted (Figures 1-4) in monotonic transformations of the number of MNase reads per nucleosome. Interestingly, :math:`E(X|n_j)` saturates in the direction of the :math:`x`-axis, meaning that as the number of supporting of MNase reads gets larger its confounding effect on histone reads gets attenuated. By taking into account this effect, :math:`E(X|n_j)` improves the sensitivity of :math:`r_j` compared to the widely used method of using :math:`n_j` as denominator for :math:`r_j`, which scale the values of :math:`x_j` in a rather linear fashion with respect to :math:`n_j`. This difference is specially important among nucleosomes with large values of :math:`n_j`, where using :math:`n_j` as denominator for :math:`r_j` would otherwise over-amplifies the difference between their histone enrichment signals. 

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/methods/17_H3K4me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/H3K4me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/n1_H3K4me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/n2_H3K4me3.expectedCounts.svg
   :width: 45%
Figure 1: Expected number of H3K4me3 reads given supporting MNase reads per nucleosome.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/methods/14_H3K27Ac.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/6_H3K27Ac.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/H3K27Ac.expectedCounts.svg
   :width: 45%
Figure 2: Expected number of H3K27Ac reads given supporting MNase reads per nucleosome.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/methods/4_H3K9me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/9_H3K9me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/H3K9me3.expectedCounts.svg
   :width: 45%
Figure 3: Expected number of H3K9me3 reads given supporting MNase reads per nucleosome.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/methods/12_H3K27me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/5_H3K27me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/H3K27me3.expectedCounts.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/methods/n3_H3K27me3.expectedCounts.svg
   :width: 45%
Figure 4: Expected number of H3K27me3 reads given supporting MNase reads per nucleosome.
