.. methods::

Methods
=======

Normalization of histone signals
--------------------------------

Position-specific properties of the genome, such as chromatin compactness or GC content, affect the sampling efficiency of MNase ChIP-seq fragments. This may introduced spurious signals of enrichment as, for instance, some region may appear enriched compared with --otherwise equal-- under-sampled regions. A widely used correction method is to normalized the counts of histone reads by the count of a control library (for instance IgG or MNase) which is expected to reflect the local level of sampling efficiency. That's, on any given nucleosome, :math:`j`, the number of histone reads, :math:`x_j`, is normalized as the ratio:

.. math::

   r_j = \frac{ x_j } { f(n_j|x,n) } 

where :math:`f(n_j|x,n)` is a parametric function that computes the sampling efficiency given the number of control reads on the current nucleosome, :math:`n_j`. :math:`f(n_j|x,n)` is parametrized by the vectors :math:`x` and :math:`n`, which contain the genome-wide count per nucleosome of signal and control libraries, respectively. 

Typically, :math:`f(n_j|x,n)` is assumed to be a linear function of :math:`n_j`:

.. math::

   f(n_j|x,n) = n_j \frac{\sum_j x_j}{\sum_j n_j}

In other words, the rate of change (sensitivity) of :math:`f(n_j|x,n)` with respect to :math:`n_j` is assumed to be constant. However, this may not be an accurate model. In particular, it's plausible to assume that :math:`n_j` is a reflexion of the functional importance of a nucleosome. Whereas nucleosomes with low :math:`n_j` may not be functionally relevant but the result of baseline coverage, nucleosomes with large values of :math:`n_j` are more likely to play position-specific functions (such as activator/repressors at TSS, enhancers, etc). This is important as the sensitivity of :math:`f(n_j|x,n)` may not be the same among baseline and selective coverage of nucleosomes. As a consequence, we hypothesized that the sensitivity of :math:`f(n_j|x,n)` with respect to :math:`n_j` is not constant.

To test our hypothesis, we estimated the relation between :math:`f(n_j|x,n)` and :math:`n_j` directly from the data by using the expected number of histone reads per nucleosome, :math:`X`, given :math:`n_j`. 

.. math::

   f(n_j|x,n) = E(X|n_j)

For the observed range of :math:`n_j` (using all genomic nucleosomes defined by iNPs) we computed :math:`E(X|n_j)` as: 

.. math::

   E(X|n_j) = \frac{1}{||J(n_j)||} \sum_{j \in J(n_j)} x_j

Here, :math:`J(n_j)` is the subset of nucleosomes with :math:`n_j` supporting MNase reads. 

By compromising all genomic nucleosomes, :math:`E(X|n_j)` is not only un-bias towards any position-specific biological function, but also tailors 1 as a reference point for :math:`r_j`: values of :math:`r_j` above and below 1 can be interpreted as enriched or depleted, respectively, of histone marks.

Additionally, the distribution of :math:`E(X|n_j)` is dependant on the total number of countes per nucleosome of both signal and control libraries. Using it as the denominator on :math:`r_j` produces a metric already normalized by library sizes.

Results
-------

We tested our hypothesis in all our MNase ChIP-seq libraries, using as control MNase. For each library after removing outliers (read counts per nucleosome over 99% quantiles) :math:`E(X|n_j)` resulted (Figures 1-4) in monotonic transformations of the number of MNase reads per nucleosome. Interestingly, the rate of change of :math:`E(X|n_j)` with regard to :math:`n_j` decreases along the :math:`x`-axis, being this direct supports for our hypothesis.

Based on the underlying assumption of our hypothesis, the decreasing sensitivity may be the result of increasing functional specificity. In this scheme, as nucleosomes with low :math:`n_j` values are considered unlikely to play position-specific biological functions, their histone tails may be indiscriminately modified, making their expected values to closely follow changes on baseline coverage. Conversely, nucleosomes with higher values of :math:`n_j` are more likely to be linked to position-specific functions. As some histone marks codes for opposite functions, these nucleosomes' histone tails are modified in a mutually exclusive fashion throughout the genome. With increasing values of :math:`n_j`, the set of nucleosomes gathering any particular histone mark is reduced. Thus, for a given histone mark, the sensitivity of :math:`E(X|n_j)` with respect of :math:`n_j` is dampened by an increasing selectivity in the modification of nucleosomes' histone tails. 

By taking into account this changes in sensitivity, :math:`E(X|n_j)` improves the measurement of :math:`r_j` when compared to linear transformation of :math:`n_j` as denominator for :math:`r_j`. This difference is specially important among nucleosomes with large values of :math:`n_j`, where using the later method would over-estimate the sensitivity of :math:`f(n_j|x,n)` with respect to :math:`n_j`, resulting in artificially larger differences in enrichment.  

The linear relationship between :math:`f(n_j|x,n)` and :math:`n_j` would only holds if the proportion of under-sampled nucleosomes remains constant with respect to :math:`n_j`. If this may be the case for a particular library, :math:`E(X|n_j)` will be simply reduce to :math:`n_j`. Thus, :math:`E(X|n_j)` can be interpreted as general formulation of :math:`f(n_j|x,n)`.

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
