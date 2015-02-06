.. _AS_exons:

Alternatively spliced exons
===========================

We tested whether histone marks and nucleosome enrichment correlates with splicing output, we separated nucleosomes into two sets:

   * nucleosomes overlapping spliced-in exons, and
   * nucleosomes overlapping spliced-out exons

For both types of exons, we requiered that each nucleosome location should be supported by at least 3 MNase reads.

Then, we counted for each nucleosome the number of MNase ChIP-seq and MNAsed reads on them. We also included sonicated ChIP-seq libraries as references. For each histone marks/nucleosome  library we tested whether the distribution on spliced-in vs spliced-out exons differ significantly using the Kolmogorov-Smirnov test.

The results, Figures 1 to 5, show that:

	* H3K4me3 are not significantly different on sonicated ChIP-seq but they are on MNase ChIP-seq libraries
	* H3K27Ac are significantly different on both MNase ChIP-seq and sonicated ChIP-seq libraries
	* H3K9me3 are not significantly different on sonicated ChIP-seq, and results are mixed for MNase ChIP-seq libraries
	* H3K27me3 are significatnly different on MNase ChIP-seq but not sonicated ChIP-seq libraries
	* MNase library is not significantly different.


.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/17_H3K4me3.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/n1_H3K4me3.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/n2_H3K4me3.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/H3K4me3.ks_test.svg
   :width: 45 %
Figure 1: H3K4me3


.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/6_H3K27Ac.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/14_H3K27Ac.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/H3K27Ac.ks_test.svg
   :width: 45 %
Figure 1: H3K27Ac

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/4_H3K9me3.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/9_H3K9me3.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/H3K9me3.ks_test.svg
   :width: 45 %
Figure 1: H3K9me3

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/5_H3K27me3.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/12_H3K27me3.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/n3_H3K27me3.ks_test.svg
   :width: 45 %
.. image:: https://132.239.135.28/public/nucChIP/files/AS_exons/H3K27me3.ks_test.svg
   :width: 45 %
Figure 1: H3K27me3


.. figure:: https://132.239.135.28/public/nucChIP/files/AS_exons/8_mnase.ks_test.svg
   :width: 45 %
Figure 5: MNase

