.. countsPerNuc::

Counts per nucleosome on canonical nucleosomes around AS exons
==============================================================

We counted the number of MNase ChIP-seq reads overlapping nucleosomes. We classify exons as spliced-in or spliced-out. For each set we look for nucleosomes overlapping 2000 bp windows center on either the 5 or 3 prime end of each exon.

We required that a nucleosome should have a p-value lower or equal to 0.05.

For each window, nucleosomes were assigned to position ranging from -4 up to +5, based on the proximity to canonical positions (determined using genome-wide positioning of nucleosomes around AS exons). Figure 1, shows the canonical nucleosome positions using MNase reads and nucleosome positions computed from the 3 and 5p end of the exons. The MNase and nucleosome peaks doesn’t differ significantly, but the MNase reads produded sharper peaks, specially when computed from the 5p end of the exons. Therefore, we used the MNase peaks computed from the 5p end of the exons to call canonical nucleosome positions.



When plotting the nucleosomes coverage around AS exons (Figure 2), it can be seen that most nucleosomes are close to the exon boundaries (either 3 or 5p). This is coherent with the genome-wide trend (Figure1).



Finally, we computed the average enrichment of all MNase ChIP-seq libraries on the canonical nucleosomes. All libraries were normalized by number of total reads. Then we plotted the ratio MNase ChIP-seq over MNase (also normalized by library size). The results are presented on Figures 3-6.

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/canonicalNuc_3p_1.nuc_avrcov.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/canonicalNuc_5p_1.nuc_avrcov.svg
   :width: 45%
Figure 1: Canonical nucleosome positions

.. figure::
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/17_H3K4me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/17_H3K4me3_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/H3K4me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/H3K4me3_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/n1_H3K4me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/n1_H3K4me3_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/n2_H3K4me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/n2_H3K4me3_5p.avrRatios.svg
   :width: 45%
Figure 2: Average enrichment of H3K4me3 on the canonical nucleosomes

.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/14_H3K27Ac_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/14_H3K27Ac_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/6_H3K27Ac_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/6_H3K27Ac_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/H3K27Ac_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/H3K27Ac_5p.avrRatios.svg
   :width: 45%
Figure 3: Average enrichment of H3K27Ac on the canonical nucleosomes

.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/4_H3K9me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/4_H3K9me3_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/9_H3K9me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/9_H3K9me3_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/H3K9me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/H3K9me3_5p.avrRatios.svg
   :width: 45%
Figure 4: Average enrichment of H3K9me3 on the canonical nucleosomes

.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/12_H3K27me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/12_H3K27me3_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/5_H3K27me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/5_H3K27me3_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/H3K27me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/H3K27me3_5p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/n3_H3K27me3_3p.avrRatios.svg
   :width: 45%
.. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/n3_H3K27me3_5p.avrRatios.svg
   :width: 45%
Figure 5: Average enrichment of H3K27me3 on the canonical nucleosomes
