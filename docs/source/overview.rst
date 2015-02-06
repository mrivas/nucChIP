.. _overview:

********
Overview
********


Unlike regular ChIP-seq data which have a resolution of ~200 bp, MNase digested ChIP-seq produces genomic footprints at single nucleosome level. This gain in resolution allow us to observe tight relations between DNA position and histone marks such as: histone profile at TSS, and alternatively spliced exons.  

To harvest this type of data, we have developed nucChIP, a Python package designed to be easy to use by both computational and experimental biologists. 

To start analyzing your data with nucChIP follow this installation instructions

	:ref:`installation`

As tutorial, we have included detailed application case of nucChIP here

	:ref:`exampleCase`

Finally, a summary list of all nucChIP tools are presented here.

	:ref:`commandLineTools` 

To do:
   * Add nucleosome based plots of coverage at 3p and 5p ends of each exon
   * Add control plots of coverage/nucleosome-based on constitutive exons
   * Add other control plots: using publicly available libraries
