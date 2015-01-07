.. _exampleCase:

************
Example case
************

Here we present the use of nucChIP to analyze our in-house data. We generated MNased digested ChIP-seq libraries of several mouse histone marks as well as MNase data alone. The samples were taken from mouse E14 and mouse fibroblast-like cells (day 4 after guided differentiation of E14 cells fibroblast).


1. Data summary
===============

We generated single-nucleosome ChIP-seq libraries of four histone marks (H4K4me3, H3K27Ac, H3K9me3, and H3K27me3), and used a publicly available MNase data set (un-spunned data from :cite:`Carone2014`). A summary of all data is presented in Table 1. In our ChIP-seq protocol, we used MNase digestion rather than sonication to produce genomic histone footprints at single-nucleosome resolution. The libraries were extracted from two cell-lines: mouse E14, and mouse fibroblast-like (See Methods).

..
   Data source:
   `Data summary of mouse E14 libraries <https://docs.google.com/spreadsheet/ccc?key=0Aueh7dagaPEZdENBUUR1Qk8tS3hhbnZFZ2NyU29CbEE#gid=4>`_
   `Data summary of mouse fibroblast-like libraries <https://docs.google.com/spreadsheet/ccc?key=0Aueh7dagaPEZdENBUUR1Qk8tS3hhbnZFZ2NyU29CbEE#gid=4>`_

.. csv-table:: Table 1: Mouse E14 libraries.
   :header: "Library (type)","Replicate ID","# Reads","Alignment rate (%)","# Fragments", "# Fragment MAPQ>=20"

   H3K27me3 (MNase),`12_H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/12_H3K27me3.pairend_rmdup.sort.bam>`_,41727332,16.17,3390346,2519561
   H3K27me3 (MNase),`5_H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/5_H3K27me3.pairend_rmdup.sort.bam>`_,84379845,82.77,35072483,30405555
   H3K9me3 (MNase),`4_H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/4_H3K9me3.pairend_rmdup.sort.bam>`_,91197916,84.16,38617758,33180391
   H3K9me3 (MNase),`9_H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/9_H3K9me3.pairend_rmdup.sort.bam>`_,63153628,75.31,23932067,19637520
   H3K27Ac (MNase),`14_H3K27Ac <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/14_H3K27Ac.pairend_rmdup.sort.bam>`_,72746990,65.41,23886474,20931875
   H3K27Ac (MNase),`6_H3K27Ac <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/6_H3K27Ac.pairend_rmdup.sort.bam>`_,36822734,75.74,13953975,11805574
   H3K4me3 (MNase),`17_H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/17_H3K4me3.pairend_rmdup.sort.bam>`_,16206340,10.72,796542,401970
   H3K4me3 (MNase),`n1_H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/n1_H3K4me3.pairend_rmdup.sort.bam>`_,23451320,95.5,11277740,9703381
   H3K4me3 (MNase),`n2_H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/n2_H3K4me3.pairend_rmdup.sort.bam>`_,14765216,94.03,6819315,5267592
   H3K27me3 (MNase),`n3_H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/n3_H3K27me3.pairend_rmdup.sort.bam>`_,16024783,95.44,7572511,6084379
   MNase (MNase),`8_mnase <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/8_mnase.sort_rmdup.bam>`_ ,203367131,96.06,120922096,44912058*
   H3K27Ac (Sonicated),`H3K27Ac <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/H3K27Ac.singend_rmdup.sort.bam>`_,18739298,71.04,13312397,NA
   H3K27me3 (Sonicated),`H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/H3K27me3.singend_rmdup.sort.bam>`_,13680637,70.16,9598335,NA
   H3K4me3 (Sonicated),`H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/H3K4me3.singend_rmdup.sort.bam>`_,6687568,63.93,4275362,NA
   H3K9me3 (Sonicated),`H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/H3K9me3.singend_rmdup.sort.bam>`_,18851693,21.43,4039918,NA

.. 
   MNase data to be added to the table
   Mnase (MNase),`all_nuc <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/all_nuc.singlend_rmdup.sort.bam>`_,203367131,59.46,120928017,120928017
   IgG data to be added to previous table  
   Sonicated_ChIP-seq,IgG,3529661,2143563,`IgG <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/IgG.singend_rmdup.sort.bam>`_

..
   .. csv-table:: Table 2: Mouse fibroblast-like libraries.
   :header: "Library type","Experiment","# Reads","# Fragments", "Replicate"
   MNased_ChIP-seq,H3K4me3,57587756,25663983,`4_H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/154-176>`_
   MNased_ChIP-seq,H3K9me3,17306377,2662586,`51_H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/152-169>`_
   MNased_ChIP-seq,H3K9me3,49969650,23075984,`52_H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/146-163>`_
   MNased_ChIP-seq,H3K27me3,28455820,5869013,`6_H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/157-172>`_
   MNase,MNase,33290631,13757353,`1_mnase <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/156-176>`_
   MNase,MNase,33464926,13270516,`2_mnase <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/156-176>`_
   MNase,MNase,57341704,25740691,`3_mnase <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/153-177>`_
   MNase,MNase,22997958,7509983,`4_mnase <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/158-176>`_
   Sonicated_ChIP-seq,H3K27me3,13638545,9242842,`H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/NA>`_
   Sonicated_ChIP-seq,H3K4me3,23088459,15974905,`H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/NA>`_
   Sonicated_ChIP-seq,H3K9me3,23839442,9764635,`H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/NA>`_
   Sonicated_ChIP-seq,IgG,16552070,11278580,`IgG <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/NA>`_

.. note:: The library types **MNase** and **Sonicated**, correspond to MNase digested and Sonicated ChIP-seq libraries, respectively. \* Reads from **8_mnase** library were also filter to represent only fragment lengths in the nucleosomal range: 135 to 165 (nt).

Nucleosome positions
====================

We used iNPS :cite:`Chen2014` to call nucleosome locations with default parameters. To improve the signal over background ratio, iNPS reduces each fragment length to 75 (nt) around their midpoint. Using 8_mnase fragments with mapping qualities above 20 as input, we found 10,468,598 nucleosome locations genome wide. iNPs classified these 4 kinds. The first, MainPeak representing isolated main nucleosome peaks, accounted for the majority of the calls (80% of the nucleosomes). This was accompanied by the categories MainPeak+Shoulder (11%), MainPeak:doublet (6%), and Shoulder (3%). Compared to MainPeak, Shoulder nucleosomes are closer to flanking nucleosomes and have shorted 'widths'. On the contrary, MainPeak+Shoulder and MainPeak:doublet nucleosomes are farther apart form nearby nucleosoms and with larger 'widths'. Unlike MainPeak, the last three categories suggest association with unstable nucleosomes positions :cite:`Chen2014`. 

All together, this amount nucleosomes positions was expected given the size of the mouse genome. The total number of nucleosomes times the combined length of each nucleosomal DNA (147 nt) and its linker sequence (38 nt as the typical distance between neighbors nucleosomes; :cite:`Jiang2009` ) covered approximately 77% of the mouse genome length (2.5 Gb; :cite:`Waterston2002`).

As shown in Figure :num:`#fig-width`, the nucleosome's widths peaks at ~75 nt, which is coherent the lenght used by iNPS to represent the enrichment signals. The sharp peaks is signal of both: well positioned nucleosomes, and isolated nucleosomes. On the other hand, the distance between adjacent nucleosomes has a tipical value of ~ 180 (nt), being this coherent with the nucleosomal DNA length (~147 nt) and the linker DNA length (~38 nt; :cite:`Jiang2009`)

.. _fig-width:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/hist_width.svg
   :width: 50%

   Distribution of nucleosomes' 'widths'.

.. _fig-dist:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/hist_dist.svg
   :width: 50%

   Distribution of distances between adjacent nucleosomes.


Reproducibility
===============

We assessed the information content of our single-nucleosome histone libraries by measuring their similarity with referencial sonicated ChIP-seq libraries to check if they cluster together.

First, we had to preprocess all our libraries. Although MNase digested ChIP-seq produces fragment lengths of 147 nu- cleotides (nt) –the length of DNA wrapped around a nucleosome– in practice fragments sizes may be larger or smaller. To improve the resolution of our data, we regularized the fragment lengths by extending +/-75 bases around each paired-end read middle point (see :num:`#fig-shift`). In the case of single-end libraries we enlarged each read up to 200 bases towards the 3’ end, and extended their length +/- 75 bases around their resulting middle points. Thus, all reads were regularized to a fragment lengths (150 bases) similar to the nucleosomal DNA size.  

.. _fig-shift:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/shift.svg

   Regularization of fragments length. (A) After MNase digestion, paired or single-end reads are produced. (B) Fragments lengths are regularized to a custom length by extending their middle point +/- 75 bases. Unlike paired-end reads which provides the fragment length, single-end fragment lengths are estimated as the average fragment length of the library. (C) Provided the nucleosome locations, reads are counted per each nucleosome, and an enrichment ratio is computed.

We did this transformation with the nucChIP's script :class:`bam2bed`.

.. program-output:: bam2bed -h

For the regularization of mouse E14, three thing are worth noticing. First, library **8_mnase** contain reads spanning lengths from 0 to 200 bp, therefore we restricted the conversion only to those reads with reads in the nucleosomal range [135, 165] bp. Second, library **all_nuc** is single-end, accordingly we assumed that the typical fragment length was 150 bp. The set of these parameters can be seen on the ``if ... else`` statements. Finally, by default :class:`bam2bed` uses only reads with mapping qualities of 20 or more, filtering out --in this case-- all pairs where at least one reads has lower qualities.

.. code-block:: bash
   
   ############################################################################################
   # Convert BAM files into regularized BED files
   nproc=0
   maxProc=10
   for bam in $(ls /data2/rivasas2/singleNucleosome/secondBatch/rmdup/byReplicates/*_rmdup.sort.bam /data2/rivasas2/singleNucleosome/Teif_data/alignments/8_mnase.sort_rmdup.bam); do
      lib=$(echo $bam | awk '{n=split($0,a,"/");split(a[n],b,".");print b[1]}' )
   
       if [ "$lib" == "8_mnase" ]; then
          fragLength=0
          upper=165
          lower=135
       elif [ "$lib" == "all_nuc" ]; then
          fragLength=150
          upper=200
          lower=0
      else
          fragLength=0
          upper=200
          lower=0
       fi
   
      echo "Converting to BED the lib frag lower upper " $lib $fragLength $lower $upper
      bam2bed \
          -b $bam \
          -l $fragLength \
          -e 75 \
          -t bed \
          -upper $upper \
          -lower $lower \
          -o ${lib}.REG.bed & 
      nproc=$(($nproc+1))
   
      if [ "$nproc" -ge "$maxProc" ]; then
          wait
          nproc=0
          echo RESET-------
      fi
          
   done

For single-nucleosomal and reference ChIP-seq libraries, coverage was computed on 50-bases windows at chromosome 1 and normalized by library size. For all histone marks, we used :cite:`Carone2014` MNase data (8_mnase) as background to compute the signal over background ratios (pseudo-values of 1 were added to both numerator and denominators). 

.. code-block:: bash

   echo "######################################################################"
   echo "Creates windows of 50 bases on chrom 1"
   fetchChromSizes mm9 | grep -w 'chr1' > mm9.chromSizes
   bedtools makewindows -g mm9.chromSizes -w 50 > windows50.bed
   
   echo "#####################################################################"
   echo "Compute coverage on each window"
   for bed in /data2/rivasas2/singleNucleosome/secondBatch/nucLocation/*bed /data2/rivasas2/singleNucleosome/secondBatch/nucLocation/regularChIP/*bed; do
      lib=$(echo $bed | awk '{n=split($0,a,"/");split(a[n],b,".");print b[1]}' )
      libSize=$(wc -l $bed | awk '{print $1}')
      scale=$(echo $libSize | awk '{print 1000000/$0}')
      echo "Computing coverage for lib libSize scale" $lib $libSize $scale
      coverageBed -a $bed -b windows50.bed | awk -v x=${scale} 'BEGIN{FS=OFS="\t"}{scaled=$4*x; print $0,scaled}' > ${lib}.coverage.bed
   #done
   
   for bed in *bed; do
       lib=${bed%%.coverage.bed}
       echo "Sorting" $lib
       sort -k1,1 -k2,2n ${lib}.coverage.bed > ${lib}.coverage_sorted.bed
   done

Alternatively, we also used the MNase data to determined the position of 10M nucleosomes genome-wide (See Methods). Over the genetic span of high quality nuclesome regions (iNPs, p-value<= 0.05) on chromosome 1, we counted the number of overlapping histone fragments, and normalized them as fragments per kilobase of nucleosome per million fragment mapped (FPKM). 

.. code-block:: bash

   #######################################################################
   # MNase digested ChIP
   nproc=0
   maxProc=6
   for bam in /data2/rivasas2/singleNucleosome/secondBatch/rmdup/byReplicates/*_rmdup.sort.bam /data2/rivasas2/singleNucleosome/Teif_data/alignments/8_mnase.sort_rmdup.bam; do
      lib=$(echo $bam | awk '{n=split($0,a,"/");split(a[n],b,".");print b[1]}' )
   
      if [ "$lib" == "all_nuc" ]; then
          fragLength=150
          upper=200
          lower=0
      elif [ "$lib" == "8_mnase" ]; then
          fragLength=0
          upper=165
          lower=135
      else
          fragLength=0
          upper=200
          lower=0
      fi
      for nucFile in /data2/rivasas2/singleNucleosome/secondBatch/nucLocation/8_mnase/8_mnase.bed /data2/rivasas2/singleNucleosome/secondBatch/nucLocation/all_nuc/all_nuc.bed; do
          nuc=$(echo $nucFile | awk '{n=split($0,a,"/");split(a[n],b,".");print b[1]}' )
          echo "Computing count for " $lib $fragLength $lower $upper $nuc
          getCounts \
              -b $bam \
              -n $nucFile \
              -pValue 1.3 \
              -l $fragLength \
              -lower $lower \
              -upper $upper \
              -e 75 \
              -o ${lib}.${nuc}.counts.bed &
          nproc=$(($nproc+1))
          
          if [ "$nproc" -ge "$maxProc" ]; then
              wait
              nproc=0
              echo RESET-------
          fi
      done
      
   done
   
   ##################################################################
   # Regular ChIP
   nproc=0
   maxProc=4
   for bam in $(ls /data2/rivasas2/singleNucleosome/secondBatch/rmdup/regularChIP/*.singend_rmdup.sort.bam | grep -v IgG ); do
      lib=$(echo $bam | awk '{n=split($0,a,"/");split(a[n],b,".");print b[1]}' )
   
      fragLength=200
      upper=200
      lower=0
      for nucFile in /data2/rivasas2/singleNucleosome/secondBatch/nucLocation/8_mnase/8_mnase.bed /data2/rivasas2/singleNucleosome/secondBatch/nucLocation/all_nuc/all_nuc.bed; do
          nuc=$(echo $nucFile | awk '{n=split($0,a,"/");split(a[n],b,".");print b[1]}' )
          echo "Computing count for " $lib $fragLength $lower $upper $nuc
          getCounts \
              -b $bam \
              -n $nucFile \
              -pValue 1.3 \
              -l $fragLength \
              -lower $lower \
              -upper $upper \
              -e 75 \
              -o ${lib}.${nuc}.counts.bed &
          nproc=$(($nproc+1))
          
          if [ "$nproc" -ge "$maxProc" ]; then
              wait
              nproc=0
              echo RESET-------
          fi
      done
      
   done

For both cases, we used two metrics to measure libraries similarities: Spearman correlation, and Euclidean distance. The clustering was done with these R's scripts:

R_script for coverage
R_script for counts per nucleosome

Regardless of the method used, we found that H3K27Ac and H3K27me3 clustered among their single-nucleosomal replicates, respectively (see Figures :num:`#fig-clustercoverage` and :num:`#fig-clustercounts`). Conversely, H3K4me3 and H3K9me3 single-nucleosomal replicates were largely dispersed on the resulting dendrograms. This is due to differences on MNase digestion intensity (see next section) and library sizes. In all cases, however, single-nucleosomal libraries didn’t cluster closed to their corresponding reference ChIP-seq libraries. We hypothesized that this is due to the resolution differences between both types of datasets.


.. _fig-clusterCoverage:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/dendrogram_ratio_coverage.svg
   
   Cluster of libraries using reads' coverage.

.. _fig-clusterCounts:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/dendrogram_ratio_counts.svg
   
   Cluster of libraries using counts per nucleosome.


H3K4me3 and H3K27Ac footprints are coherent with known epigenetic functions
===========================================================================

Then, we asked: are the single-nucleosome ChIP-seq libraries coherent with the known biological associations of regular ChIP-seq? H3K4me3 and H3K27Ac are known to be preferentially enriched on the transcription starting site (TSS) of active genes. On the contrary, H3K9me3, and H3K27me3 are preferentially enriched at inactive genes' TSS. Concordantly, we classify genes according to their expression levels. For each expression level we computed the single-nucleosomal histone enrichment (see Figures :num:`#fig-h3k4me3-tss`, :num:`#fig-h3k27ac-tss`, :num:`#fig-h3k9me3-tss`, and :num:`#fig-h3k27me3-tss`).

H3K4me3 and H3K27Ac footprints are coherent with literature :cite:`Carone2014`. They were depleted immediately before TSS, but enriched thereafter in a manner proportional to gene expression. Conversely, among H3K9me3 replicates only one (9\_H3K9me3) was inversely correlated with gene expression whereas all H3K27me3's replicates were positively correlated with gene expression.

This may be the result of TSS specific biases. According to REF, regions flanked by nucleosome depleted regions are easier to digest by MNase and in consequence their presence may be overrepresented on the ChIP-seq library. This problem is significant if the MNase digestion is poor. To explore this possibility, we plotted the fragment densities around the TSS. Fragment densities about the nucleosomal size (150 nt) mean poor MNase digestion. Figures see Figures :num:`#fig-h3k4me3-tss`, :num:`#fig-h3k27ac-tss`, :num:`#fig-h3k9me3-tss`, and :num:`#fig-h3k27me3-tss` show that this is the case for all but n1\_H3K4me3.  

.. _fig-h3k4me3-tss:
   
.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K4me3_tss.svg
   :width: 90 %   
   
   Coverage and fragment lendth of H3K4me3 reads at TSS.

.. _fig-h3k27ac-tss:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K27Ac_tss.svg
   :width: 60 %   
   
   Coverage and fragment lendth of H3K27Ac reads at TSS.

.. _fig-h3k9me3-tss:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K9me3_tss.svg
   :width: 60 %   
   
   Coverage and fragment lendth of H3K9me3 reads at TSS.
   
.. _fig-h3k27me3-tss:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K27me3_tss.svg
   :width: 90 %   
   
   Coverage and fragment lendth of H3K27me3 reads at TSS.
   
.. _fig-8-mnase-tss:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/8_mnase_tss.svg
   :width: 30 %   

   Coverage and fragment lendth of MNase reads at TSS.
   

Position-specific nucleosomes control gene expression
=====================================================

The positive correlation between gene expression and the enrichment levels of H3K4me3 at the +1 nucleosome, leaded us to hypothesized that gene expression is controlled by position-specific nucleosome epigenetic marks. We tested this hypothesis, by computing the linear correlation between the average histone enrichment at eight gene expression quantiles (inactive genes were not included). 

Among H3K4me3 replicates, enrichment on +1 nucleosome positively correlates with gene expression (:num:`fig-h3k4me3-lm`) in a consistent manner.  In the case of H3K27Ac, both nucleosome +1 and +2 positively correlates with gene expression (:num:`fig-h3k27ac-lm`). Surprisingly, H3K9me3 and H3K27me3 are also positively correlated to gene expression at nucleosomes +1 and +2 which is the opposite as expected (Figures :num:`#fig-h3k9me3-lm` and :num:`#fig-h3k27me3-lm`). As discussed before, we believe that this is an artifact due to MNase under-digestion.

.. _fig-H3K4me3-lm:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K4me3_lm.svg
   :width: 90 %
   
   Correlation between gene expression and nucleosomes-specific enrichment of H3K4me3.

.. _fig-H3K27Ac-lm:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K27Ac_lm.svg
   :width: 60 %
   
   Correlation between gene expression and nucleosomes-specific enrichment of H3K27Ac.

.. _fig-H3K9me3-lm:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K9me3_lm.svg
   :width: 60 %
   
   Correlation between gene expression and nucleosomes-specific enrichment of H3K9me3.

.. _fig-H3K27me3-lm:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K27me3_lm.svg
   :width: 90 %
   
   Correlation between gene expression and nucleosomes-specific enrichment of H3K27me3.

Inclusion of alternatively spliced exons is signaled by histone marks
=====================================================================

Alternative splicing is know to be affected by the epigenome [REF], yet current ChIP-seq resolutions have hampered our ability to understand their effect. We used our single-nucleosomal ChIP-seq libraries to compute the histone enrichment at alternatively spliced exons. We used mouse E14 RNA-seq data classify alternatively spliced exons as included or excluded on the the mRNAs. On each set we then computed the histones' enrichment as signal over background rations. Again we used MNase as a background.  The results (Figures :num:`#fig-h3k4me3-as`, :num:`#fig-h3k27ac-as`, :num:`#fig-h3k9me3-as`, and :num:`#fig-h3k27me3-as` show preferential enrichment of all histone marks on included exons. 

.. _fig-H3K4me3-as-ratios:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K4me3_as_ratios.svg
   :width: 90 %
   
   Normalized enrichement of H3K4me3 on alternatively spliced exons.

.. _fig-H3K4me3-as:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K4me3_as.svg
   :width: 90 %
   
   Enrichement of H3K4me3 on alternatively spliced exons.

.. _fig-H3K27Ac-as-ratios:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K27Ac_as_ratios.svg
   :width: 60 %
   
   Normalized enrichement of H3K27Ac on alternatively spliced exons.

.. _fig-H3K27Ac-as:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K27Ac_as.svg
   :width: 60 %
   
   Enrichement of H3K27Ac on alternatively spliced exons.

.. _fig-H3K9me3-as-ratios:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K9me3_as_ratios.svg
   :width: 60 %
   
   Normalized enrichement of H3K9me3 on alternatively spliced exons.

.. _fig-H3K9me3-as:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K9me3_as.svg
   :width: 60 %
   
   Enrichement of H3K9me3 on alternatively spliced exons.

.. _fig-H3K27me3-as-ratios:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K27me3_as_ratios.svg
   :width: 90 %
   
   Normalzied enrichement of H3K27me3 on alternatively spliced exons.

.. _fig-H3K27me3-as:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/H3K27me3_as.svg
   :width: 90 %

.. _fig-8-mnase-as:

.. figure:: https://132.239.135.28/public/nucChIP/files/exampleCase/8_mnase_as.svg
   :width: 30 %
   
   Enrichement of MNase on alternatively spliced exons.

Bibliography
============

.. bibliography:: Mendeley.bib
   :style: plain
