.. _exampleCase:

************
Example case
************

Here we present the use of nucChIP to analyze our in-house data. We generated MNased digested ChIP-seq libraries of several mouse histone marks as well as MNase data alone. The samples were taken from mouse E14 and mouse E14-derived cells (day 4 after fibroblast inducing factors were added to the medium).


1. Data summary
===============

We generated single-nucleosome ChIP-seq libraries of four histone marks (H4K4me3, H3K27Ac, H3K9me3, and H3K27me3), and MNase (See Table 1). In our ChIP-seq protocol, we used MNase digestion rather than sonication to produce genomic histone footprints at single-nucleosome resolution. The libraries were extracted from two cell-lines: mouse E14, and mouse fibroblast-like (See Methods).


`Data summary of mouse E14 libraries <https://docs.google.com/spreadsheet/ccc?key=0Aueh7dagaPEZdENBUUR1Qk8tS3hhbnZFZ2NyU29CbEE#gid=4>`_

`Data summary of mouse fibroblast-like libraries <https://docs.google.com/spreadsheet/ccc?key=0Aueh7dagaPEZdENBUUR1Qk8tS3hhbnZFZ2NyU29CbEE#gid=4>`_

but for the sake of brevity, you can see a quick description of on the following tables.

.. csv-table:: Mouse E14 libraries.
   :header: "Library type","Experiment","# Reads","# Fragments", "Replicate"

   MNased_ChIP-seq, H3K27me3,41727332,3390346,`12_H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/12_H3K27me3.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq, H3K27me3,84379845,35072483,`5_H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/5_H3K27me3.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq,H3K27me3,16024783,7572511,`n3_H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/n3_H3K27me3.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq, H3K9me3,91197916,38617758,`4_H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/4_H3K9me3.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq, H3K9me3,63153628,23932067,`9_H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/9_H3K9me3.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq, H3K27Ac,72746990,23886474,`14_H3K27Ac <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/14_H3K27Ac.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq, H3K27Ac,36822734,13953975,`6_H3K27Ac <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/6_H3K27Ac.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq, H3K4me3,16206340,796542,`17_H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/17_H3K4me3.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq,H3K4me3,23451320,11277740,`n1_H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/n1_H3K4me3.pairend_rmdup.sort.bam>`_
   MNased_ChIP-seq,H3K4me3,14765216,6819315,`n2_H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/n2_H3K4me3.pairend_rmdup.sort.bam>`_
   MNase,MNase,203367131,120922096,`all_nuc <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/all_nuc.singlend_rmdup.sort.bam>`_
   MNase,MNase,203367131,120922096,`8_mnase <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/8_mnase.sort_rmdup.bam>`_
   Sonicated_ChIP-seq,H3K27Ac,18739298,13312397,`H3K27Ac <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/H3K27Ac.singend_rmdup.sort.bam>`_
   Sonicated_ChIP-seq,H3K27me3,13680637,9598335,`H3K27me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/H3K27me3.singend_rmdup.sort.bam>`_
   Sonicated_ChIP-seq,H3K4me3,6687568,4275362,`H3K4me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/H3K4me3.singend_rmdup.sort.bam>`_
   Sonicated_ChIP-seq,H3K9me3,18851693,4039918,`H3K9me3 <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/H3K9me3.singend_rmdup.sort.bam>`_
   Sonicated_ChIP-seq,IgG,3529661,2143563,`IgG <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/IgG.singend_rmdup.sort.bam>`_


.. csv-table:: Mouse E14 derived libraries: day 4 after differentiation induction.
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

.. note:: The library types **MNase_ChIP-seq**, **MNase**, and **Sonicated_ChIP-seq** correspond to MNase digested, MNase alone, and Sonicated ChIP-seq libraries, respectively. All but Library **8_mnase** of day 0 are in-housed data. Library **8_mnase** was downloaded from:
   http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1400766

Reproducibility
===============

We assessed the information content of our single-nucleosome histone libraries by measuring their similarity with referencial sonicated ChIP-seq libraries to check if they cluster together.

First, we had to preprocess all our libraries. Although MNase digested ChIP-seq produces fragment lengths of 147 nu- cleotides (nt) –the length of DNA wrapped around a nucleosome– in practice fragments sizes may be larger or smaller. To improve the resolution of our data, we regularized the fragment lengths by extending +/-75 bases around each paired-end read middle point (see Figure 1). In the case of single-end libraries we enlarged each read up to 200 bases towards the 3’ end, and extended their length +/- 75 bases around their resulting middle points. Thus, all reads were regularized to a fragment lengths (150 bases) similar to the nucleosomal DNA size.  

.. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/shift.svg

   Figure1: Regularization of reads.

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

For single-nucleosomal and reference ChIP-seq libraries, coverage was computed on 50-bases windows at chromosome 1 and normalized by library size. For all histone marks, we used `Carone <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1400766>_`\’s MNase data (8_mnase) as background to compute the signal over background ratios (pseudo-values of 1 were added to both numerator and denominators). 

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

Regardless of the method used, we found that H3K27Ac and H3K27me3 clustered among their single-nucleosomal replicates, respectively (see Figures 2 and 3). Conversely, H3K4me3 and H3K9me3 single-nucleosomal replicates were largely dispersed on the resulting dendrograms. This is due to differences on MNase digestion intensity (see next section) and library sizes. In all cases, however, single-nucleosomal libraries didn’t cluster closed to their corresponding reference ChIP-seq libraries. We hypothesized that this is due to the resolution differences between both types of datasets.


.. _fig_clusterCoverage:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/dendrogram_ratio_coverage.svg
   
   Figure 2. Cluster of libraries using counts per nucleosome.

.. _fig_clusterCounts:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/dendrogram_ratio_counts.svg
   
   Figure 3. Cluster of libraries using reads' coverage.

.. _fig_H3K4me3_TSS:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n1_H3K4me3.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n2_H3K4me3.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/17_H3K4me3.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n1_H3K4me3_on_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n2_H3K4me3_on_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/17_H3K4me3_on_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n1_H3K4me3_off_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n2_H3K4me3_off_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/17_H3K4me3_off_genes_1500.vplot.svg
      :width: 30 %
   
   Figure 4. Coverage and fragment lendth of H3K4me3 reads at TSS.

.. _fig_H3K27Ac_TSS:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/6_H3K27Ac.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/14_H3K27Ac.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/6_H3K27Ac_on_genes_1500.vplot.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/14_H3K27Ac_on_genes_1500.vplot.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/6_H3K27Ac_off_genes_1500.vplot.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/14_H3K27Ac_off_genes_1500.vplot.svg
      :width: 40 %
   
   Figure 5. Coverage and fragment lendth of H3K27Ac reads at TSS.

.. _fig_H3K9me3_TSS:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/4_H3K9me3.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/9_H3K9me3.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/4_H3K9me3_on_genes_1500.vplot.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/9_H3K9me3_on_genes_1500.vplot.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/4_H3K9me3_off_genes_1500.vplot.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/9_H3K9me3_off_genes_1500.vplot.svg
      :width: 40 %
   
   Figure 6. Coverage and fragment lendth of H3K9me3 reads at TSS.
   
.. _fig_H3K27me3_TSS:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/5_H3K27me3.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/12_H3K27me3.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n3_H3K27me3.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/5_H3K27me3_on_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/12_H3K27me3_on_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n3_H3K27me3_on_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/5_H3K27me3_off_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/12_H3K27me3_off_genes_1500.vplot.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n3_H3K27me3_off_genes_1500.vplot.svg
      :width: 30 %
   
   Figure 7. Coverage and fragment lendth of H3K27me3 reads at TSS.
   
.. _fig_8_mnase_TSS:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/8_mnase.nineTiles_rmdup_5perc_1500.avrcov.svg
      :width: 51 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/8_mnase_on_genes.3tercile.vplot.svg
      :width: 51 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/8_mnase_off_genes.vplot.svg
      :width: 51 %
   
   Figure 8. Coverage and fragment lendth of MNase reads at TSS.
   
.. _fig_H3K4me3_lm:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n1_H3K4me3.lm.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n2_H3K4me3.lm.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/17_H3K4me3.lm.svg
      :width: 30 %
   
   Figure 9. Correlation between gene expression and nucleosomes-specific enrichment of H3K4me3.

.. _fig_H3K27Ac_lm:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/6_H3K27Ac.lm.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/14_H3K27Ac.lm.svg
      :width: 40 %
   
   Figure 10. Correlation between gene expression and nucleosomes-specific enrichment of H3K27Ac.

.. _fig_H3K9me3_lm:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/4_H3K9me3.lm.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/9_H3K9me3.lm.svg
      :width: 40 %
   
   Figure 11. Correlation between gene expression and nucleosomes-specific enrichment of H3K9me3.

.. _fig_H3K27me3_lm:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/5_H3K27me3.lm.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/12_H3K27me3.lm.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n3_H3K27me3.lm.svg
      :width: 30 %
   
   Figure 12. Correlation between gene expression and nucleosomes-specific enrichment of H3K27me3.

.. _fig_H3K4me3_as:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n1_H3K4me3_exonStartDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n2_H3K4me3_exonStartDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/17_H3K4me3_exonStartDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n1_H3K4me3_exonEndDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n2_H3K4me3_exonEndDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/17_H3K4me3_exonEndDB.avrcov.svg
      :width: 30 %
   
   Figure 13. Enrichement of H3K4me3 on alternatively spliced exons.

.. _fig_H3K27Ac_as:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/6_H3K27Ac_exonStartDB.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/14_H3K27Ac_exonStartDB.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/6_H3K27Ac_exonEndDB.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/14_H3K27Ac_exonEndDB.avrcov.svg
      :width: 40 %
   
   Figure 14. Enrichement of H3K27Ac on alternatively spliced exons.

.. _fig_H3K9me3_as:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/4_H3K9me3_exonStartDB.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/9_H3K9me3_exonStartDB.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/4_H3K9me3_exonEndDB.avrcov.svg
      :width: 40 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/9_H3K9me3_exonEndDB.avrcov.svg
      :width: 40 %
   
   Figure 15. Enrichement of H3K9me3 on alternatively spliced exons.

.. _fig_H3K27me3_as:

   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/5_H3K27me3_exonStartDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/12_H3K27me3_exonStartDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n3_H3K27me3_exonStartDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/5_H3K27me3_exonEndDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/12_H3K27me3_exonEndDB.avrcov.svg
      :width: 30 %
   .. image:: https://132.239.135.28/public/nucChIP/files/exampleCase/n3_H3K27me3_exonEndDB.avrcov.svg
      :width: 30 %
   
   Figure 16. Enrichement of H3K27me3 on alternatively spliced exons.









3. Call nucleosome positions from Mnase data
============================================

The first step is to use Mnase-seq data to determine the nucleosome positions. We used `iNPS <http://www.picb.ac.cn/hanlab/iNPS.html>`_, to produce . 


DANPOS will create a folder called danpos_mnaseData where there is a file with extension ``xls`` containing the nucleosome positions along their statistics. We'll use this on section 4, when we'll count the number of Mnased ChIP-seq reads per nucleosome. 

3. Plot coverage and heatmaps 
=============================

Explain quantiles distribution

4. Count ChIP-seq reads per nucleosome
======================================

