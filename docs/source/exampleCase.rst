.. _exampleCase:
.. role:: bash(code)
   :language: bash


************
Example case
************

To exemplify the use of nucChIP, we'll use it to analyze MNChIP-seq data from Rivas-Astroza et al `(GEO accession number GSM1880581) <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?token=odsnwwionxcnbqv&acc=GSE73004>`_, and MNase data from Carone et al `(GEO accession number GSM1400766) <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1400766>`_. Overall, the procedure is divided in three steps:

1. Preprocessing: mapping of the MNChIP-seq and MNase-seq reads and data quality assesment.
2. Creation of nucleosome maps: MNChIP-seq and MNase-seq data are pooled together and use to call nucleosome positions genome-wide.
3. Tracing histone marks to individual nucleosomes: MNChIP-seq and MNase-seq reads are counted over every nucleosome they happen to overlap, and statistical methods are use to determine whether the enrichment of histone marks is strong enough to call a nuclesome as marked by a histone mark.  

In what follows, we'll describe step-by-step instructions. These instructions are based on a Linux-GNU system.

1. Preprocessing
================

To store the data we'll create the folder :bash:`data`.

.. code-block:: bash

   $mkdir data

and download into it the following data

* MNChIP-seq data containing information for four histone marks H3K4me3, H3K27Ac, H3K9me3, and H3K27me3: `GEO accession number GSM1880581 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?token=odsnwwionxcnbqv&acc=GSE73004>`_, 

* MNase data: `GEO accession number GSM1400766 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1400766>`_. 

After unziping the fastq files the :bash:`data` folder should look like this 

.. code-block:: bash

   $ gunzip * # unziping all files in the current directory
   $ ls data  # list all files in the current directory
   H3K4me3_rep1_1.fastq
   H3K4me3_rep1_2.fastq
   H3K4me3_rep2_1.fastq
   H3K4me3_rep2_2.fastq
   H3K27Ac_rep1_1.fastq
   H3K27Ac_rep1_2.fastq
   H3K27Ac_rep2_1.fastq
   H3K27Ac_rep2_2.fastq
   H3K9me3_1.fastq
   H3K9me3_2.fastq
   H3K27me3_rep1_1.fastq
   H3K27me3_rep1_2.fastq
   H3K27me3_rep2_1.fastq
   H3K27me3_rep2_2.fastq


Then, we'll map the file to the mouse genome `mm9 <ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm9.zip>`_  using `bowtie <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_, and conver the output files from SAM to BAM format using `samtools <http://www.htslib.org/>`_.

.. code-block:: bash

   # Map libraries with replicates
   for library in H3K4me3 H3K27Ac H3K27me3; do # Iterates over all fastq files
   for replicate in rep1 rep2; do
      name=${library}_${replicate} 
      # Map fastq file to mm9
      bowtie2 -x mm9 -1 ${name}_1.fastq -2 ${name}_2.fastq -p 10 -s ${name}.sam
      # Convert sam to sorted bam file
      samtools view -Sb ${name}.sam > ${name}.bam # Convert sam to bam format
      samtools sort ${name}.bam ${name}.sort      # Sort bam file
      rm ${name}.bam ${name}.sam                  # Remove unsorted bam and sam files
   done; done
   # Map H3K9me3, the only library without replicates
   name=H3K9me3 
   # Map fastq file to mm9
   bowtie2 -x mm9 -1 ${name}_1.fastq -2 ${name}_2.fastq -p 10 -s ${name}.sam
   # Convert sam to sorted bam file
   samtools view -Sb ${name}.sam > ${name}.bam # Convert sam to bam format
   samtools sort ${name}.bam ${name}.sort      # Sort bam file
   rm ${name}.bam ${name}.sam                  # Remove unsorted bam and sam files

To remove duplicates we'll use MarkDuplicates from `picard tools <http://broadinstitute.github.io/picard/>`_

.. code-block:: bash

   for file in *bam; do
      name=${$file%.bam}
      java -jar MarkDuplicates.jar \
         REMOVE_DUPLICATES = true \
         ASSUME_SORTED     = true \
         METRICS_FILE      = metricFile \
         INPUT             = $file.pairend.sort.bam \
         OUTPUT            = $file.pairend_rmdup.sort.bam
   done

Determining nucleosome positions
================================

Then we checked that the insert size of the MNChIP-seq and MNase-seq libraries were in the mono-nucleosomeal range: \~ 147 b.  For this we used the tool which generated the following figures

.. code-block:: bash

   fragDistribution \
      -b file.pairend_rmdup.sort.bam \
      -t title_in_figure \
      -bins 30 \
      -o prefix_output_file

From Figures ... it can be seen that all MNChIP-seq libraries are in the mono-nucleosomal range. The MNase-seq data, on the other hand, also showed protection of sub-nucleosomal molecules. This was the intended result of by Carone et al protocol in order to explore the positioning of a wider range of molecules. To avoid biasing effect of these reads when discovering nucleosomes we followed Carone et al procedure to discover nucleomes, namely filtering out any reads-pairs with insert sizes outside the 135-165 b range (red vertical lines on Figure ...).   


.. code-block:: bash

   bam2bed \
       -b $bam \
       -l $fragLength \
       -e 75 \
       -t bed \
       -upper 165 \
       -lower 135 \
       -o ${lib}.REG.bed
   bam2bed \
       -b $bam \
       -l $fragLength \
       -e 75 \
       -t bed \
       -o ${lib}.REG.bed

All MNChIP-seq reads plus the filtered MNase-seq reads were pooled into a single library (pooled-dataset) and used to call nucleosomes using the software iNPs.

.. code-block:: bash

   # Merge BED files by chromosomes
   for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY; do
   for lib in H3K4me3_rep1 H3K4me3_rep2 H3K27Ac_rep1 H3K27Ac_rep2 H3K9me3 H3K27me3_rep1 H3K27me3_rep2 mnase; do
      # Print only lines from chromosome chrom 
      awk -v x=$chrom '$1==x' ${lib}.REG.bed >> ${chrom}.all.bed
   done; done 

Then we run iNPs 

.. code-block:: bash

   declare -A chromSizes
   chromSizes["chr1"]=197195432
   chromSizes["chr2"]=181748087
   chromSizes["chrX"]=166650296
   chromSizes["chr3"]=159599783
   chromSizes["chr4"]=155630120
   chromSizes["chr5"]=152537259
   chromSizes["chr7"]=152524553
   chromSizes["chr6"]=149517037
   chromSizes["chr8"]=131738871
   chromSizes["chr10"]=129993255
   chromSizes["chr14"]=125194864
   chromSizes["chr9"]=124076172
   chromSizes["chr11"]=121843856
   chromSizes["chr12"]=121257530
   chromSizes["chr13"]=120284312
   chromSizes["chr15"]=103494974
   chromSizes["chr16"]=98319150
   chromSizes["chr17"]=95272651
   chromSizes["chr18"]=90772031
   chromSizes["chr19"]=61342430
   chromSizes["chrY_random"]=58682461
   chromSizes["chrY"]=15902555
   chromSizes["chrUn_random"]=5900358
   chromSizes["chrX_random"]=1785075
   chromSizes["chr1_random"]=1231697
   chromSizes["chr8_random"]=849593
   chromSizes["chr17_random"]=628739
   chromSizes["chr9_random"]=449403
   chromSizes["chr13_random"]=400311
   chromSizes["chr7_random"]=362490
   chromSizes["chr5_random"]=357350
   chromSizes["chr4_random"]=160594
   chromSizes["chr3_random"]=41899
   chromSizes["chrM"]=16299
   chromSizes["chr16_random"]=3994
   # Compute nuc locations
   nproc=0
   maxProc=7
   for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY; do
      bedChrom=${chrom}.all.bed
      libChrom=${chrom}
      echo $libChrom ==============================================================
      python3 /home/rivasas2/tools/iNPS_V1.0/iNPS_V1.0.py \
          -i $bedChrom \
          -o $libChrom \
          -c $libChrom \
          -l ${chromSizes[$libChrom]} &
      nproc=$(($nproc+1))
      if [ "$nproc" -ge "$maxProc" ]; then
         wait
         nproc=0
         echo RESET-------
      fi
   done
   #######################################################################
   # Summarize chr files into one single file
   for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY; do
      bedChrom=${chrom}.like_bed
      echo $bedChrom
      if [ "${chrom}" == "chr1" ]; then
          cat $bedChrom > allDataNuc.bed
      else # avoid header lines
          awk 'BEGIN{FS=OFS="\t"} NR>2' $bedChrom >> allDataNuc.bed
      fi
   done

This resulted in 10,292,810 nucleosomes.

Assigning histone marks to individual nucleosomes
=================================================

Once the nucleosomal map is available we'll determine whether a nucleosome is marked by a histone mark based on its enrichment of MNChIP-seq reads. We'll measure enrichment as the counts of MNChIP-seq reads overlappin a nucleosome. However, care has to be taken to avoid the co-founding effect between the co-localization level of nucleosomes and nucleosomal enrichemnt of histone marks. Namely, nucleosomes highly-colocalized have in general higher counts of MNChIP-seq reads. More details on this and to how to avoid it will be cover in the last step (Calling histone marks over nucleosomes). Let's start by counting the MNase-seq and MNase-seq reads per nucleosomes.

Count MNChIP-seq MNase-seq reads per nucleosome
-----------------------------------------------

Then, we used our genome-wide nucleosomal map to count, for each nucleosome, the overlapping MNChIP-seq reads. 

.. code-block:: bash

  getCounts \
     -b $bam \
     -n $nucFile \
     -pValue 0 \
     -l $fragLength \
     -lower $lower \
     -upper $upper \
     -e 75 \
     -o ${lib}.${nuc}.counts.bed &



Calling histone marks over nucleosomes
--------------------------------------

Expected counts were calculated with the R script `expectedCounts.R <https://szbio.ucsd.edu/public/nucChIP/files/exampleCase/expectedCounts.R>`_:

.. code-block:: bash
   
   Rscript expectedCounts.R

Determining nucleosomes marked by histone marks.

.. code-block:: bash

   getEnrichedRegions \
      -signal $signal \
      -control $control \
      -expV $expV \
      -prefix $lib

.. Bibliography
.. ============
.. 
.. .. bibliography:: Mendeley.bib
..    :style: plain
