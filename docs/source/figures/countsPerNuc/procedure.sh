echo ".. figure::"
for fileName in /data2/rivasas2/singleNucleosome/secondBatch/counts/counts_per_nucleosome/distance_from_exon_overExpected/canonicalNuc_*_1.nuc_avrcov.svg; do
ln -s $fileName .
name=$(echo $fileName | awk '{n=split($0,a,"/");print a[n]}')
echo ".. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/$name"
echo "   :width: 45%"
done
echo "Figure 1: Canonical nucleosome positions"
echo ""


##echo ".. figure::"
##for fileName in /data2/rivasas2/singleNucleosome/alternativeSplicing/SE/summary_day0/summary/ASexonCov/nucleosomesHighQual_exon*DB.avrcov.svg; do
###ln -s $file .
##name=$(echo $fileName | awk '{n=split($0,a,"/");print a[n]}')
##echo ".. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/$name"
##echo "   :width: 45%"
##done
##echo "Figure 2: Nucleosomes coverage around AS exons"
##echo ""

i=2
echo ".. figure::"

for histone in H3K4me3 H3K27Ac H3K9me3 H3K27me3; do
for fileName in $(ls /data2/rivasas2/singleNucleosome/secondBatch/counts/counts_per_nucleosome/ref_nuc_pos_overExpected/*.svg | grep $histone );do
ln -s $fileName .
name=$(echo $fileName | awk '{n=split($0,a,"/");print a[n]}')
echo ".. image:: https://132.239.135.28/public/nucChIP/files/countsPerNuc/$name"
echo "   :width: 45%"
done
echo "Figure $i: Average enrichment of $histone on the canonical nucleosomes"
echo ""
i=$(($i+1))
done

