i=0
for histone in 8_mnase H3K4me3 H3K27Ac H3K9me3 H3K27me3; do
i=$(($i+1))
echo ".. figure::"
for file in $(ls /data2/rivasas2/singleNucleosome/secondBatch/counts/counts_per_nucleosome/distance_from_exon_overExpected/tss/*svg | grep $histone); do
	ln -f -s $file .
	name=$(echo $file | awk '{n=split($0,a,"/");print a[n]}')
	echo ".. image:: https://132.239.135.28/public/nucChIP/files/coverage/tss/${name}"
	echo "   :width: 45%"
done
echo "Figure $i: Normalized coverage of $histone on TSS."
echo ""
done
