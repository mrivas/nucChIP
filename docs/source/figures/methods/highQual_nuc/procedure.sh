i=0
for histone in H3K4me3 H3K27Ac H3K9me3 H3K27me3; do
i=$(($i+1))
echo ".. figure::"
for file in $(ls /data2/rivasas2/singleNucleosome/secondBatch/counts/expected_counts/highQual_nuc/*svg | grep $histone); do
	ln -s $file .
	name=$(echo $file | awk '{n=split($0,a,"/");print a[n]}')
	echo ".. image:: https://132.239.135.28/public/nucChIP/files/methods/highQual_nuc/${name}"
	echo "   :width: 45%"
done
echo "Figure $i: Expected number of $histone reads given supporting MNase reads per nucleosome."
echo ""
done
