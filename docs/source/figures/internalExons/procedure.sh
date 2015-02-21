i=0
for histone in 8_mnase H3K4me3 H3K27Ac H3K9me3 H3K27me3; do
i=$(($i+1))
echo ".. figure::"
for file in $(ls /data2/rivasas2/singleNucleosome/alternativeSplicing/SE/summary_day0/summary/internal_exons/*svg | grep $histone); do
	ln -s $file .
	name=$(echo $file | awk '{n=split($0,a,"/");print a[n]}')
	echo ".. image:: https://132.239.135.28/public/nucChIP/files/internalExons/${name}"
	echo "   :width: 45%"
done
echo "Figure $i: Coverage of $histone on internal exons."
echo ""
done
