rm *svg *rst

i=0
for hist in mnase H3K4me3 H3K27Ac H3K9me3 H3K27me3; do
	echo ".. figure::"
	(( i++ ))
	for file in $( ls /data2/rivasas2/singleNucleosome/secondBatch/fragDistribution/*svg | grep $hist) ; do
		ln -s $file .
		file_name=$(echo $file | awk '{n=split($0,a,"/");print a[n]}')
		echo ".. image:: https://132.239.135.28/public/nucChIP/files/fragmentDist/hist/${file_name}"
		echo "   :width: 45%"
	done	
	
	echo "Figure $i: Genome-wide fragment length distribution of $hist."
	echo "" 
done

