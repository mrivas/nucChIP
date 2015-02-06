for file in /data2/rivasas2/singleNucleosome/secondBatch/counts/AS_exons/*svg; do
	name=$(echo $file | awk '{n=split($0,a,"/");print a[n]}')
	#ln -s $file .
	echo ".. figure:: https://132.239.135.28/public/nucChIP/files/AS_exons/${name}"
	echo "   :width: 33 %"
done


