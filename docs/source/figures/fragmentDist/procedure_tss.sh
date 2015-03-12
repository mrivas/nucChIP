rm *svg *rst

for hist in mnase H3K4me3 H3K27Ac H3K9me3 H3K27me3; do

for file in $( ls /data2/rivasas2/singleNucleosome/secondBatch/v-plots/tss/*svg | grep $hist) ; do
	
	ln -s $file .
	file_name=$(echo $file | awk '{n=split($0,a,"/");print a[n]}')
	name=$(echo $file | awk '{n=split($0,a,"/");split(a[n],b,"_");print b[1]"_"b[2]}')
	
	if [ -f ${hist}_fragmentDist_tss.rst ]; then
		nLine=$(wc -l ${hist}_fragmentDist_tss.rst | awk '{print $1}')
		mod=$(($nLine%7))
		nFig=$(( $nLine/7 + 1 )) # every 3 lines (constitutive, spliced_in, spliced_out) add a new figure
	else
		mod=0
		nFig=0
	fi

	if [ "$mod" == 0 ]; then 
		echo ".. figure::" >> ${hist}_fragmentDist_tss.rst
	fi

	echo ".. image:: https://132.239.135.28/public/nucChIP/files/fragmentDist/tss/${file_name}" >> ${hist}_fragmentDist_tss.rst
	echo "   :width: 45%" >> ${hist}_fragmentDist_tss.rst

	if [ "$mod" == 3 ]; then
		echo "Figure ${nFig}: Fragment distribution of $name on TSS." >> ${hist}_fragmentDist_tss.rst
		echo "" >> ${hist}_fragmentDist_tss.rst
	fi

done; done

# Add title lines to rst files
for file in *rst; do
	hist=$(echo $file | awk '{split($0,a,"_");print a[1]}')
	title="Fragment distribution of $hist libraries"
	nChar=${#title}
	line=$(echo $nChar |  awk '{n=$0; l = sprintf("%*s",n,"") ; gsub(/ /,"=",l) ; print l }')
	echo "" | cat - $file > tmp && mv tmp $file
	echo $line | cat - $file > tmp && mv tmp $file
	echo "Fragment distribution of $hist libraries" | cat - $file > tmp && mv tmp $file
done

