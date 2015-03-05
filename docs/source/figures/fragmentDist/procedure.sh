#rm *svg *rst
#
#for hist in mnase H3K4me3 H3K27Ac H3K9me3 H3K27me3; do
#
#for file in $( ls /data2/rivasas2/singleNucleosome/secondBatch/v-plots/AS_exons/*svg | grep $hist) ; do
#	
#	ln -s $file .
#	file_name=$(echo $file | awk '{n=split($0,a,"/");print a[n]}')
#	name=$(echo $file | awk '{n=split($0,a,"/");split(a[n],b,"_");print b[1]"_"b[2]}')
#	
#	if [ -f ${hist}_fragmentDist.rst ]; then
#		nLine=$(wc -l ${hist}_fragmentDist.rst | awk '{print $1}')
#		mod=$(($nLine%15))
#		nFig=$(( $nLine/15 + 1 )) # every 3 lines (constitutive, spliced_in, spliced_out) add a new figure
#	else
#		mod=0
#		nFig=0
#	fi
#
#	if [ "$mod" == 0 ]; then 
#		echo ".. figure::" >> ${hist}_fragmentDist.rst
#	fi
#
#	echo ".. image:: https://132.239.135.28/public/nucChIP/files/fragmentDist/${file_name}" >> ${hist}_fragmentDist.rst
#	echo "   :width: 45%" >> ${hist}_fragmentDist.rst
#
#	if [ "$mod" == 11 ]; then
#		echo "Figure ${nFig}: Fragment distribution of $name on internal exons." >> ${hist}_fragmentDist.rst
#		echo "" >> ${hist}_fragmentDist.rst
#	fi
#
#done; done
#
## Add title lines to rst files
#for file in *rst; do
#	hist=$(echo $file | awk '{split($0,a,"_");print a[1]}')
#	title="Fragment distribution of $hist libraries"
#	nChar=${#title}
#	line=$(echo $nChar |  awk '{n=$0; l = sprintf("%*s",n,"") ; gsub(/ /,"=",l) ; print l }')
#	echo "" | cat - $file > tmp && mv tmp $file
#	echo $line | cat - $file > tmp && mv tmp $file
#	echo "Fragment distribution of $hist libraries" | cat - $file > tmp && mv tmp $file
#done

##############################################################################
# Add ratios to MNase tab
hist=mnase

echo ".. figure::" >> ${hist}_fragmentDist.rst
for exon_type in all constitutive spliced_in spliced_out; do
for file in $( ls /data2/rivasas2/singleNucleosome/secondBatch/v-plots/AS_exons/ratios/*svg | grep $exon_type) ; do
	
	ln -s $file .
	file_name=$(echo $file | awk '{n=split($0,a,"/");print a[n]}')
	echo ".. image:: https://132.239.135.28/public/nucChIP/files/fragmentDist/${file_name}" >> ${hist}_fragmentDist.rst
	echo "   :width: 45%" >> ${hist}_fragmentDist.rst

done; done

echo "Figure 2: Ratios of fragment distribution of $name on internal exons." >> ${hist}_fragmentDist.rst
echo "" >> ${hist}_fragmentDist.rst

