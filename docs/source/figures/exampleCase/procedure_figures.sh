##################################################################
# 1. Link to BAM files
##################################################################

# Data summaries files obtained from https://docs.google.com/spreadsheet/ccc?key=0Aueh7dagaPEZdENBUUR1Qk8tS3hhbnZFZ2NyU29CbEE#gid=4

#mkdir bam_rmdup_d0
#mkdir bam_rmdup_d4
#
## Create csv-tables
#awk 'BEGIN{OFS=",";FS="\t"} NR>1{n=split($13,a,"/");link=a[n];$13="`"$2" <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day0/"link">`_";print $1,$13,$2,$6,$9,$7}' dataSummary_day0_tmp.txt
#awk 'BEGIN{OFS=",";FS="\t"} (NR>1 && NR<10){n=split($13,a,"/");link=a[n];$14="`"$3" <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/"link">`_";print $1,$3,$7,$8,$14} NR>=10{n=split($13,a,"/");link=a[n];$14="`"$2" <https://132.239.135.28/public/nucChIP/files/bam_rmdup_day4/"link">`_";print $1,$2,$7,$8,$14}' dataSummary_day4.txt
#
## Create soft links to rmdup BAM files
#cd bam_rmdup_day0
#for file in $(awk 'BEGIN{FS="\t"} NR>1{print $13}' ../dataSummary_day0.txt); do ln -s $file .; done
#ln -s /data2/rivasas2/singleNucleosome/Teif_data/alignments/8_mnase.sort_rmdup.bam .
#cd -
#cd bam_rmdup_day4
#for file in $(awk 'BEGIN{FS="\t"} NR>1{print $14}' ../dataSummary_day4.txt); do ln -s $file .; done
#cd -


###########################################################################
# Shift diagram and dendrograms 


#ln -s /data2/rivasas2/singleNucleosome/tools/thirdparty/docs/source/_static/shift.svg .


#ln -s /data2/rivasas2/singleNucleosome/secondBatch/counts/dendrogram_ratios.svg dendrogram_ratio_counts.svg
#ln -s /data2/rivasas2/singleNucleosome/secondBatch/coverage/dendrograms_ratios.svg dendrogram_ratio_coverage.svg


###########################################################################
# TSS coverage and fragment length 

#for file in /data2/rivasas2/singleNucleosome/secondBatch/nucProfile/*_rmdup_5perc_1500.avrcov.svg; do
#	ln -s $file .
#done

#for file in /data2/rivasas2/singleNucleosome/secondBatch/v-plots/*_1500.vplot.svg; do
#	ln -s $file .
#done
#ln -s /data2/rivasas2/singleNucleosome/Teif_data/v-plots/8_mnase_off_genes.vplot.svg .
#ln -s /data2/rivasas2/singleNucleosome/Teif_data/v-plots/8_mnase_on_genes.3tercile.vplot.svg .

## Combine figures
#declare -a analyses
#analyses["1"]=.nineTiles_rmdup_5perc_1500.avrcov
#analyses["2"]=_on_*vplot
#analyses["3"]=_off_*vplot
#for hm in H3K4me3 H3K27Ac H3K9me3 H3K27me3 8_mnase; do
#	# Horizontal stacking
#	for idx in 1 2 3; do
#		analysis=${analyses[${idx}]}
#		svg_list=""
#		for file in $(ls *${hm}${analysis}.svg);do
#			svg_list=$svg_list" "$file
#		done
#		svg_stack.py --direction=h ${svg_list} > tmp_${idx}.svg
#	done	
#	# Vertical stacking
#	tmp_list=""
#	for file in $(ls tmp*);do
#		tmp_list=${tmp_list}" "${file}
#	done
#	svg_stack.py --direction=v ${tmp_list} > ${hm}_tss.svg
#	rm tmp*
#done

###########################################################################
# lm linear models per nucleosome around TSS
#for file in /data2/rivasas2/singleNucleosome/secondBatch/nucExpProfile/8_mnase_nucLoc/asRatios/*.lm.svg; do
#	ln -s $file .
#done

## Combine figures
#declare -a analyses
#analyses["1"]=.lm
#for hm in H3K4me3 H3K27Ac H3K9me3 H3K27me3; do
#	# Horizontal stacking
#	for idx in 1; do
#		analysis=${analyses[${idx}]}
#		svg_list=""
#		for file in $(ls *${hm}${analysis}.svg);do
#			svg_list=$svg_list" "$file
#		done
#		svg_stack.py --direction=h ${svg_list} > tmp_${idx}.svg
#	done
#	mv tmp_${idx}.svg ${hm}_lm.svg
#done

###########################################################################
# as alternative splicing around included and excluded exons

# Combine figures with MNase normalization #####################################

#mkdir ASexonCovBack # normalized using as background MNase data
#cd ASexonCovBack
#for file in /data2/rivasas2/singleNucleosome/alternativeSplicing/SE/summary_day0/summary/ASexonCovBack/*.svg; do
#	ln -s $file .
#done
#cd -

#folder=ASexonCovBack
#declare -a analyses
#analyses["1"]=_exonStartDB.avrcov
#analyses["2"]=_exonEndDB.avrcov
#for hm in H3K4me3 H3K27Ac H3K9me3 H3K27me3 8_mnase; do
#	# Horizontal stacking
#	for idx in 1 2; do
#		analysis=${analyses[${idx}]}
#		svg_list=""
#		for file in $(ls ${folder}/*${hm}${analysis}.svg);do
#			svg_list=$svg_list" "$file
#		done
#		svg_stack.py --direction=h ${svg_list} > tmp_${idx}.svg
#	done
#	# Vertical stacking
#	tmp_list=""
#	for file in $(ls tmp*);do
#		tmp_list=${tmp_list}" "${file}
#	done
#	svg_stack.py --direction=v ${tmp_list} > ${hm}_as_ratios.svg
#	rm tmp*
#done

# Combine figures without normalization ######################################

#mkdir ASexonCov # without normalization

#cd ASexonCov
#for file in /data2/rivasas2/singleNucleosome/alternativeSplicing/SE/summary_day0/summary/ASexonCov/*.svg; do
#	ln -s $file .
#done
#for file in /data2/rivasas2/singleNucleosome/alternativeSplicing/SE/summary_day0/summary/v-plots/*svg; do
#	ln -s $file .
#done

#cd -

#folder=ASexonCov
#declare -a analyses
#analyses["1"]=_exonStartDB.avrcov
#analyses["2"]=_on_exons_Start.vplot
#analyses["3"]=_off_exons_Start.vplot
#analyses["4"]=_exonEndDB.avrcov
#analyses["5"]=_on_exons_End.vplot
#analyses["6"]=_off_exons_End.vplot
#for hm in H3K4me3 H3K27Ac H3K9me3 H3K27me3 8_mnase; do
#	# Horizontal stacking
#	for idx in 1 2 3 4 5 6; do
#		analysis=${analyses[${idx}]}
#		svg_list=""
#		for file in $(ls ${folder}/*${hm}${analysis}.svg);do
#			svg_list=$svg_list" "$file
#		done
#		svg_stack.py --direction=h ${svg_list} > tmp_${idx}.svg
#	done
#	# Vertical stacking
#	tmp_list=""
#	for file in $(ls tmp*);do
#		tmp_list=${tmp_list}" "${file}
#	done
#	svg_stack.py --direction=v ${tmp_list} > ${hm}_as.svg
#	rm tmp*
#done

#########################################################################3
# Nucleosome figures
#for file in /data2/rivasas2/singleNucleosome/secondBatch/nucLocation/8_mnase/*svg; do
#	ln -s $file .
#done
