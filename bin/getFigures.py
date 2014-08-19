import argparse
import sys
import pysam
import HTSeq
import numpy
import cPickle
#import matplotlib.pyplot as pyplot

############################################################
## Get command line arguments
#
#parser = argparse.ArgumentParser(description='Count reads per nucleosome.')
#parser.add_argument('-b',type=str,dest="bFile",help="BAM file.")
#parser.add_argument('-n',type=str,dest="nFile",help="BED file with nucleosome positions.")
#parser.add_argument('-f',type=int,dest="fragL",help="Fragment length if single-end library. If library is paired-end, provide 0 as fragment length.")
#parser.add_argument('-hr',type=int,dest="halfR",help="Half length of reads extended at their midpoint.")
#parser.add_argument('-hw',type=int,dest="halfW",help="Half length of the windows.")
#parser.add_argument('-o',type=str,dest="oFile",help="Output file.")
#
#if len(sys.argv) == 1:
#	print >> sys.stderr,parser.print_help()
#	exit(0)
#
#args = parser.parse_args()
#bamName = args.bFile
#nuclFile = args.nFile
#fragLength  = args.fragL
#halfRead  = args.halfR
#halfwin  = args.halW
#oFile = args.oFile


##############################################################
# Functions

def getLibSize(bamName):
	size = pysam.view("-c","-F","4",bamName)
	size = float(size[0])
	return (size)

def getNuclIV(nuclFile,halfwin):
	nucls = []
	nLine=0
	for line in open(nuclFile,'r'):
		line = line.strip('\n')
		fields = line.split('\t')
		chrom = fields[0]
		start = max(int(fields[1])-halfwin,0)
		end = int(fields[2])+halfwin
		iv = HTSeq.GenomicInterval(chrom,start,end,'.')
		nucls.append(iv)
		nLine += 1
		if nLine==4000: break
	return nucls

def getAvrLength(bamFile):
	
	totalLength,nReads = 0,0
	first,second = None,None

	for almnt in bamFile:
		if not( almnt.aligned and almnt.paired_end and almnt.proper_pair and almnt.mate_aligned):
			continue
		if first == None:
			first = almnt
		else:
			second = almnt
			# Check that first and second have same read name
			if first.read.name != second.read.name:
#				print "Pair names are different"
#				print first.read.name, second.read.name
				continue

			start = min(first.iv.start,second.iv.start)
			end = max(first.iv.end,second.iv.end)
			length = end-start
			# Check length size
			if length<=0:
				print "negative length"
				print first.read.name
				print second.rad.name
				first = None
				second = None
				continue

			totalLength += length
			nReads +=1

			first = None
			second = None
	
	return round(float(totalLength)/float(nReads))

def getFullCoverage(avrLength,bamFile,halfRead):

	coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )

	first,second = None,None
	for almnt in bamFile:
		# Discard not aligned reads
		if not almnt.aligned: continue
		# Pair end reads
		if almnt.paired_end and almnt.proper_pair and almnt.mate_aligned:
			if first == None:    
				first = almnt
			else: 
				second = almnt
				# Check that first and second have same read name
				if first.read.name != second.read.name:
					#print "Pair names are different"
					#print first.read.name, second.read.name
					first, second = None, None
					continue
				# Define coordenates of new interval
				start = min(first.iv.start,second.iv.start)
				end = max(first.iv.end,second.iv.end)
				midpoint = round((end+start)/2.0)
				start = midpoint -  halfRead
				end = midpoint + halfRead
				if start < 0: start = 0
				iv = HTSeq.GenomicInterval(first.iv.chrom,start,end,'.')
				# Add to coverage
				coverage[iv] += 1
				# Reset first and second as None
				first, second = None, None
		# Orphan reads (only one read in a pair is aligned), or single-end reads
		elif not almnt.mate_aligned:
			# Define coordenates of new interval
			iv = almnt.iv
			if iv.strand == "+":
				start = iv.start
				end   = iv.start + avrLength
			elif iv.strand == "-":
				end = iv.end
				start = iv.end - avrLength
			midpoint=round((end+start)/2.0)
			start = midpoint - halfRead
			end   = midpoint + halfRead
			if start < 0: start=0
			iv=HTSeq.GenomicInterval(iv.chrom,start,end,'.')
			# Add to coverage
			coverage[iv] += 1
			# Reset first and second as None
			first, second = None, None

	return coverage

def getProfiles(coverage,halfwin,nucl_ivs,libSize):
	# Computes average and matrix profiles of coverage
	average=numpy.zeros(2*halfwin, dtype='d')
	matrix=None
	for nucl_iv in nucl_ivs:
		wincvg = numpy.fromiter(coverage[nucl_iv],dtype="i",count=2*halfwin)*float(1e6)/libSize	
		average += wincvg
		if matrix==None:
			matrix=numpy.matrix( wincvg )
		else:
			matrix=numpy.vstack([matrix,wincvg] )
	
	average= average/len(nucl_ivs)

	return [average, matrix]

		
#def getPlots(out,average,matrix,halfwin):
#
#	fig, ax = pyplot.subplots()
#	heatmap = ax.pcolor(matrix,cmap=pyplot.cm.Blues)
#	pyplot.savefig(out+"_heatmap.pdf")
#	pyplot.close()
#
#	pyplot.plot(numpy.arange(-halfwin,halfwin),profile,label="300 nt")
#	pyplot.plot(numpy.arange(-halfwin,halfwin),profile200,label="200 nt")
#	pyplot.title(mark)
#	pyplot.legend()
#	pyplot.savefig(out+"_average.pdf")
#	pyplot.close()
	
####################################################
# Excecution


bamName = "/data2/rivasas2/singleNucleosome/enriched_histones/secondBatch/pooled/all_H3K27Ac_d0.sortedByName.bam"
nuclFile = "/data2/rivasas2/singleNucleosome/enriched_histones/secondBatch/pooled/H3K27Ac_enriched.bed"
fragLength  = 150
halfRead  = 74
halfwin  = 3000
out = "test"


bamFile = HTSeq.BAM_Reader(bamName)

print "Getting library size"
#libSize=getLibSize(bamName)
libSize=20000000
print "Getting nucleosomes"
nucl_ivs = getNuclIV(nuclFile,halfwin)
print "Getting average read length"
if fragLength!=0: # if single-end
	avrLength = fragLength
else:
	avrLength =  getAvrLength(bamFile)
print avrLength
print "Producing coverage"
coverage = getFullCoverage(avrLength,bamFile,halfRead)
print "Producing average and matrix profiles"
average,matrix = getProfiles(coverage,halfwin,nucl_ivs,libSize)
print "Plotting results"
#getPlots(out,average,matrix,halfwin)
print "Saving results to output file"
outFile=open(out+".cpickle",'wb')
cPickle.dump(obj,outFile)
outFile.close()
##To load the objects:
#outFile=file(out+".cpickle","rb")
#loaded_obj=[]
#for i in range(1):
#	loaded_objects.append(cPickle,load(f) )
#f.close()
