import argparse
import sys
import pysam
import HTSeq
import numpy

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Transform BAM reads to resemble the true position of a ChIP seq fragment.')
parser.add_argument('-b',type=str,dest="bFile",help="BAM file. Aligned single or paried-end ChIP-seq reads")
parser.add_argument('-l',type=int,dest="lType",help="INT. Library type: 0 if library is paired-end, fragment length is estimated from data; different from 0 if library is single-end, the value assigned here is used as fragment length.",default=200)
parser.add_argument('-e',type=int,dest="exten",help="INT. Half length of each read around their fragment midpoing.",default=74)
parser.add_argument('-o',type=str,dest="oFile",help="Output file.")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
bamName = args.bFile
libType  = args.lType
extension  = args.exten
oFile = args.oFile

##############################################################
# Functions

def getAvrLength(bamFile):
	# If BAM file is paired-end, outputs the average fragment length
	# This value is then used to extend orphans reads

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

def printBED(avrLength,bamFile,extension,oFile):
	# Take reads in BAM format and print them in BED format: chrom start end 

	first,second = None,None
	out = open(oFile,'w')

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
					print "Pair names are different"
					print first.read.name, second.read.name
					first, second = None, None
					continue
				# Define coordenates of new interval
				start = min(first.iv.start,second.iv.start)
				end = max(first.iv.end,second.iv.end)
				midpoint = round((end+start)/2.0)
				start = midpoint -  extension
				end = midpoint + extension
				if start < 0: start = 0
				iv = HTSeq.GenomicInterval(first.iv.chrom,start,end,'.')
				print >>out, iv.chrom+"\t"+str(iv.start)+"\t"+str(iv.end)+"\t"+almnt.read.name

				# Reset first and second as None
				first = None
				second = None
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
			start = midpoint - extension
			end   = midpoint + extension
			if start < 0: start=0
			iv=HTSeq.GenomicInterval(iv.chrom,start,end,'.')
			print >>out, iv.chrom+"\t"+str(iv.start)+"\t"+str(iv.end)+"\t"+almnt.read.name
	out.close()
	return 1

####################################################
# Excecution

bamFile = HTSeq.BAM_Reader(bamName)

print "Getting average read length"
if libType==0: # if paried-end
	avrLength =  getAvrLength(bamFile)
else:
	avrLength = libType
print avrLength
print "Printing BED"
printBED(avrLength,bamFile,extension,oFile)
