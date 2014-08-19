import argparse
import sys
import pysam
import HTSeq
import numpy

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Count reads per nucleosome.')
parser.add_argument('-b',type=str,dest="bFile",help="BAM file. ChIP-seq reads.")
parser.add_argument('-n',type=str,dest="nFile",help="BED file. Nucleosome data, output of Danpos. Only nucleosomes with p-values < 1e-5 are used")
parser.add_argument('-l',type=int,dest="lType",help="INT. Library type: if equal to 0, data is assumed to be paired-end and read lengths are estimated from data itself; if different from 0, data is assumed to be single-end and read lengths use this value extend every read.")
parser.add_argument('-e',type=int,dest="exten",help="INT. Half length of reads extended at their midpoint.")
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file.")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
bamName = args.bFile
nuclFile = args.nFile
libType  = args.lType
extension  = args.exten
oFile = args.oFile

##############################################################
# Functions

def getNucl(nuclFile):
	#Assigns an ID to each nucleosome

	nucl = HTSeq.GenomicArrayOfSets('auto',stranded=False)
	nLine=1
	for line in open(nuclFile,'r'):
		if nLine==1: # skip header line 
			nLine=2
			continue
		line = line.strip('\n')
		fields = line.split('\t')
		# Skip nucleosomes with P-values greater than 5e-5
		if float(fields[5])>0.00005: continue
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])
		iv = HTSeq.GenomicInterval(chrom,start,end,'.')
		nucl[iv] += str(chrom)+"_"+str(start)+"_"+str(end)
	return nucl

def getAvrLength(bamFile):
	# Gets average length of paired-end fragments 

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

def getMaxRegion(nucl,iv):
	# Gets the nucleosome ID with the largest overlap with iv

	# Dictionary, keys: nucleosome iv; values: length of nucleosome iv
	regionLength={}
	for nucl_iv, nucl_ids in nucl[iv].steps():
		length=abs(nucl_iv.end - nucl_iv.start)
		for nucl_id in nucl_ids:
			if nucl_id == None: continue
			if nucl_id in regionLength:
				regionLength[nucl_id] += length
			else:
				regionLength[nucl_id] = length
	# Get largest region of overlapp
	maxLength=0
	maxRegion=None
	for nucl_id in regionLength:
		if regionLength[nucl_id]>maxLength:
			maxLength=regionLength[nucl_id]
			maxRegion=nucl_id
	return maxRegion

def getCounts(avrLength,nucl,bamFile,extension):
	# Counts the number of reads overlapping each nucleosome 

	first,second = None,None
	counts = {}

	for almnt in bamFile:
		# Discard not aligned reads
		if not almnt.aligned: continue
		# Pair end reads
		if almnt.paired_end and almnt.proper_pair and almnt.mate_aligned:
			if first == None:    
				first = almnt
			else: 
				second = almnt
				# Check that first and second reads have same read name
				if first.read.name != second.read.name:
					print "Pair names are different"
					print first.read.name, second.read.name
					first, second = None, None
					continue
				# Define coordenates of new interval
				start = min(first.iv.start,second.iv.start)
				end = max(first.iv.end,second.iv.end)
				midpoint = round((end+start)/2.0)
				start = max(midpoint -  extension,0)
				end = midpoint + extension
				iv = HTSeq.GenomicInterval(first.iv.chrom,start,end,'.')
				# Get ID of the largest region overlapping iv	
				maxRegion=getMaxRegion(nucl,iv)
				# Count the read only in the nucleosome with the largest overlap
				if maxRegion!=None and maxRegion in counts:
					counts[maxRegion] += 1
				elif maxRegion!=None:
					counts[maxRegion] = 1
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
			start = max(midpoint - extension,0)
			end   = midpoint + extension
			iv=HTSeq.GenomicInterval(iv.chrom,start,end,'.')
			# Get the ID of the largest region overlapping iv	
			maxRegion=getMaxRegion(nucl,iv)
			# Count the read only on the nucleosome with the largest overlap
			if maxRegion!=None and maxRegion in counts:
				counts[maxRegion] += 1
			elif maxRegion!=None:
				counts[maxRegion] = 1
	return counts

####################################################
# Excecution

bamFile = HTSeq.BAM_Reader(bamName)

print "Getting nucleosomes"
nucl = getNucl(nuclFile)
print "Getting average read length"
if libType==0: # for paired-end libraries computes fragment length form data itself
	avrLength =  getAvrLength(bamFile)
else:
	avrLength = libType
print avrLength
print "Counting reads per nucleosome"
counts = getCounts(avrLength,nucl,bamFile,extension)
print "Saving results to output file"
out = open(oFile,'w')
nLine=1
for line in open(nuclFile,'r'):
	if nLine==1: #skip header line
		nLine=2
		print >>out, "Avr length = "+str(avrLength)
		continue
	fields=line.split("\t")
	if float(fields[5])>0.00005: continue
	nucID=fields[0]+"_"+str(fields[1])+"_"+str(fields[2])
	if nucID in counts:
		print >>out, nucID,"\t",counts[nucID]
	else:
		print >>out, nucID,"\t0"
out.close()
