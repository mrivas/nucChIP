import pysam, HTSeq, numpy

#########################################################################
# bam2bed.py
#########################################################################
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

#########################################################################
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
#########################################################################
# getCounts.py
#########################################################################
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

#########################################################################
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

#########################################################################
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

#########################################################################
# getFigures.py 
#########################################################################
def getLibSize(bamName):
	size = pysam.view("-c","-F","4",bamName)
	size = float(size[0])
	return (size)

#########################################################################
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

#########################################################################
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

#########################################################################
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
#########################################################################
def getPlots(out,average,matrix,halfwin):

	fig, ax = pyplot.subplots()
	heatmap = ax.pcolor(matrix,cmap=pyplot.cm.Blues)
	pyplot.savefig(out+"_heatmap.pdf")
	pyplot.close()

	pyplot.plot(numpy.arange(-halfwin,halfwin),profile,label="300 nt")
	pyplot.plot(numpy.arange(-halfwin,halfwin),profile200,label="200 nt")
	pyplot.title(mark)
	pyplot.legend()
	pyplot.savefig(out+"_average.pdf")
	pyplot.close()
