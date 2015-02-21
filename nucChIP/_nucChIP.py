#!/usr/bin/env python
import pysam, HTSeq, numpy, pickle, itertools, operator, pybedtools, wLib, pandas
from numpy import inf
from numpy import nan
from scipy import signal
numpy.random.seed(seed=10)
#########################################################################
# bam2bed.py
#########################################################################
# Deprecated. Schedule for DELETION
##def getAvrLength(bamFile):
##	# If BAM file is paired-end, outputs the average fragment length
##	# This value is then used to extend orphans reads
##
##	
##	totalLength,nReads = 0,0
##	first,second = None,None
##
##	for almnt in bamFile:
##		if not( almnt.aligned and almnt.paired_end and almnt.proper_pair and almnt.mate_aligned):
##			continue
##		if first == None:
##			first = almnt
##		else:
##			second = almnt
##			# Check that first and second have same read name
##			if first.read.name != second.read.name:
###				print "Pair names are different"
###				print first.read.name, second.read.name
##				continue
##
##			start = min(first.iv.start,second.iv.start)
##			end = max(first.iv.end,second.iv.end)
##			length = end-start
##			# Check length size
##			if length<=0:
##				print "negative length"
##				print first.read.name
##				print second.rad.name
##				first = None
##				second = None
##				continue
##
##			totalLength += length
##			nReads +=1
##
##			first = None
##			second = None
##	
##	return round(float(totalLength)/float(nReads))

#########################################################################
def printBED(avrLength,bamName,extension,oFile,lower,upper):
	# Take reads in BAM format and print them in BED format: chrom start end 

	bamFile = pysam.Samfile(bamName,'rb')
	out = open(oFile,'w')
	for almnt in bamFile:
		# Red filter
		if avrLength==0: # Paired-end
			if almnt.is_unmapped or not(almnt.is_read1) or not(almnt.is_proper_pair) or abs(almnt.isize)>upper or abs(almnt.isize)<lower or almnt.mapq<20: continue
		else: # single-end
			if almnt.is_unmapped or almnt.mapq<20: continue

		chrom=bamFile.getrname(almnt.rname)
		if not(almnt.is_reverse):
			strand="+"
		else:
			strand="-"
		
		if avrLength==0: # Paired-end
			fragLength = abs(almnt.isize)
		else:           # if Single-end
			fragLength = avrLength
		if not(almnt.is_reverse):
			almntStart = almnt.pos
			almntEnd   = almnt.pos + fragLength
		else:
			almntStart = almnt.aend - fragLength
			almntEnd   = almnt.aend
		
		midpoint=numpy.mean([almntStart,almntEnd])
		start = int( max(midpoint - extension,0) )
		end   = int( midpoint  + extension       )
		
		output=[ chrom, start,end, almnt.qname, almnt.mapq, strand ]
		print >>out, "\t".join(map(str,output))

	out.close()
	return 1
#########################################################################
def printBAM(avrLength,bamName,extension,oFile):
	# Take reads in BAM format and print them in BED format: chrom start end 

	bamFileRef = pysam.Samfile(bamName,'rb')
	out = pysam.Samfile(oFile,'wbu',header=bamFileRef.header)
	
	bamFile = HTSeq.BAM_Reader(bamName)
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
					print "Pair names are different"
					print first.read.name, second.read.name
					first, second = None, None
					continue
				# Define coordenates of new interval
				start = min(first.iv.start,second.iv.start)
				end = max(first.iv.end,second.iv.end)
				midpoint = round((end+start)/2.0)
				start = max(0,midpoint -  extension)
				a = pysam.AlignedRead()
				a.qname = almnt.read.name
				a.flag  = 0
				a.rname = bamFileRef.gettid(iv.chrom)
				a.pos   = start
				a.mapq  = 254
				a.cigar = [(0,2*extension)]
				a.mrnm  = -1
				a.mpos  = -1
#				a.tlen = avrLength
#				a.seq   = "*"
#				a.qual  = "*"
				a.setTag('AL', avrLength, value_type='i', replace=True)
				out.write(a)

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
			start = max(0,midpoint - extension)
			a = pysam.AlignedRead()
			a.qname = almnt.read.name
			a.flag  = 0
			a.rname = bamFileRef.gettid(iv.chrom)
			a.pos   = start
			a.mapq  = 254
			a.cigar = [(0,2*extension)]
			a.mrnm  = -1
			a.mpos  = -1
#			a.tlen = avrLength
#			a.seq   = "*"
#			a.qual  = "*"
			a.setTag('AL', avrLength, value_type='i', replace=True)
			out.write(a)
	bamFileRef.close()
	out.close()
	return 1
#########################################################################
# getCounts.py
#########################################################################
def getNucl(nuclFile,pValue):
	#Assigns an ID to each nucleosome

	nucl = HTSeq.GenomicArrayOfSets('auto',stranded=False)
	nLine=1
	for line in open(nuclFile,'r'):
		fields = line.strip().split('\t')
		if len(fields)==1 or not( fields[1].isdigit()): continue # skip header lines
#		# Skip nucleosomes with P-values greater than 5e-5
		if float(fields[8])<pValue: continue
		chrom = fields[0]
		start_original = int(fields[1])
		end_original   = int(fields[2])
		midpoint = numpy.mean([start_original,end_original])
		start = int(max(midpoint -  75,0))
		end   = int(midpoint + 75)
		iv = HTSeq.GenomicInterval(chrom,start,end,'.')
		nucl[iv] += str(chrom)+"_"+str(start_original)+"_"+str(end_original)
	return nucl

#########################################################################
def getMaxRegion(nucl,iv):
	# Gets the nucleosome ID with the largest iv overlap

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
def getCounts(avrLength,nucl,bamName,extension,lower,upper):
	# Counts the number of reads overlapping each nucleosome 
	bamFile = pysam.Samfile( bamName, 'rb')
	counts = {}

	for almnt in bamFile:
		if avrLength==0: # library is paired-end
			if almnt.is_unmapped or almnt.mapq<20 or not(almnt.is_read1) or not(almnt.is_proper_pair) or abs(almnt.isize)>upper or abs(almnt.isize)<lower: continue
			fragLength=abs(almnt.isize)
		else: # library is single-end
			if almnt.is_unmapped or almnt.mapq<20: continue
			fragLength = avrLength

		if not almnt.is_reverse: # If read is on '+' strand
			almntStart = almnt.pos 
			almntEnd  = almnt.pos + fragLength
		else:                    # If read is on '-' strand
			almntStart = almnt.aend - fragLength 
			almntEnd  = almnt.aend

		# Define coordenates of new interval
		chrom=bamFile.getrname(almnt.rname)
		midpoint = numpy.mean([almntEnd,almntStart])
		start = max(midpoint -  extension,0)
		end = midpoint + extension
		iv = HTSeq.GenomicInterval(chrom,start,end,'.')
		# Get ID of the largest region overlapping iv	
		maxRegion=getMaxRegion(nucl,iv)
		# Count the read only in the nucleosome with the largest overlap
		if maxRegion!=None and maxRegion in counts:
			counts[maxRegion] += 1
		elif maxRegion!=None:
			counts[maxRegion] = 1
	return counts

#########################################################################
def lineToIv(line):
	fields = line.split('\t')
	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	midpoint = numpy.mean([start,end])
	start = max(midpoint -  75,0)
	end = midpoint + 75
	iv = HTSeq.GenomicInterval(chrom,start,end,'.')
	return iv

#########################################################################
def maxOverlap(almnt_iv,regions_dic):
    
	regionsSet = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	for key in regions_dic:
		chrom=regions_dic[key].chrom
		start=regions_dic[key].start
		end  =regions_dic[key].end
		region_iv=HTSeq.GenomicInterval( chrom, start, end)
		regionsSet[ region_iv ] += key

	keyLength = {}
	for iv, key_set in regionsSet[ almnt_iv ].steps():
		if len(key_set)==0: continue
		for key in key_set:
			if key in keyLength:
				keyLength[key] += iv.end-iv.start
			else:
				keyLength[key] = iv.end-iv.start
	
	if len(keyLength)==0:
		maxKey="none"
	else:
		maxKey=max(keyLength.iteritems(), key=operator.itemgetter(1))[0]				
	if maxKey=='curr':
		return True
	else:
		return False				

#########################################################################
def getCounts2(avrLength,nucName,bamName,extension,lower,upper,oFile,headerLine):
	# Counts the number of reads overlapping each nucleosome 
	
	# Determine number of lines on nucFile
	nucFile=open(nucName,'r')
	nLines=0
	for line in nucFile: nLines += 1
	nucFile.seek(0)

	bamFile = pysam.Samfile( bamName, 'rb')
	out=open(oFile,'w')
	# Count reads per nucleosome
	for idx in range(nLines):
		# Print header line
		if (idx+1)<headerLine:
			line=nucFile.readline().strip()
			continue
		elif (idx+1)==headerLine:
			line=nucFile.readline().strip()
			fields = line.strip().split('\t')
			output = fields + [ "count"  ]
			print >>out, "\t".join(map(str,output))
			continue		
		# Print counts
		regions={}
		if idx==(headerLine):
			lineCurr = nucFile.readline().strip()
			lineDown = nucFile.readline().strip()
			regions['curr']	= lineToIv( lineCurr )
			regions['down'] = lineToIv( lineDown )
		elif idx==(nLines-1):
			lineUpst = lineCurr
			lineCurr = lineDown
			regions['upst'] = lineToIv( lineUpst )
			regions['curr'] = lineToIv( lineCurr )
		else:
			lineUpst = lineCurr
			lineCurr = lineDown
			lineDown = nucFile.readline().strip()
			regions['upst'] = lineToIv( lineUpst )
			regions['curr']	= lineToIv( lineCurr )
			regions['down'] = lineToIv( lineDown )

		counts = 0
		for almnt in bamFile.fetch( regions['curr'].chrom,max(regions['curr'].start-200,0),regions['curr'].end+200 ):

			if avrLength==0: # library is paired-end
				if almnt.is_unmapped or almnt.mapq<20 or not(almnt.is_read1) or not(almnt.is_proper_pair) or abs(almnt.isize)>upper or abs(almnt.isize)<lower: continue
				fragLength=abs(almnt.isize)
			else: # library is single-end
				if almnt.is_unmapped or almnt.mapq<20: continue
				fragLength = avrLength

			if not almnt.is_reverse: # If read is on '+' strand
				almntStart = almnt.pos 
				almntEnd  = almnt.pos + fragLength
			else:                    # If read is on '-' strand
				almntStart = almnt.aend - fragLength 
				almntEnd  = almnt.aend

			# Define coordenates of new interval
			chrom=bamFile.getrname(almnt.rname)
			midpoint = numpy.mean([almntEnd,almntStart])
			start = max(midpoint -  extension,0)
			end = midpoint + extension
			almnt_iv = HTSeq.GenomicInterval(chrom,start,end,'.')
			# Get ID of the largest region overlapping iv
			if not( maxOverlap(almnt_iv,regions) ): continue
			# Count the read only in the nucleosome with the largest overlap
			counts += 1
		output = lineCurr.split('\t') + [counts]
		print >>out, "\t".join(map(str,output))
	out.close()
	nucFile.close()

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

##############################################################
# getProfile
##############################################################
def getRegions(geneList,database,halfwinwidth): 
	tss = pickle.load(open(database,'rb'))
	regions = []
	noPresent = 0
	for gene_id in open(geneList,'r'):
		gene_id = gene_id.strip()
		if not gene_id in tss: 
			noPresent += 1 
			continue
		chrom = tss[gene_id][0]
		start = tss[gene_id][1]
		end   = tss[gene_id][2]
		strand = tss[gene_id][5]
		midpoint = round((start+end)/2.0)
		start = int(max( midpoint-halfwinwidth, 0))
		end   = int(midpoint+halfwinwidth)
		window = HTSeq.GenomicInterval(chrom,start,end,strand)
		regions.append( window )
	print str(noPresent)+" genes not present on database of TSS"
	return regions
################################################################
def getProfile(halfwinwidth,regions,bamName,fragLength,lower,upper,pcount):
	bamFile = pysam.Samfile( bamName, 'rb')
	libSize = 0
	for almnt in bamFile: 
		if fragLength==0: # if lib is paired-end
			if not( almnt.is_unmapped) and almnt.is_read1 and almnt.is_proper_pair and abs(almnt.isize)<=upper and abs(almnt.isize)>=lower and almnt.mapq>=20:
				libSize += 1
		else: # if lib is single-end
			if not( almnt.is_unmapped) and almnt.mapq>=20: 
				libSize += 1

	if type(regions)==dict:
		nucArrayCoverage = {}
		for gene_id in regions:
			regionsList = regions[gene_id]
			nRegions = len( regionsList )
			profileMatrix = numpy.zeros( ( nRegions, 2*halfwinwidth ), dtype="d" )
			for row,region in enumerate(regionsList):
				if bamFile.gettid(region.chrom)==-1: continue
				for almnt in bamFile.fetch( region.chrom,max(region.start-200,0),region.end+200 ):
#					if almnt.is_unmapped or almnt.pos < 0: continue
					if fragLength==0: # If lib is paried-end
						if almnt.is_unmapped or almnt.pos < 0 or not(almnt.is_read1) or not(almnt.is_proper_pair) or abs(almnt.isize)>upper or abs(almnt.isize)<lower or bamFile.gettid(region.chrom)==-1 or almnt.mapq<20: continue
						if not almnt.is_reverse: # If read is on '+' strand
							almntStart = almnt.pos 
							almntEnd   = almnt.pos + abs(almnt.isize)
						else:                    # If read is on '-' strand
							almntStart = almnt.aend - abs(almnt.isize) 
							almntEnd   = almnt.aend
					else: # If lib is single-end
						if almnt.is_unmapped or bamFile.gettid(region.chrom)==-1 or almnt.mapq<20: continue
						if not almnt.is_reverse: # If read is on '+' strand
							almntStart = almnt.pos 
							almntEnd   = almnt.pos + fragLength
						else:                    # If read is on '-' strand
							almntStart = almnt.aend - fragLength 
							almntEnd   = almnt.aend
						
					if almntStart > region.end or almntEnd < region.start: continue
					if region.strand == "+":
						start_in_window = almntStart - region.start
						end_in_window   = almntEnd   - region.start
					else:
						start_in_window = -almntEnd   + region.end
						end_in_window   = -almntStart + region.end 
					start_in_window = max( start_in_window, 0 )
					end_in_window = min( end_in_window, 2*halfwinwidth )
					profileMatrix[ row, start_in_window : end_in_window ] += 1
		
			# Normalized by library size
			profileMatrix = (profileMatrix + pcount ) / libSize * 1e6 
			meanCoverage  = numpy.mean( profileMatrix, axis=1 )
			nucArrayCoverage[ gene_id ] = meanCoverage
		bamFile.close()
		return nucArrayCoverage

	elif type(regions)==list:
		nRegions = len( regions )
		profileMatrix = numpy.zeros( ( nRegions, 2*halfwinwidth ), dtype="d" )
		for row,region in enumerate(regions):
			if bamFile.gettid(region.chrom)==-1: continue
			for almnt in bamFile.fetch( region.chrom,max(region.start-200,0),region.end+200 ):
#				if almnt.is_unmapped or almnt.pos < 0: continue
				# If library is PAIRED-END
				if fragLength==0:
					if almnt.is_unmapped or almnt.pos < 0  or not(almnt.is_read1) or not(almnt.is_proper_pair) or abs(almnt.isize)>upper or abs(almnt.isize)<lower or bamFile.gettid(region.chrom)==-1 or almnt.mapq<20: continue
					if not almnt.is_reverse: # If read is on '+' strand
						almntStart = almnt.pos 
						almntEnd  = almnt.pos + abs(almnt.isize)
					else:                    # If read is on '-' strand
						almntStart = almnt.aend - abs(almnt.isize) 
						almntEnd  = almnt.aend
				# If library is SINGLE-END
				else:
					if almnt.is_unmapped or almnt.pos < 0 or  bamFile.gettid(region.chrom)==-1 or almnt.mapq<20: continue
					if not almnt.is_reverse: # If read is on '+' strand
						almntStart = almnt.pos 
						almntEnd  = almnt.pos + fragLength
					else:                    # If read is on '-' strand
						almntStart = almnt.aend - fragLength 
						almntEnd  = almnt.aend
					
				if almntStart > region.end or almntEnd < region.start: continue
				if region.strand == "+":
					start_in_window = almntStart - region.start
					end_in_window   = almntEnd   - region.start
				else:
					start_in_window = -almntEnd   + region.end
					end_in_window   = -almntStart + region.end 
				start_in_window = max( start_in_window, 0 )
				end_in_window = min( end_in_window, 2*halfwinwidth )
				profileMatrix[ row, start_in_window : end_in_window ] += 1
		bamFile.close()
		# Normalized by library size
		profileMatrix = ( profileMatrix + pcount ) / libSize * 1e6 
		return profileMatrix
#######################################################
# exprHistCorr
#######################################################
def getGeneList( geneListData):
	geneList={}
	nLine = 0
	for line in open(geneListData,'r'): 
		nLine += 1
		gene_id=line.strip()
		geneList[ gene_id ] = 1
	return geneList


def getExpression( expressionData ):
	expression={}
	nLine=0
	for line in open(expressionData,'r'): 
		nLine += 1
		if nLine==1: continue # skips header line
		line=line.strip().split('\t')
		gene_id = line[0]
		fpkm	= line[6]
		expression[ gene_id ] = fpkm
	return expression

# Compute nucleosomes IVs and counts per gene
def getNucsCounts( countsData ):
	nucleosomes = HTSeq.GenomicArrayOfSets('auto', stranded=False)
	counts = {}
	nLine=0
	for line in open(countsData,'r'): 
		nLine += 1
		if nLine==1: continue # skips header line
		line=line.strip().split('\t')
		chrom, start, end = line[0].split('_')
		region = HTSeq.GenomicInterval( chrom, int(start), int(end), '.' ) 
		nucleosomes[ region ] += region
		counts[ region ] = int( line[1] )
	return [ nucleosomes, counts ]

# Determine nucleosome positions around each gene TSS
def getTssNucs( geneList, tssData, nucleosomes):
	# Determine nucleosomes around TSS
	tss=pickle.load( open(tssData, 'r') )
	
	tssNucs = {}
	for gene_id in geneList:
		if not gene_id in tss: continue
		chrom  = tss[gene_id][0]
		pos	= int( tss[gene_id][2] )
		strand = tss[gene_id][5]
		if strand == "+":
			region = HTSeq.GenomicInterval( chrom, pos-1000, pos+2000, '.' )
		else:
			region = HTSeq.GenomicInterval( chrom, pos-2000, pos+1000, '.' )
		# Save all nucleosomes into a set to avoid redundancies
		nucleosomeSet = set()
		for iv, nucleosome in nucleosomes[ region ].steps():
			nucleosomeSet |= nucleosome
		# Save nucleosomes as arrays
		minusNucDist, minusNucRegion = [],[]
		plusNucDist,  plusNucRegion  = [],[]
		for nucleosome in nucleosomeSet:	
			distance = numpy.mean([nucleosome.start,nucleosome.end]) - pos
			if strand == "-": distance = -distance
			if distance < 0:
				minusNucDist.append( distance )
				minusNucRegion.append( nucleosome )
			elif distance >= 0:
				plusNucDist.append( distance )
				plusNucRegion.append( nucleosome )
		# Determine relative position of nucleosomes around TSS	
		minusNucIdx, plusNucIdx = [], []
		if len(minusNucDist) >= 1:
			minusNucIdx = numpy.argsort( minusNucDist )[-1::]
		if len(plusNucDist)>= 1 and len(plusNucDist)<=4:
			plusNucIdx  = numpy.argsort( plusNucDist )
		elif len(plusNucDist)>=1 and len(plusNucDist)>4:
			plusNucIdx = numpy.argsort(plusNucDist)[4::]
		# Extract -1, +1, +2, +3, +4 nucleosomes
		nucRegions = {}
		if len(minusNucIdx)  > 0: nucRegions[ 'neg1' ] = minusNucRegion[ minusNucIdx ]
		if len(plusNucIdx)  >= 1:
			labels = [ 'pos1', 'pos2', 'pos3', 'pos4' ]
			for Nidx,idx in enumerate( plusNucIdx ):
				if Nidx > 3: break
				label = labels[Nidx]
				nucRegions[ label ] = plusNucRegion[ idx ]

		tssNucs[ gene_id ] = nucRegions
	
	return tssNucs	   

##################################################################################	
# getPromoterCounts
##################################################################################	
def getNuclPerPromoter(tssData, nucRelPosFile ):
	
	names, nucRelPos = [],[]
	for line in open( nucRelPosFile, 'r'):
		line=line.strip().split("\t")
		ID = line[0]
		pos = round(float( line[1] ))
		names.append( ID )
		nucRelPos.append( pos )
	
	tss=pickle.load( open(tssData, 'r') )
	nuclArray = {}
	for gene_id in tss:
		chrom  = tss[gene_id][0]
		tssPos	= int( tss[gene_id][2] )
		strand = tss[gene_id][5]
		
		nuclArray[ gene_id ] = []
		for nucPos in nucRelPos:
			if strand == "+":
				newPos = tssPos + nucPos
			else:
				newPos = tssPos - nucPos
			start = max( 0, newPos-74 ) 
			end   = newPos + 74
			nuclArray[gene_id].append( HTSeq.GenomicInterval( chrom, start, end, '+' ) )

	return [names,nuclArray]
####################################################################################3
def getNuclPerPromoter2(tssData, nucPosFileName ):
	
	#nucPosFile=pybedtools.BedTool(nucPosFileName)
	nucPosFile=wLib.BigWigFile(nucPosFileName)
	tss=pickle.load( open(tssData, 'r') )
	nuclArray, nuclArrayDists = {},{}
	for gene_id in tss:
		chrom  = tss[gene_id][0]
		tssPos	= int( tss[gene_id][2] )
		strand = tss[gene_id][5]
		
		halfwin=2000
		nuclArray[ gene_id ] = numpy.array([])
		nuclArrayDists[ gene_id ] = numpy.array([])
		while len(nuclArray[gene_id])<5 and halfwin<=5000:
			#iv = pybedtools.Interval( chrom, max(tssPos-halfwin,0), tssPos+halfwin )
			#for nucRegion in nucPosFile.tabix_intervals( iv ):
			for nucRegion in nucPosFile.fetch( chrom, max(tssPos-halfwin,0), tssPos+halfwin ):
				#nucPos = round(float(numpy.mean([nucregion.start,nucregion.end]))) 
				nucPos = round(nucRegion.score) 
				start = nucPos - 74
				end   = nucPos + 74
				nuclArray[gene_id]=numpy.append( nuclArray[gene_id], HTSeq.GenomicInterval( chrom, start, end, '+' ) )
				nuclArrayDists[gene_id]=numpy.append( nuclArrayDists[gene_id], nucPos-tssPos )
			if strand == "-":
				nuclArray[gene_id] = nuclArray[gene_id][::-1]
				nuclArrayDists[gene_id]   = -nuclArrayDists[gene_id][::-1]
			if sum((nuclArrayDists[gene_id]>=0))<3 or sum((nuclArrayDists[gene_id]<0))<2:
				nuclArray[ gene_id ] = numpy.array([])
				nuclArrayDists[ gene_id ] = numpy.array([])
				halfwin += 1000
				continue
			posNuc = nuclArray[gene_id][ (nuclArrayDists[gene_id]>=0) ][:3]
			negNuc = nuclArray[gene_id][ (nuclArrayDists[gene_id]< 0) ][-2::]
			posDists = nuclArrayDists[gene_id][ (nuclArrayDists[gene_id]>=0) ][:3]
			negDists = nuclArrayDists[gene_id][ (nuclArrayDists[gene_id]< 0) ][-2::]
			nuclArray[ gene_id ] = list(numpy.append(negNuc,posNuc)) 
			nuclArrayDists[ gene_id ] = list(numpy.append(negDists, posDists))
		if halfwin>5000: # discard genes with non-well stablished promoter nucleosomes
			nuclArray.pop(gene_id,None)
			nuclArrayDists.pop(gene_id,None)
			
	return [ nuclArray, nuclArrayDists ]
####################################################################################	
# nucLocPrediction
####################################################################################	
def getSumLines(File,chrom,end,pseudocount):
	nLines, sumTotal = 0.0, 0.0
	for line in File.tabix_intervals(pybedtools.Interval(chrom, 3000000,end)):
		nLines = nLines + 1
		sumTotal += float(line[3])+ pseudocount
	return [sumTotal,nLines]

def getSumLog(File,chrom,end,mu,pseudocount):
	sumLogTotal = 0.0
	for line in File.tabix_intervals(pybedtools.Interval(chrom, 3000000,end)):
		sumLogTotal += numpy.log( (float(line[3])+pseudocount)/mu )
	return sumLogTotal
####################################################################################	
# simulateRegularChIP
####################################################################################	
def printNewAlmnt(bamInName,bamOutName,mu,sd):

	bamIn  = pysam.Samfile(bamInName,'rb')
	bamOut = pysam.Samfile(bamOutName,'wbu',header=bamIn.header)
	
	for almnt in bamIn:
		# Filter out seconds, unmapped, and low quality reads
		if not(almnt.is_read1) or almnt.is_unmapped or almnt.mapq<20 or not(almnt.is_proper_pair): continue

		# New fragment length
		fragLength = int(numpy.random.normal(mu,sd))
		while fragLength < 100: fragLength = numpy.random.normal(mu,sd) # makes sure new fragment lenght >= 100
		# Redefine reads
		diff=int(round((fragLength - abs(almnt.isize))/2))
		if not(almnt.is_reverse): # + strand
			pos = almnt.pos - diff
		else:
			pos = (almnt.pos + almnt.cigar[0][1] - abs(almnt.isize)) - diff
		a = almnt
		a.pos = max(pos,0)
		a.cigar = [(0,fragLength)]
		bamOut.write(a)

	bamIn.close()
	bamOut.close()
########################################################################################
# fragProfile
########################################################################################

def getFragProfile(halfwinwidth,regions,bamName):
	bamFile = pysam.Samfile( bamName, 'rb')

	nRegions = len( regions )
	profileMatrix = numpy.zeros( ( nRegions, 2*halfwinwidth ), dtype="d" )
	for row,region in enumerate(regions):
		profileCounts = numpy.zeros( 2*halfwinwidth , dtype="d" )
		if bamFile.gettid(region.chrom)==-1: continue
		for almnt in bamFile.fetch( region.chrom,region.start,region.end ):
			if almnt.is_unmapped or almnt.pos < 0 or not( almnt.is_read1) or not(almnt.is_proper_pair): continue
			if not almnt.is_reverse:
				almntStart = almnt.pos 
				almntEnd  = almnt.pos + abs(almnt.isize)
			else:
				almntStart = almnt.aend - abs(almnt.isize) 
				almntEnd  = almnt.aend
			if region.strand == "+":
				start_in_window = almntStart - region.start
				end_in_window   = almntEnd   - region.start
			else:
				start_in_window = -almntEnd   + region.end
				end_in_window   = -almntStart + region.end 
			start_in_window = max( start_in_window, 0 )
			end_in_window = min( end_in_window, 2*halfwinwidth )
			profileMatrix[ row, start_in_window : end_in_window ] += abs(almnt.isize)
			profileCounts[      start_in_window : end_in_window ] += 1
		profileMatrix[row,:] = numpy.divide(profileMatrix[row,:],profileCounts)
		profileMatrix[row,:][ profileMatrix[row,:]== inf ] = nan
		profileMatrix[row,:][ profileMatrix[row,:]==-inf ] = nan
	bamFile.close()
	
	return profileMatrix
########################################################################################
# fragDistribution
########################################################################################

def getFragDistribution(bamName):
	bamFile = pysam.Samfile( bamName, 'rb')

	lengths = []
	for almnt in bamFile:
		if almnt.is_unmapped or almnt.pos < 0 or not( almnt.is_read1) or not(almnt.is_proper_pair): continue
		lengths.append(	abs(almnt.isize) )
	bamFile.close()
	
	return lengths
########################################################################################
# vPlot
########################################################################################

def getSizeProfile(halfwinwidth,regions,bamName,lower,upper):
	bamFile = pysam.Samfile( bamName, 'rb')
	
	nRows=upper-lower+1
	profileMatrix = numpy.zeros( ( nRows, 2*halfwinwidth ), dtype="d" )
	# Iterate over regions
	for region in regions:
		# Check if chromosome is on the BAM file
		if bamFile.gettid(region.chrom)==-1: continue
		for almnt in bamFile.fetch( region.chrom,region.start-200,region.end+200 ):
			if almnt.is_unmapped or not( almnt.is_read1) or almnt.mapq<20: continue
			# Discard fragments larger than 200 nt
			if abs(almnt.isize)>200 or abs(almnt.isize)<lower or abs(almnt.isize)>upper: continue
			if not almnt.is_reverse: # If read is on '+' strand
				almntStart = almnt.pos 
				almntEnd  = almnt.pos + abs(almnt.isize)
			else:                    # If read is on '-' strand
				almntStart = almnt.aend - abs(almnt.isize) 
				almntEnd  = almnt.aend

			if almntStart>region.end or almntEnd<region.start: continue 
			if region.strand == "+":
				start_in_window = almntStart - region.start
				end_in_window   = almntEnd   - region.start
			else:
				start_in_window = -almntEnd   + region.end
				end_in_window   = -almntStart + region.end 
			start_in_window = max( start_in_window, 0 )
			end_in_window = min( end_in_window, 2*halfwinwidth )
			row = upper-int(abs(almnt.isize))
			profileMatrix[ row, start_in_window : end_in_window ] += 1 
	
	bamFile.close()
	
	return profileMatrix

#######################################################################################	
# exonProfile
#######################################################################################	
def getExons( exonList, misoSummary, halfwinwidth ):
	miso = pickle.load(open(misoSummary,'rb') )
	nRows=0
	for line in open(exonList,'r'): nRows+=1
	exonMatrix = numpy.zeros( (nRows,6), dtype='int' )
	regions = []
	for row,line in enumerate( open(exonList,'r') ):
		exon_id = line.strip()
		if not(exon_id in miso): continue
		exonMatrix[row,:] = miso[ exon_id ][0:6]
		start =  int(miso[exon_id][0] - halfwinwidth)
		end   =  int(miso[exon_id][5] + halfwinwidth)
		chrom  = miso[ exon_id ][6]
		strand = miso[ exon_id ][7]
		iv = HTSeq.GenomicInterval( chrom, start, end, strand )
		regions.append( iv )
	return regions, exonMatrix	
			

def getSamplePoints( exons, halfwinwidth, numSamples ):
    start, end = 0,halfwinwidth
    samplePoints = numpy.linspace(start,end, numSamples,endpoint=False).astype('int')
    n = len( exons )
    for i in range( n-1 ):
        start = exons[i]   - exons[0] + halfwinwidth
        end   = exons[i+1] - exons[0] + halfwinwidth
        points = numpy.linspace( start, end, numSamples,endpoint=False).astype('int')
        samplePoints = numpy.append(samplePoints, points)
    start = exons[n-1] -                 exons[0] + halfwinwidth
    end   = exons[n-1] + halfwinwidth  - exons[0] + halfwinwidth
    points = numpy.linspace( start, end, numSamples,endpoint=False ).astype('int')
    samplePoints = numpy.append(samplePoints, points )

    return samplePoints

def scaleProfile( profile, exons, halfwinwidth, numSamples ):
	samplePoints = getSamplePoints( exons, halfwinwidth, numSamples )
	scaledProfile = numpy.zeros( len(samplePoints), dtype='float')
	points = profile[samplePoints]
	smoothPoints = signal.savgol_filter( points,147,3)
	scaledProfile = smoothPoints
	
	return scaledProfile

def getExonProfile(halfwinwidth,regions,exons,bamName,fragLength,lower,upper,numSamples):
	bamFile = pysam.Samfile( bamName, 'rb')
	libSize = 0
	for almnt in bamFile: 
		if fragLength==0: # if lib is paired-end
			if not( almnt.is_unmapped) and almnt.is_read1 and almnt.is_proper_pair and abs(almnt.isize)<=upper and abs(almnt.isize)>=lower and almnt.mapq>=20:
				libSize += 1
		else: # if lib is single-end
			if not( almnt.is_unmapped) and almnt.mapq>=20: 
				libSize += 1

	nRegions = len( regions )
	profileMatrix = numpy.zeros( ( nRegions, 7*numSamples ), dtype="d" )
	for row,region in enumerate(regions):
		if bamFile.gettid(region.chrom)==-1: continue
		profile = numpy.zeros( (region.end-region.start), dtype="d" )
		for almnt in bamFile.fetch( region.chrom,max(region.start-200,0),region.end+200 ):
			# If library is PAIRED-END
			if fragLength==0:
				if almnt.is_unmapped or almnt.pos < 0  or not(almnt.is_read1) or not(almnt.is_proper_pair) or abs(almnt.isize)>upper or abs(almnt.isize)<lower or almnt.mapq<20: continue
				if not almnt.is_reverse: # If read is on '+' strand
					almntStart = almnt.pos 
					almntEnd  = almnt.pos + abs(almnt.isize)
				else:                    # If read is on '-' strand
					almntStart = almnt.aend - abs(almnt.isize) 
					almntEnd  = almnt.aend
			# If library is SINGLE-END
			else:
				if almnt.is_unmapped or almnt.pos < 0 or  bamFile.gettid(region.chrom)==-1 or almnt.mapq<20: continue
				if not almnt.is_reverse: # If read is on '+' strand
					almntStart = almnt.pos 
					almntEnd  = almnt.pos + fragLength
				else:                    # If read is on '-' strand
					almntStart = almnt.aend - fragLength 
					almntEnd  = almnt.aend
			
			if almntStart > region.end or almntEnd < region.start: continue
			if region.strand == "+":
				start_in_window = almntStart - region.start
				end_in_window   = almntEnd   - region.start
			else:
				start_in_window = -almntEnd   + region.end
				end_in_window   = -almntStart + region.end 
			start_in_window = max( start_in_window, 0 )
			end_in_window = min( end_in_window, (region.end-region.start) )
			profile[ start_in_window : end_in_window ] += 1
		
		profileScaled = scaleProfile( profile, exons[row,:], halfwinwidth, numSamples )
		profileMatrix[ row, : ] = profileScaled
	bamFile.close()
	# Normalized by library size
	profileMatrix = profileMatrix / libSize * 1e6 
	return profileMatrix

##############################################################	
# getCanonicalNucProfile
##############################################################	
def getNuclCount(nuclFile,cutoff):
	#Assigns an ID to each nucleosome

	nucl = HTSeq.GenomicArrayOfSets('auto',stranded=False)
	nLine=1
	for line in open(nuclFile,'r'):
		fields = line.strip().split('\t')
		if len(fields)==1 or not( fields[1].isdigit()): continue # skip header lines
		# Skip nucleosomes with cutoff below given threshold
		if float(fields[8])<cutoff: continue # using p-value
		#if float(fields[5])<cutoff: continue # using number of supporting reads
		chrom = fields[0]
		start = int(fields[1])
		end   = int(fields[2])
		count          = float(fields[10])
		iv = HTSeq.GenomicInterval(chrom,start,end,'.')
		nucl[iv] += str(chrom)+"_"+str(start)+"_"+str(end)+"_"+str(count)
	return nucl
##############################################################	
def getCountsPerNuc(exonList,nucleosomes,nucPos,halfWin,exon_position,labels):
	# Make labels dict
	labelsDict={}
	for idx,position in enumerate(nucPos):
		labelsDict[position]=labels[idx]
	# Initialized nucPositions
	exonPositions, exonCounts = {}, {}
	for nucleosome in labels: 
		exonPositions[nucleosome]=[]
		exonCounts[nucleosome]=[]
	# Find nuc list per exon
	for line in open(exonList,"r"):
		exon_name = line
		line = line.strip().split("\t")
		chrom,start,end,exon_strand = line[0],int(line[1]),int(line[2]),line[3]
		if exon_strand=="+":
			if exon_position=="5p": exon_ref = start
			else:    				exon_ref = end
		else:
			if exon_position=="5p": exon_ref = end
			else:    				exon_ref = start
		window = HTSeq.GenomicInterval( chrom,exon_ref-halfWin,exon_ref+halfWin)
		# Intersect nucleosomes with current exon
		distances, counts = [], []
		for nuc_iv,nuc_name in nucleosomes[ window ].steps():
			if len(nuc_name)==0: continue # if no nuc, skip
			if len(nuc_name)>1: 
				print "Overlapping nuc on exon ", exon_name
				for x in nuc_name:
					print nuc_iv,x
			chrom, start, end, count = list(nuc_name)[0].split("_")
			nuc_midpoint = numpy.mean([int(start),int(end)])
			if exon_strand=="+": distance = nuc_midpoint - exon_ref
			if exon_strand=="-": distance = exon_ref - nuc_midpoint
			if abs(distance)>halfWin: continue # skip distance out of range
			distances.append(distance)
			counts.append(float(count))
		# Fetch nucleosome by nominal positions
		#countsSort, distancesSort  = sortNucs(counts,distances) 
		countsSort, distancesSort  = sortNucs2(counts,distances,nucPos,labelsDict) 
		#for nucleosome in nucPos:
		for nucleosome in labels:
			if not( nucleosome in countsSort): continue
			exonPositions[nucleosome] += distancesSort[nucleosome]
			exonCounts[nucleosome] += countsSort[nucleosome]
	return exonCounts, exonPositions
##############################################################	
def sortNucs(counts,distances):
	# Sort nucleosomes names by their distances to exon start
	counts=[y for (x,y) in sorted(zip(distances,counts))]
	distances.sort()
	counts=numpy.array(counts)
	distances=numpy.array(distances)
	# Assign nucl position labels
	n=len(distances) # number of elements 
	n_neg=sum(distances<0) # number of negative elements
	nucLabels=numpy.array( range(-n_neg,n-n_neg) )
	nucLabels[nucLabels>=0] += 1
	
	countsSort = {}
	distancesSort = {}
	for idx,label in enumerate(nucLabels):	
		countsSort[label]    =  counts[idx]
		distancesSort[label] =  distances[idx]
	# Output results
	return countsSort, distancesSort
##############################################################	
def sortNucs2(counts,distances,nucPos,labelsDict):
	# Sort nucleosomes names by their distances to exon start
	counts=[y for (x,y) in sorted(zip(distances,counts))]
	distances.sort()
	counts=numpy.array(counts)
	distances=numpy.array(distances)
	
	countsSort = {}
	distancesSort = {}
	# Start from the positive distances
	for idx,distance in enumerate(distances):
		min_idx = numpy.argmin( abs(nucPos - distance) )
		min_dist = nucPos[min_idx]
		label = labelsDict[min_dist]
		if label in distancesSort:
			distancesSort[ label ] += [ distance ]
			countsSort[ label ]    += [ counts[idx] ]
		else:
			distancesSort[ label ] = [ distance ]
			countsSort[ label ] = [ counts[idx] ]
	return countsSort, distancesSort	

########################################################################
# getNucProfile
########################################################################
def getNucCoverage(halfwinwidth,regions,signalFile,controlFile,pvalue,expValues):

	nuc_signals = HTSeq.GenomicArrayOfSets('auto',stranded=False)
	with open(signalFile,'r') as sFile, open(controlFile,'r') as cFile:
		for sLine,nLine in itertools.izip(sFile,cFile):
			if float(cLine[8])<-numpy.log10(pvalue): continue
			sLine = sLine.strip().split("\t")
			nLine = nLine.strip().split("\t")
			chrom, start, end = sLine[0], int(sLine[1]),int(sLine[2])
			nuc_iv = HTSeq.GenomicInterval(chrom, start,end)
			x,n = float(sLine[10]),float(cLine[10])
			r = x/expValues[n]
			nuc_signals[ nuc_iv ] += r

	nRegions = len( regions )
	profileMatrix = numpy.zeros( ( nRegions, 2*halfwinwidth ), dtype="d" )
	for row,region in enumerate(regions):
		for nuc_iv,nuc_signal in nuc_signals[region].steps():

			if region.strand == "+":
				start_in_window = nuc_iv.start - region.start
				end_in_window   = nuc_iv.end   - region.start
			else:
				start_in_window = -nuc_iv.end   + region.end
				end_in_window   = -nuc_iv.start + region.end

			start_in_window = max( start_in_window, 0 )
			end_in_window = min( end_in_window, 2*halfwinwidth )
			profileMatrix[ row, start_in_window : end_in_window ] += nuc_signal
		
	return profileMatrix
