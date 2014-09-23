#!/usr/bin/env python
import pysam, HTSeq, numpy, pickle

#########################################################################
# bam2bed.py
#########################################################################
def getAvrLength(bamName):
    # If BAM file is paired-end, outputs the average fragment length
    # This value is then used to extend orphans reads

    bamFile = HTSeq.BAM_Reader(bamName)
    
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
#                print "Pair names are different"
#                print first.read.name, second.read.name
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
def printBED(avrLength,bamName,extension,oFile):
    # Take reads in BAM format and print them in BED format: chrom start end 

    bamFile = HTSeq.BAM_Reader(bamName)
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
#                a.tlen = avrLength
#                a.seq   = "*"
#                a.qual  = "*"
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
#            a.tlen = avrLength
#            a.seq   = "*"
#            a.qual  = "*"
            a.setTag('AL', avrLength, value_type='i', replace=True)
            out.write(a)
    bamFileRef.close()
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

##############################################################
# getProfile
##############################################################
def getRegions(geneList,database,halfwinwidth): #regionsFile,halfwinwidth):
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

def getProfile(halfwinwidth,regions,bamName):
    bamFile = pysam.Samfile( bamName, 'rb')
    libSize = 0
    for almnt in bamFile: 
        if not( almnt.is_unmapped) : libSize += 1

    if type(regions)==dict:
        nucArrayCoverage = {}
        for gene_id in regions:
            regionsList = regions[gene_id]
            nRegions = len( regionsList )
            profileMatrix = numpy.zeros( ( nRegions, 2*halfwinwidth ), dtype="d" )
            for row,region in enumerate(regionsList):
                for almnt in bamFile.fetch( region.chrom,region.start,region.end ):
                    if almnt.is_unmapped or almnt.pos < 0: continue
                    if region.strand == "+":
                        start_in_window = almnt.pos  - region.start
                        end_in_window   = almnt.aend - region.start
                    else:
                        start_in_window = -almnt.aend + region.end
                        end_in_window   = -almnt.pos  + region.end 
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*halfwinwidth )
                    profileMatrix[ row, start_in_window : end_in_window ] += 1
        
            # Normalized by library size
            profileMatrix = profileMatrix / libSize * 1e6 
            meanCoverage  = numpy.mean( profileMatrix, axis=1 )
            nucArrayCoverage[ gene_id ] = list(meanCoverage)
        bamFile.close()
        return nucArrayCoverage
    elif type(regions)==list:
        nRegions = len( regions )
        profileMatrix = numpy.zeros( ( nRegions, 2*halfwinwidth ), dtype="d" )
        for row,region in enumerate(regions):
            for almnt in bamFile.fetch( region.chrom,region.start,region.end ):
                if almnt.is_unmapped or almnt.pos < 0: continue
                if region.strand == "+":
                    start_in_window = almnt.pos  - region.start
                    end_in_window   = almnt.aend - region.start
                else:
                    start_in_window = -almnt.aend + region.end
                    end_in_window   = -almnt.pos  + region.end 
                start_in_window = max( start_in_window, 0 )
                end_in_window = min( end_in_window, 2*halfwinwidth )
                profileMatrix[ row, start_in_window : end_in_window ] += 1
        bamFile.close()
        # Normalized by library size
        profileMatrix = profileMatrix / libSize * 1e6 
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
        fpkm    = line[6]
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
        pos    = int( tss[gene_id][2] )
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
    
    nucRelPos = []
    for line in open( nucRelPosFile, 'r'):
        line = round(float( line.strip() ))
        nucRelPos.append( line )
    
    tss=pickle.load( open(tssData, 'r') )
    nuclArray = {}
    for gene_id in tss:
        chrom  = tss[gene_id][0]
        tssPos    = int( tss[gene_id][2] )
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

    return nuclArray
