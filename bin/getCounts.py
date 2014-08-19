import argparse, sys, pysam, HTSeq, numpy
import mnaseChIP
###########################################################
def getParser():
	parser = argparse.ArgumentParser(description='Count reads per nucleosome.')
	parser.add_argument('-b',type=str,dest="bFile",help="BAM file. ChIP-seq reads.")
	parser.add_argument('-n',type=str,dest="nFile",help="BED file. Nucleosome data, output of Danpos. Only nucleosomes with p-values < 1e-5 are used")
	parser.add_argument('-l',type=int,dest="lType",help="INT. Library type: if equal to 0, data is assumed to be paired-end and read lengths are estimated from data itself; if different from 0, data is assumed to be single-end and read lengths use this value extend every read.")
	parser.add_argument('-e',type=int,dest="exten",help="INT. Half length of reads extended at their midpoint.")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file.")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser
###########################################################
def main():
	args = getParser().parse_args()
	bamName = args.bFile
	nuclFile = args.nFile
	libType  = args.lType
	extension  = args.exten
	oFile = args.oFile

	####################################################
	# Excecution

	bamFile = HTSeq.BAM_Reader(bamName)

	print "Getting nucleosomes"
	nucl = mnaseChIPgetNucl(nuclFile)
	print "Getting average read length"
	if libType==0: # for paired-end libraries computes fragment length form data itself
		avrLength =  mnaseChIP.getAvrLength(bamFile)
	else:
		avrLength = libType
	print avrLength
	print "Counting reads per nucleosome"
	counts = mnaseChIP.getCounts(avrLength,nucl,bamFile,extension)
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
################################################################
if __name__ == '__main__':
	main()
