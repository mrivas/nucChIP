import argparse
import sys

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Counts the number of overlapping between two BED files.')
parser.add_argument('-f1',type=str,dest="file1",help="BED enrichment file replicate 1.")
parser.add_argument('-f2',type=str,dest="file2",help="BED enrichment file replicate 2.")
parser.add_argument('-o',type=str,dest="oFile",help="BED enrichment file with the intersection of replicates 1 and 2.")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
file1 = args.file1
file2 = args.file2
oFile = args.oFile

nucl1={}
for line in open(file1,'r'):
	line=line.strip("\n")
	fields=line.split("\t")
	key=fields[0]+"_"+fields[1]+"_"+fields[2]
	nucl1[key]=1

nucl2={}
for line in open(file2,'r'):
	line=line.strip("\n")
	fields=line.split("\t")
	key=fields[0]+"_"+fields[1]+"_"+fields[2]
	nucl2[key]=1
# Counts matches
out=open(oFile,'w')
count=0
for key1 in nucl1:
	if key1 in nucl2:
		count += 1
		fields=key1.split("_")
		print >> out, fields[0]+"\t"+fields[1]+"\t"+fields[2]

out.close()
print file1,file2,count
