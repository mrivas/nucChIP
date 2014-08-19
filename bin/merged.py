import argparse
import sys

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Merge the rows of two BED files.')
parser.add_argument('-f1',type=str,dest="f1File",help="BED file. File 1.")
parser.add_argument('-f2',type=str,dest="f2File",help="BED file. File 2.")
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file.")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)
args = parser.parse_args()
file1 = args.f1File
file2 = args.f2File
outFile = args.oFile

################################################################
# Excecution

out=open(outFile,"w")
rows1={}
for line in open(file1,"r"):
	line=line.strip("\n")
	fields=line.split("\t")
	key=fields[0]+"_"+fields[1]+"_"+fields[2]
	rows1[key]=1
	print >>out, fields[0]+"\t"+fields[1]+"\t"+fields[2]

for line in open(file2,"r"):
	line=line.strip("\n")
	fields=line.split("\t")
	key=fields[0]+"_"+fields[1]+"_"+fields[2]
	if key in rows1: continue
	print >>out, fields[0]+"\t"+fields[1]+"\t"+fields[2]
out.close()
