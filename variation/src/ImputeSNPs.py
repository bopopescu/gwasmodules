#!/usr/bin/env python
"""
Usage: InputeSnps.py [OPTIONS] -o OUTPUT_FILE INPUT_FILE

Option:

        -o ...,                     output file
	-w ..., --windowSize=...    NPUTE window size, default is 30. 
	-d ..., --delim=...         default is \", \"      
        -m ..., --missingval=...    default is \"NA\"
	-a ..., --withArrayId=...   0 for no array ID info (default), 1 if file has array ID info.
	--monomorphic               Filter monomorphic SNPs (after imputation).
	-h, --help	show this help

Examples:
	ImputeSnps.py -o /tmp/file.csv 2010.csv
	
Description:
       Uses NPUTE to impute data.  Currently only NPUTE imputation is implemented.

"""
import env

path_NPUTE = env.home_dir+"Projects/NPUTE/"

import sys, getopt, traceback
import snpsdata
import tempfile,os
import dataParsers


def writeAsNputeFile(snpsd,filename):
	decoder = {'A':'A','C':'C','G':'G','T':'T','-':'-','NA':'?'}
	f = open(filename,"w")
	for snp in snpsd.snps:
		for i in range(0,len(snp)):
			#snp[i] = decoder.get(snp[i],'?')
                        if snp[i]=='NA':
                                snp[i]='?'                        
		outStr = ",".join(snp)+"\n"
		f.write(outStr)
	f.close()

def checkNputeFile(filename):
	f = open(filename,"r")
	lines = f.readlines()
	linelen = len(lines[0])
	for line in lines:
		if len(line) != linelen:
			print "An error was found in the NPUTE file"
			return False
		snp  = (line.strip()).split(',')
		for nt in snp:
			if not nt in ['A','C','G','T','-','?']:
				print nt
				print "A nucleotide error was found in the NPUTE file"
				return False
						
	print "The NPUTE file appears good?"
	return True
		
	
def readNputeFile(filename,accessions,positions, arrayIds=None):
	decoder = {'A':'A','C':'C','G':'G','T':'T','-':'-','?':'NA'}
	f = open(filename,'r')
	snps = []
	print len(positions), len(f.readlines())
	f.close()
	f = open(filename,'r')

	for i in range(0, len(positions)):
		line = (f.readline()).strip()
		oldline = line
		line = line.upper()
		snp = line.split(',')
		#oldsnp = snp
		#snp = snp[1:]
		#print snp
		if len(snp) != len(accessions):
			print "line nr.",i,", len(snp)=",len(snp),", snp=",snp
			print "line=",oldline
			#print oldsnp[0]
			print len(accessions),accessions
			raise Exception("accessions are messed up")
		for j in range(0,len(snp)):
			if snp[j]=='':
				print ' i , j :',i,j,", snp=",snp
                        #snp[j] = decoder[snp[j]]
                        if snp[j]=='?':
                                snp[j]='NA'
		snps.append(snp)
	f.close()
	print "Read the NPUTE file.."
	return snpsdata.RawSnpsData(snps=snps,positions=positions,accessions=accessions,arrayIds=arrayIds)
		




def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["monomorphic", "monomorphic", "delim=", "missingval=", "withArrayId=", "windowSize=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:d:w:m:a:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	inputFile = args[0]
	output_fname = None
	delim = ","
	missingVal = "NA"
	monomorphic = False
	help = 0
	withArrayIds = 0
	windowSize = 30

	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-a","--withArrayId"):
			withArrayIds = int(arg)
		elif opt in ("--monomorphic"):
			monomorphic = True
		elif opt in ("--windowSize"):
			windowSize = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if not output_fname:
		output_fname
		if help==0:
			print "Output file missing!!\n"
			print __doc__
		sys.exit(2)

	waid1 = withArrayIds==1 or withArrayIds==2
	waid2 = withArrayIds==2

	(snpsds,chromosomes) = dataParsers.parseCSVData(inputFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1, returnChromosomes=True)
	
	accessions = snpsds[0].accessions
	arrayIds = snpsds[0].arrayIds
	positionsList = []
	tmpFiles = []
	#tempfile.tempdir='/tmp'
	i = 1
	for snpsd in snpsds:
		tmpFile1 = tempfile.mkstemp()
		os.close(tmpFile1[0])
		tmpFile2 = tempfile.mkstemp()
		os.close(tmpFile2[0])
		tmpFiles.append((tmpFile1[1],tmpFile2[1]))
		positionsList.append(snpsd.positions)
		print "Preparing data in",tmpFile1[1]
		writeAsNputeFile(snpsd,tmpFile1[1])
		checkNputeFile(tmpFile1[1])
		del snpsd.snps	
		nputeCmd = "python "+path_NPUTE+"NPUTE.py -m 0 -w "+str(windowSize)+" -i "+str(tmpFile1[1])+" -o "+str(tmpFile2[1])
		print "Imputing chromosome",i
		i += 1
		print nputeCmd
		os.system(nputeCmd)  


	for i in range(0,len(tmpFiles)):
		print "Reading chromosome",i+1
		snpsds[i] = readNputeFile(tmpFiles[i][1],accessions,positionsList[i],arrayIds=arrayIds)
		os.remove(tmpFiles[i][0])
		os.remove(tmpFiles[i][1])

	snpsDataSet = snpsdata.SnpsDataSet(snpsds,[1,2,3,4,5])

        #Filtering monomorphic
	if monomorphic:
		print "Filtering monomorphic SNPs"
		for snpsd in snpsds:
			print "Removed", str(snpsd.filterMonoMorphicSnps()),"Snps"


	snpsDataSet.writeToFile(output_fname, deliminator=delim, missingVal = missingVal, withArrayIds = waid1)


if __name__ == '__main__':
	_run_()


