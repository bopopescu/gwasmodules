#!/usr/bin/env python2.5
"""
Usage: SnpsStats.py [OPTIONS] -o OUTPUT_R_FILE INPUT_FILE

Option:

        -o ...,	output file
	-d ..., --delim=...         default is \", \"      
        -m ..., --missingval=...    default is \"NA\"
	-a ..., --withArrayId=...   0 for no array ID info (default), 1 if file has array ID info, 2 if comparison file also.
        --comparisonFile=...
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	FilterSnps.py --maxMissing=0.5 -o /tmp/2010_filtered.csv 2010.csv
	
Description:
	Filter a csv formatted file with respect to various criteria.  
	Note that this script only removes SNPs not accessions. 

	Requires MySQLdb to be installed, as well as util.py, rfun.py, snpsdata.py and dataParsers.py.
"""

import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["delim=", "missingval=", "withArrayId=", "comparisonFile=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:d:m:a:brh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	inputFile = args[0]
	output_fname = None
	delim = ", "
	missingVal = "NA"
	comparisonFile = None
	debug = None
	report = None
	help = 0
	withArrayIds = 0

	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-a","--withArrayId"):
			withArrayIds = int(arg)
		elif opt in ("--comparisonFile"):
			comparisonFile = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
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

	import dataParsers
        import snpsdata
        snpsds = dataParsers.parseCSVData(inputFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	

	#Calculating Error rates
	#if comparisonFile:
	#	snpsds2 = dataParsers.parseCSVData(comparisonFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
	#	for i in range(0,len(snpsds)):
                        #Compare ... and record relevant information...
                        #snpsds[i].compare filterBadSnps(snpsds2[i],maxError)
        #            pass

	#Calculating NA rates..
	print "Calculating NA rates"
	snpsNARates = []
	for i in range(0,len(snpsds)):
		snpsNARates += snpsds[i].getSnpsNArates()
	import util
	rstr = ""
	rstr += "snpsNARates <- c("+",".join(util.valListToStrList(snpsNARates))+")\n"
	rstr += 'hist(snpsNARates, xlab="NA rates", ylab="SNP frequency", breaks=60)'
	
	f = open(output_fname,"w")
	f.write(rstr)
	f.close()


if __name__ == '__main__':
	_run_()



