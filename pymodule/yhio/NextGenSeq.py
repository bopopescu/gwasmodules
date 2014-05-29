#!/usr/bin/env python
"""
2011-8-28
	class to store functions related to next-gen sequencing
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import re
from Bio import SeqIO
from pymodule import utils


def getPEInputFiles(input_dir, isPE=True):
	"""
	2011-8-28
		copied from MpiBWA.py
	2011-8-5
		add argument isPE, which flags whether input_dir contains PE or single-end reads
		become a classmethod
	2011-2-7
		for paired-end files, sequence_628BWAAXX_1_1.fastq.gz and sequence_628BWAAXX_1_2.fastq.gz
			are regarded as one pair of two files.
	"""
	sys.stderr.write("Pair input files from %s ..."%input_dir)
	pairedEndPrefix2FileLs = {}
	files = os.listdir(input_dir)
	no_of_fastq_files = 0
	for fname in files:
		fname_prefix, fname_suffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(fname)
		if fname_suffix!='.fastq':		#skip non-fastq files
			continue
		no_of_fastq_files += 1
		if isPE==True:
			pairedEndPrefix = fname_prefix[:-2]
			pairedEndOrder = fname_prefix[-2:]
			
			if pairedEndPrefix not in pairedEndPrefix2FileLs:
				pairedEndPrefix2FileLs[pairedEndPrefix] = ['', '']
			
			if pairedEndOrder=='_1':	#the first file
				pairedEndPrefix2FileLs[pairedEndPrefix][0] = fname
			else:
				pairedEndPrefix2FileLs[pairedEndPrefix][1] = fname
		else:
			pairedEndPrefix2FileLs[fname_prefix] = [fname]	#single End
	no_of_files = len(files)
	no_of_pairedEndPrefix = len(pairedEndPrefix2FileLs)
	if no_of_pairedEndPrefix>0:
		avg_no_of_files_per_prefix = no_of_fastq_files/float(no_of_pairedEndPrefix)
	else:
		avg_no_of_files_per_prefix = 0.0
	sys.stderr.write("%.2f files per one pairedEnd prefix. %s fastq files. %s total files. Done.\n"%\
					(avg_no_of_files_per_prefix, no_of_fastq_files, no_of_files))
	return pairedEndPrefix2FileLs

def isFileNameVCF(inputFname, includeIndelVCF=False):
	"""
	2011-11-11
	"""
	isVCF = False
	if (inputFname[-3:]=='vcf' or inputFname[-6:]=='vcf.gz'):
		isVCF=True
		if not includeIndelVCF and inputFname.find('indel')!=-1:	#exclude indel vcf
			isVCF =False
	return isVCF

def isVCFFileEmpty(inputFname, checkContent=False):
	"""
	2011-11-11
		function to test if the input VCF has any locus.
		empty VCF file could still have headers.
	"""
	import os
	if not os.path.isfile(inputFname):
		return True
	fileSize = os.path.getsize(inputFname)
	if fileSize==0:
		return True
	if checkContent:
		if inputFname[-2:]=='gz':
			import gzip
			inf = gzip.open(inputFname)
		else:
			inf = open(inputFname)
		
		fileIsEmpty = True
		for line in inf:
			if line[0]!='#':
				fileIsEmpty=False
				break
		del inf
		return fileIsEmpty
	else:
		return False

def countNoOfChromosomesBasesInFastQFile(inputFname=None):
	"""
	2013.2.16 add the try...except around the parser
	2013.2.9 count the #chromosomes, #bases of inputFname
	"""
	sys.stderr.write("Counting #chromosomes, #bases of %s ..."%(inputFname))
	no_of_chromosomes = 0
	no_of_bases = 0
	inf = utils.openGzipFile(inputFname)
	try:
		for seq_record in SeqIO.parse(inf, 'fastq'):
			no_of_chromosomes += 1
			no_of_bases += len(seq_record)
	except:
		sys.stderr.write("Except after handling %s chromosomes & %s bases.\n"%(no_of_chromosomes, no_of_bases))
		sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
		import traceback
		traceback.print_exc()
		raise
		
	inf.close()
	sys.stderr.write("%s chromosomes, %s bases\n"%(no_of_chromosomes, no_of_bases))
	return utils.PassingData(no_of_chromosomes=no_of_chromosomes, no_of_bases=no_of_bases)

if __name__ == '__main__':
	pass