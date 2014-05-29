#!/usr/bin/env python
"""
2011-7-11
	a class wrapper for SAM/BAM file. It is an extension of pysam.Samfile
	
	http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html#pysam.Samfile
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule.ProcessOptions import  ProcessOptions
from pymodule.utils import dict_map, importNumericArray, figureOutDelimiter, PassingData
import copy
import pysam

num = importNumericArray()
numpy = num


class BamFile(pysam.Samfile):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): [None, 'o', 1, 'a bam input file.'],\
						('openMode', 1, ): ['rb', '', 1, 'rb: bam file. r: sam file.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, inputFname, openMode, **keywords):
		"""
		2011-7-11
		"""
		pysam.Samfile.__init__(self, inputFname, openMode, **keywords)
		self.inputFname = inputFname
		self.openMode = openMode
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		"""
	def traverseBamByRead(self, processor=None):
		"""
		2011-7-10
			add samfile to param_obj
		2011-2-8
			a traverser used by other functions
		"""
		self.seek(0)
		it = self.fetch()
		counter = 0
		real_counter = 0
		qname2count = {}
		param_obj = PassingData(real_counter=real_counter, counter=counter, qname2count=qname2count, samfile=self)
		for read in it:
			counter += 1
			exitCode = processor.run(read, param_obj=param_obj)
				
			if counter%10000==0:
				sys.stderr.write("%s\t%s\t\t%s"%('\x08'*80, param_obj.counter, param_obj.real_counter))
			if exitCode:	#2011-7-8
				break
		processor.qname2count = param_obj.qname2count	#2011-2-9 pass it to the processor
		max_redundant_read_count = max(param_obj.qname2count.values())
		sys.stderr.write("\n %s unique reads among %s mapped reads, max redundant read count=%s. Done.\n"%\
						(len(param_obj.qname2count), param_obj.real_counter, max_redundant_read_count))

from pysam import AlignedRead
class YHAlignedRead(object):
	def __init__(self, alignedRead=None):
		self.read = alignedRead
		#AlignedRead(self, *keywordList, **keywords)
		
		cigar_code2no_of_bases = {}
		for cigar_tuple in self.read.cigar: #presented as a list of tuples (operation,length).
			#For example, the tuple [ (0,3), (1,5), (0,2) ] refers to an alignment with 3 matches, 5 insertions and another 2 matches.
			# CIGAR operation MIDNSHP=X =>012345678
			cigar_code = cigar_tuple[0]
			if cigar_code not in cigar_code2no_of_bases:
				cigar_code2no_of_bases[cigar_code] = 0
			cigar_code2no_of_bases[cigar_code] += cigar_tuple[1]
		self.no_of_matches = cigar_code2no_of_bases.get(0, 0)	#0 represents match.
		self.match_fraction = self.no_of_matches/float(self.read.alen)
		self.no_of_aligned_bases = self.read.alen
		self.no_of_insertions = cigar_code2no_of_bases.get(1, 0)
		self.no_of_deletions = cigar_code2no_of_bases.get(2, 0)
		self.no_of_mismatches = self.read.alen - self.no_of_matches
		
		if 8 in cigar_code2no_of_bases:
			if cigar_code2no_of_bases[8]!=self.no_of_mismatches:
				sys.stderr.write("Error: No of mismatches in cigar (%s) doesn't match the number (%s) calculated.\n"%\
										(cigar_code2no_of_bases[8], self.no_of_mismatches))
				import pdb
				pdb.set_trace()
		self.cigar_code2no_of_bases = cigar_code2no_of_bases
	
	def getMatchFraction(self):
		"""
		2012.10.14
		"""
		return self.match_fraction
	
	def getNoOfMatches(self):
		"""
		2012.10.14
		"""
		return self.no_of_matches
	
	def getNoOfMismatches(self):
		"""
		2012.10.14
		"""
		return self.no_of_mismatches
	
	def getNoOfInsertions(self):
		"""
		2012.10.14
		"""
		return self.no_of_insertions
	
	def getNoOfDeletions(self):
		"""
		2012.10.14
		"""
		return self.no_of_deletions
	
	def getNoOfIndels(self):
		"""
		2012.10.14
		"""
		return self.no_of_deletions + self.no_of_insertions