#!/usr/bin/env python
"""
2011-7-11
	a class wrapper for SAM/BAM file. It is an extension of pysam.Samfile
	
	http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html#pysam.Samfile
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
from utils import dict_map, importNumericArray, figureOutDelimiter, PassingData
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
