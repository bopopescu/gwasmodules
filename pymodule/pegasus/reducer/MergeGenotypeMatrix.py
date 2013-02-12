#!/usr/bin/env python
"""
Examples:
	#testing merge three identical genotype files
	%s -o /tmp/ccc.tsv /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
	
	%s 
	
Description:
	2011-7-12
		This program merges all files with the same header into one while retaining the header.
	2012.7.31 the input file could be gzipped as well.

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions

class MergeGenotypeMatrix(object):
	__doc__ = __doc__
	option_default_dict = {('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('noHeader', 0, int): [0, 'n', 0, 'all input has no header'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		header = None
		outf = open(self.outputFname, 'w')
		for inputFname in self.inputFnameLs:
			if not os.path.isfile(inputFname):
				continue
			suffix = os.path.splitext(inputFname)[1]
			if suffix=='.gz':
				import gzip
				inf = gzip.open(inputFname, 'r')
			else:
				inf = open(inputFname, 'r')
			if self.noHeader==0:	#in the case that every input has a common header
				if not header:	#2012.7.26 bugfix: empty file will return an empty string, which "is not None". 
					try:
						header = inf.readline()
						outf.write(header)
					except:	#in case something wrong (i.e. file is empty)
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
						import traceback
						traceback.print_exc()
						print sys.exc_info()
				else:
					#skip the header for other input files
					try:
						inf.readline()
					except:	#in case something wrong (i.e. file is empty)
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
						import traceback
						traceback.print_exc()
						print sys.exc_info()
			for line in inf:
				outf.write(line)
		

if __name__ == '__main__':
	main_class = MergeGenotypeMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
