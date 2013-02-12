#!/usr/bin/env python
"""
Examples:
	#testing merge three identical genotype files
	%s -o /tmp/ccc.tsv /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
	
	%s --exitNonZeroIfAnyInputFileInexistent  -o /tmp/ccc.tsv /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
	
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

import copy
from pymodule import ProcessOptions
from AbstractReducer import AbstractReducer

class MergeGenotypeMatrix(AbstractReducer):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractReducer.option_default_dict)
	option_default_dict.update({
					})

	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2011-7-12
		"""
		AbstractReducer.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		header = None
		outf = open(self.outputFname, 'w')
		for inputFname in self.inputFnameLs:
			sys.stderr.write("File %s ... "%(inputFname))
			if not os.path.isfile(inputFname):
				if self.exitNonZeroIfAnyInputFileInexistent:
					sys.stderr.write(" doesn't exist. Exit 3.\n")
					sys.exit(3)
				else:
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
			sys.stderr.write(".\n")

if __name__ == '__main__':
	main_class = MergeGenotypeMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
