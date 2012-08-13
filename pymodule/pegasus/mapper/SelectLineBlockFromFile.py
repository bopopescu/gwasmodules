#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i input.tped.gz -o /tmp/output.tped

Description:
	2012.7.30
		This program selects lines specified on the commandline and output them.
		The input file could be plain text or gzipped.
		The line number starts from 1 (not 0).
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from AbstractMapper import AbstractMapper

class SelectLineBlockFromFile(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('startLineNumber', 1, int): [1, 's', 1, 'the 1st line of the chosen block ', ],\
						('stopLineNumber', 1, int): [1, 't', 1, 'the last line to be included', ],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		inf = utils.openGzipFile(self.inputFname)
		outf= open(self.outputFname, 'w')
		lineNumber = 0
		real_counter = 0
		for line in inf:
			lineNumber += 1
			if lineNumber>=self.startLineNumber and lineNumber<=self.stopLineNumber:
				outf.write(line);
				real_counter += 1
			elif lineNumber>self.stopLineNumber:	#stop here
				break
			
		inf.close()
		outf.close()
		sys.stderr.write("%s lines chosen.\n"%(real_counter))

if __name__ == '__main__':
	main_class = SelectLineBlockFromFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()