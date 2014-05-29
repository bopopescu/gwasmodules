#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /tmp/input.fasta.gz -a 2 -O /tmp/output -f .fasta 

Description:
	2012.5.24
		split a big fasta input file (gzipped or not).
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from Bio import SeqIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class SplitFastaFile(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.update({
							('noOfSequences', 1, int): [1000, 'l', 1, 'number of sequences in each split file', ],\
							('suffixLength', 1, int): [3, 'a', 1, 'length of suffix that is used to distinguish all split files. i.e. 001,002,003', ],\
							('filenameSuffix', 0, ): ['', 'f', 1, 'the suffix attached to the final split filename. i.e. .fasta, .txt', ],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2012.5.23
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def splitFastaFile(self, inputFname=None, outputFnamePrefix=None, noOfSequences=1000, suffixLength=3, filenameSuffix=""):
		"""
		2012.5.24
		"""
		sys.stderr.write("Splitting fasta file %s ..."%(inputFname))
		inf = utils.openGzipFile(inputFname)
		counter = 0 
		real_counter = 0
		outputFname = utils.comeUpSplitFilename(outputFnamePrefix=outputFnamePrefix, suffixLength=suffixLength, fileOrder=real_counter,\
											filenameSuffix=filenameSuffix)
		outputHandle = open(outputFname, 'w')
		for seq_record in SeqIO.parse(inf, "fasta"):
			counter += 1
			SeqIO.write([seq_record], outputHandle, "fasta")
			if counter%noOfSequences==0:
				outputHandle.close()
				real_counter += 1
				outputFname = utils.comeUpSplitFilename(outputFnamePrefix=outputFnamePrefix, suffixLength=suffixLength, fileOrder=real_counter,\
											filenameSuffix=filenameSuffix)
				outputHandle = open(outputFname, 'w')
		#close the last handle
		outputHandle.close()
		sys.stderr.write(" into %s files.\n"%(real_counter+1))	#real_counter starts from 0
		
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.splitFastaFile(inputFname=self.inputFname, outputFnamePrefix=self.outputFnamePrefix, noOfSequences=self.noOfSequences, \
						suffixLength=self.suffixLength, filenameSuffix=self.filenameSuffix)
		
if __name__ == '__main__':
	main_class = SplitFastaFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()