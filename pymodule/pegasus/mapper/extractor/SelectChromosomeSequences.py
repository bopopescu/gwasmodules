#!/usr/bin/env python
"""
Examples:
	%s 
	# input & output could be plain or gzipped files.
	%s -i /tmp/input.fasta.gz -o /tmp/output.fasta.gz --chromosomeList Contig0,Contig1 

Description:
	2013.2.15 program selects certain chromosome sequences out of fasta/fastq files.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Bio import SeqIO
from pymodule import ProcessOptions, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class SelectChromosomeSequences(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
							('chromosomeList', 1, ): [None, '', 1, 'coma-separated list of chromosome IDs', ],\
							('inputFileFormat', 1, int): [1, '', 1, "1: fasta; 2: fastq", ],\
							('outputFileFormat', 1, int): [1, '', 1, "1: fasta; 2: fastq", ],\
							('defaultBaseQuality', 1, ): ["z", '', 1, "if input format is not fastq and output is. base quality to be assigned to each base.\n\
		Assuming Sanger format. z is ascii no. 122, corresponding to quality 87 ",],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2012.5.23
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		self.chromosomeList = utils.getListOutOfStr(self.chromosomeList, data_type=str, separator2=None)
		self.chromosomeSet = set(self.chromosomeList)
		
		self.fileFormatDict = {1: 'fasta', 2:'fastq'}
		
		self.inputFileFormat = self.fileFormatDict.get(self.inputFileFormat)
		self.outputFileFormat = self.fileFormatDict.get(self.outputFileFormat)
		
	
	def selectSequences(self, inputFname=None, outputFname=None, inputFileFormat='fasta', outputFileFormat='fasta', chromosomeSet=None):
		"""
		2012.5.24
		"""
		sys.stderr.write("Choosing %s chromosome sequences from %s ..."%(len(chromosomeSet), inputFname))
		inf = utils.openGzipFile(inputFname, 'r')
		counter = 0 
		real_counter = 0
		outputHandle = utils.openGzipFile(outputFname, 'w')
		for seq_record in SeqIO.parse(inf, inputFileFormat):
			counter += 1
			if seq_record.id in chromosomeSet:
				SeqIO.write([seq_record], outputHandle, outputFileFormat)
				real_counter += 1
			elif real_counter==len(chromosomeSet):	#got enough chromosomes
				break
		#close the last handle
		outputHandle.close()
		sys.stderr.write(" %s records chosen into %s.\n"%(real_counter, outputFname))
		
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.selectSequences(inputFname=self.inputFname, outputFname=self.outputFname, inputFileFormat=self.inputFileFormat, \
						outputFileFormat=self.outputFileFormat, chromosomeSet=self.chromosomeSet)
		
if __name__ == '__main__':
	main_class = SelectChromosomeSequences
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()