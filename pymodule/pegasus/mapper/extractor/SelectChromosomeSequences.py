#!/usr/bin/env python
"""
Examples:
	%s  -i /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta -o /tmp/anything.fastq.gz
		--chromosomeList Contig0,Contig1 --inputFileFormat 1 --outputFileFormat 2
	
	# input & output could be plain or gzipped files.
	%s -i /tmp/input.fasta.gz -o /tmp/output.fasta.gz --chromosomeList Contig0,Contig1 
		--inputFileFormat 1 --outputFileFormat 2

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
							('defaultBasePhredQuality', 1, int): [87, '', 1, "if input format is not fastq and output is fastq. need this base quality.\n\
		Assuming Sanger format for the quality score (87=x).",],\
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
		
	
	def selectSequences(self, inputFname=None, outputFname=None, inputFileFormat='fasta', outputFileFormat='fasta', chromosomeSet=None,\
					defaultBasePhredQuality=87):
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
				if outputFileFormat=='fastq' and 'phred_quality' not in seq_record.letter_annotations:
					#fake quality for fastq output
					seq_record.letter_annotations['phred_quality'] = [defaultBasePhredQuality]*len(seq_record.seq)
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
						outputFileFormat=self.outputFileFormat, chromosomeSet=self.chromosomeSet,\
						defaultBasePhredQuality=self.defaultBasePhredQuality)
		
if __name__ == '__main__':
	main_class = SelectChromosomeSequences
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()