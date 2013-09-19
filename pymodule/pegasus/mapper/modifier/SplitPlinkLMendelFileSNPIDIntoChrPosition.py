#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i input.tped.gz -o /tmp/output.tped

Description:
	2013.08 program that splits SNP-ID in plink .lmendel output into chromosome position , two columns.
		all input files could be gzipped or not.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, MatrixFile
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class SplitPlinkLMendelFileSNPIDIntoChrPosition(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('run_type', 1, int): [1, 'y', 1, 'which run type. 1:'],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
	
	def processRow(self, row):
		"""
		2012.8.9
		"""
		oldChromosome, snp_id, count = row[:3]
		snp_id_split = snp_id.split('_')
		if len(snp_id_split)>2:
			chromosome, start, stop = snp_id_split[:3]
		else:
			chromosome, start = snp_id_split[:2]
			stop = ""
		new_row = [snp_id, oldChromosome, chromosome, start, stop, count]
		return new_row
	
	
	def run(self):
		"""
		2013.07.24
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#inf = utils.openGzipFile(self.inputFname)
		reader = MatrixFile(inputFname=self.inputFname)
		reader.constructColName2IndexFromHeader()
		writer = MatrixFile(inputFname=self.outputFname, openMode='w', delimiter='\t')
		header = ["SNPID", "oldChromosome", "Chromosome", "Start", "Stop", "N"]
		writer.writeHeader(header)
		
		counter = 0
		for row in reader:
			new_row = self.processRow(row)
			writer.writerow(new_row)
			counter += 1
		sys.stderr.write("%s lines processed.\n"%(counter))
		
		del reader
		del writer
	
	
if __name__ == '__main__':
	main_class = SplitPlinkLMendelFileSNPIDIntoChrPosition
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()