#!/usr/bin/env python
"""
Examples:
	%s
	
	%s -i  LDPrunedMerged.bim -o LDPrunedMerged.tsv 

Description:
	2012.9.13
		This script converts plink's bim file (where markers are stored) into a tab delimited output.
		Input file could be gzipped or not.
		Output has one line header, tab-delimited. could be used by vcftools to keep sites.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, MatrixFile
from AbstractMapper import AbstractMapper

class ConvertPlinkBIM(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('run_type', 1, int): [1, 'y', 1, 'which run_type to run. \n\
		1: convert to tab-delimited, with one line header. two column: chr, position.\n'
		],\
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
		
		reader = MatrixFile(self.inputFname)
		#reader.constructColName2IndexFromHeader()	#no header
		#noOfMendelErrorColumnIndex = reader.getColIndexGivenColHeader(colHeader='N')
		SNPIDColumnIndex = 1
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ['chromosome', 'position']
		writer.writerow(header)
		
		counter = 0
		real_counter = 0
		for row in reader:
			SNPID = row[SNPIDColumnIndex]
			chr, pos = SNPID.split('_')
			data_row = [chr, pos]
			writer.writerow(data_row)
			real_counter += 1
			counter += 1
			
		del reader
		del writer
		sys.stderr.write("%s/%s lines outputted.\n"%(real_counter, counter))

if __name__ == '__main__':
	main_class = ConvertPlinkBIM
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()