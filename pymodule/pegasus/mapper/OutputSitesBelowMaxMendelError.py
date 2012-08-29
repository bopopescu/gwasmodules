#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i merged_lmendel.tsv.gz -o merged_lmendelMax6MendelError.tsv -m 6

Description:
	2012.8.28
		This script outputs sites below maxMendelError. The inputFname is output of plink mendel: _lmendel.tsv 
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

class OutputSitesBelowMaxMendelError(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('maxNoOfMendelError', 1, int): [None, 'm', 1, 'the max number number of Mendel errors tolerated.'],\
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
		reader.constructColName2IndexFromHeader()
		noOfMendelErrorColumnIndex = reader.getColIndexGivenColHeader(colHeader='N')
		SNPIDColumnIndex = reader.getColIndexGivenColHeader(colHeader='SNP')
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ['chromosome', 'position', 'noOfMendelErrors']
		writer.writerow(header)
		
		counter = 0
		real_counter = 0
		for row in reader:
			SNPID = row[SNPIDColumnIndex]
			noOfMendelErrors = int(row[noOfMendelErrorColumnIndex])
			if noOfMendelErrors <=self.maxNoOfMendelError:
				chr, pos = SNPID.split('_')
				data_row = [chr, pos, noOfMendelErrors]
				writer.writerow(data_row)
				real_counter += 1
			counter += 1
			
		del reader
		del writer
		sys.stderr.write("%s/%s lines outputted.\n"%(real_counter, counter))

if __name__ == '__main__':
	main_class = OutputSitesBelowMaxMendelError
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()