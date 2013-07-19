#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i input.tped.gz -o /tmp/output.tped --tfamFname ...

Description:
	2013.07.19
		this program reads in the tfam (pedigree) file and expands input tped file to include the individuals that 
			are in the pedigree, but not tped file. Assign them with missing genotype (0).
		This program assumes that the order of individuals in .tfam and .tped is same.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, MatrixFile
from pymodule.pegasus.mapper.filter.ModifyTPED import ModifyTPED

parentClass = ModifyTPED
class AppendExtraPedigreeIndividualsToTPED(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.pop(('run_type', 1, int))
	option_default_dict.pop(('mendelErrorFname', 0, ))
	option_default_dict.pop(('newChr', 0, ))
	option_default_dict.pop(('positionStartBase', 0, int))
	option_default_dict.update({
						})
	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
	
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#inf = utils.openGzipFile(self.inputFname)
		reader = MatrixFile(inputFname=self.inputFname)	#a TPED file
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		counter = 0
		individualID2Index = self.getIndividualID2IndexFromTFAMFile(tfamFname=self.tfamFname)
		noOfIndividuals = len(individualID2Index)
		
		noOfExtraIndividuals = None
		for row in reader:
			#chromosome, snp_id, genetic_distace, physical_distance = row[:4]
			noOfExistingIndividuals = len(row[4:])/2
			noOfExtraIndividuals = noOfIndividuals - noOfExistingIndividuals
			writer.writerow(row+ [0]*2*noOfExtraIndividuals)
			counter += 1
			
		del reader
		del writer
		sys.stderr.write("%s rows (loci) and added %s extra individuals.\n"%(counter, noOfExtraIndividuals))

if __name__ == '__main__':
	main_class = AppendExtraPedigreeIndividualsToTPED
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()