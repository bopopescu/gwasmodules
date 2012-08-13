#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i input.tped.gz -o /tmp/output.tped

Description:
	2012.7.20
		input file could be gzipped or not.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, Genome
from AbstractMapper import AbstractMapper
import h5py

class ModifyTPED(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('run_type', 1, int): [1, 'y', 1, 'which run type. 1: modify snp_id (2nd-column) = chr_phyiscalPosition,\
						2: change chr (1st column) to X (chromosome X, for sex check by plink) , snp_id=X_physicalPosition.', ],\
						('positionStartBase', 0, int): [0, 'p', 1, 'the number to be added to position of every SNP'],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
	
	def processRow(self, row):
		"""
		2012.8.9
		"""
		chr, snp_id, genetic_distace, physical_distance = row[:4]
		chr = Genome.getContigIDFromFname(chr)
		snp_id = '%s_%s'%(chr, physical_distance)
		new_row = [chr, snp_id, genetic_distace, physical_distance] + row[4:]
		return new_row
	
	def processRow_Chr2X(self, row):
		"""
		2012.8.9
		"""
		chr, snp_id, genetic_distace, physical_distance = row[:4]
		chr = "X"
		snp_id = '%s_%s'%(chr, physical_distance)
		new_row = [chr, snp_id, genetic_distace, physical_distance] + row[4:]
		return new_row
	
	def processRow_addPositionStartBase(self, row):
		"""
		2012.8.9
		"""
		chr, snp_id, genetic_distace, physical_distance = row[:4]
		physical_distance = int(physical_distance) + self.positionStartBase
		#chr = Genome.getContigIDFromFname(chr)
		chr = 1	#everything is on one giant chromosome.
		snp_id = '%s_%s'%(chr, physical_distance)
		new_row = [chr, snp_id, genetic_distace, physical_distance] + row[4:]
		return new_row
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		inf = utils.openGzipFile(self.inputFname)
		reader = csv.reader(inf, delimiter='\t')
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		counter = 0
		for row in reader:
			if self.run_type==2:
				new_row = self.processRow_Chr2X(row)
			elif self.run_type==3:
				new_row = self.processRow_addPositionStartBase(row)
			else:
				new_row = self.processRow(row)
			writer.writerow(new_row)
			counter += 1
			
		del reader
		del writer
		sys.stderr.write("%s lines modified.\n"%(counter))

if __name__ == '__main__':
	main_class = ModifyTPED
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()