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
import h5py
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, Genome
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class ModifyTPED(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
						('run_type', 1, int): [1, 'y', 1, 'which run type. 1: modify snp_id (2nd-column) = chr_phyiscalPosition,\
		2: snp_id=chr_physicalPosition (original data), chr (1st column) = X (chromosome X, for sex check by plink), pos += positionStartBase.,\
		3: snp_id=chr_physicalPosition (original data), chr (1st column) = newChr, pos += positionStartBase', ],\
						('newChr', 0, ): ["1", 'n', 1, 'the new chromosome for the TPED file (replace the old one)'],\
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
		#chr = Genome.getContigIDFromFname(chr)	# 2012.8.16 getting rid of the string part of chromosome ID doesn't help.
			#   non-human chromosome numbers would still be regarded as 0 by plink.
		snp_id = '%s_%s'%(chr, physical_distance)
		new_row = [chr, snp_id, genetic_distace, physical_distance] + row[4:]
		return new_row
	
	def processRow_Chr2X(self, row):
		"""
		2012.8.9
		"""
		chr, snp_id, genetic_distace, physical_distance = row[:4]
		snp_id = '%s_%s'%(chr, physical_distance)	#the snp_id is the original contig & position
		#chr = self.newChr	#new chromosome, new position
		physical_distance = int(physical_distance) + self.positionStartBase
		chr = "X"
		new_row = [chr, snp_id, genetic_distace, physical_distance] + row[4:]
		return new_row
	
	def processRow_addPositionStartBase(self, row):
		"""
		2012.8.9
		"""
		chr, snp_id, genetic_distace, physical_distance = row[:4]
		snp_id = '%s_%s'%(chr, physical_distance)	#the snp_id is the original contig & position
		chr = self.newChr	#new chromosome, new position
		physical_distance = int(physical_distance) + self.positionStartBase
		#chr = Genome.getContigIDFromFname(chr)
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