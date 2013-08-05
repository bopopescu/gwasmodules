#!/usr/bin/env python
"""
Examples:
	%s 
	
	#
	%s 

Description:
	2013.07.31
		abstract walker class for genome files
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, random, numpy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.db import GenomeDB
from pymodule.yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker

parentClass = AbstractMatrixFileWalker

class AbstractGenomeFileWalker(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(parentClass.genome_db_option_dict.copy())
	#option_default_dict.update(parentClass.db_option_dict.copy())
	
	genome_option_dict = {
					('chromosomeHeader', 1, ): ['chromosome', '', 1, 'header of the column that designates chromosome info' ],\
					('positionHeader', 1, ): ['position', '', 1, 'header of the column that designates position of each locus on chromosome' ],\
					
					('tax_id', 0, int): [3702, '', 1, 'taxonomy ID of the organism from which to retrieve the chromosome info', ],\
					('sequence_type_id', 0, int): [1, '', 1, 'sequence_type_id (annot_assembly) of the chromosomes, to retrieve the chromosome info. 1: assembledChromosome, 9: Scaffolds', ],\
					('chrOrder', 0, int): [1, '', 1, 'how to order chromosomes. 1: column genome_order in db; 2: by chromosome size, descending', ],\
					
						}
	option_default_dict.update(genome_option_dict)
	option_default_dict.update({
					
					
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
	
	def _loadGenomeStructureFromDB(self):
		"""
		2013.07.31
			so that both AbstractGenomeFileWalker and PlotGenomeWideData could call it
		"""
		db_genome = GenomeDB.GenomeDatabase(drivername=self.genome_drivername, username=self.genome_db_user,
						db_passwd=self.genome_db_passwd, hostname=self.genome_hostname, dbname=self.genome_dbname, \
						schema=self.genome_schema)
		db_genome.setup(create_tables=False)
		self.db_genome = db_genome
		
		#chrOrder=1 is to order chromosomes alphabetically
		self.oneGenomeData = db_genome.getOneGenomeData(tax_id=self.tax_id, chr_gap=0, chrOrder=self.chrOrder, \
													sequence_type_id=self.sequence_type_id)
		self.chr_id2cumu_start = self.oneGenomeData.chr_id2cumu_start
		#2013.2.18 get centromere locations
		self.chr_id2centromere = self.oneGenomeData.chr_id2centromere
		
	
	def setup(self, **keywords):
		"""
		2013.07.31
			construct an RBTree dictionary map between windows and their data
		"""
		parentClass.setup(self, **keywords)
		self._loadGenomeStructureFromDB()
	

if __name__ == '__main__':
	main_class = AbstractGenomeFileWalker
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()