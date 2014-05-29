#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -o /Network/Data/250k/tmp-yh/genomeRBDict_tax3702_padding20kb.pickle

Description:
	2012.3.8
		abstract mapper for variation mappers.

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from pymodule import GenomeDB

class PickleGenomeRBDict(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
							('hostname', 1, ): ['uclaOffice', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['vervetdb', 'd', 1, 'genome database name', ],\
							('schema', 0, ): ['genome', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): ['yh', 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, '', 1, 'database port number'],\
							
							('genomeRBDictPickleFname', 1, ): ['', 'o', 1, 'The file to contain pickled genomeRBDict.'],\
							('genePadding', 0, int): [20000, 'x', 1, "the extension around a gene on both sides to allow association between a locus and a gene. Proper distance is LD-decay."],\
							('tax_id', 0, int): [3702, 't', 1, 'Taxonomy ID to get gene position and coordinates.'],\
							
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		self.inputFnameLs = inputFnameLs	#useless
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#need to setup a different db setting
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, \
						schema=self.schema)
		db_genome.setup(create_tables=False)
		
		genomeRBDict = db_genome.dealWithGenomeRBDict(self.genomeRBDictPickleFname, tax_id=self.tax_id, \
													max_distance=self.genePadding, debug=self.debug)


if __name__ == '__main__':
	main_class = PickleGenomeRBDict
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()