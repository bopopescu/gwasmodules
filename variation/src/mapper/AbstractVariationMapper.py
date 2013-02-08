#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 

Description:
	2012.3.8
		abstract mapper for variation mappers.

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import AbstractMapper, AbstractDBInteractingJob
from variation.src import Stock_250kDB

class AbstractVariationMapper(AbstractDBInteractingJob):
	__doc__ = __doc__
	option_default_dict = AbstractDBInteractingJob.option_default_dict.copy()
	option_default_dict.update({
							('drivername', 1,):['postgresql', '', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['vervetdb', '', 1, 'stock_250k database name', ],\
							('schema', 0, ): ['stock_250k', '', 1, 'database schema name', ],\
							('min_MAF', 0, float): [0.1, '', 1, 'minimum Minor Allele Frequency.'],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractDBInteractingJob.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
	def connectDB(self):
		"""
		2012.11.18
		"""
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, db_user=self.db_user, db_passwd=self.db_passwd, \
									hostname=self.hostname, dbname=self.dbname, schema=self.schema, port=self.port)
		db_250k.setup(create_tables=False)
		self.db_250k = db_250k
	
	def run(self):
		pass
	
if __name__ == '__main__':
	main_class = AbstractVariationMapper
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()