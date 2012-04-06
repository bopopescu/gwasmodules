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
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from variation.src import Stock_250kDB

class AbstractVariationMapper(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
							('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['banyan', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'stock_250k database name', ],\
							('schema', 0, ): ['', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, '', 1, 'database port number'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							('min_MAF', 0, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		self.inputFnameLs = inputFnameLs
		
		
		db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname, schema=self.schema, port=self.port)
		db_250k.setup(create_tables=False)
		self.db_250k = db_250k