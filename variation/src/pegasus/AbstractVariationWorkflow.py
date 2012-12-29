#!/usr/bin/env python
"""
2012.2.15
	an abstract class for pegasus workflows that work on arabidopsis variation data
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import re
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from pymodule.pegasus.AbstractWorkflow import AbstractWorkflow
from variation.src import Stock_250kDB

class AbstractVariationWorkflow(AbstractWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractWorkflow.option_default_dict.copy()
	
	option_default_dict.update({
						('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['banyan', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'stock_250k database name', ],\
						('genome_dbname', 1, ): ['genome', '', 1, 'genome database name', ],\
						('schema', 0, ): ['', '', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, '', 1, 'database password', ],\
						('port', 0, ):[None, '', 1, 'database port number'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_alignment and/or individual_alignment.path'],\
						
						('data_dir', 0, ):[None, 't', 1, 'The top folder of db-affiliated storage. If not given, use db.data_dir.'],\
						('tabixPath', 1, ): ["%s/bin/tabix", '', 1, 'path to the tabix binary', ],\
						("variationSrcPath", 1, ): ["%s/script/variation/src", 'S', 1, 'variation source code folder'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.2.15
		"""
		AbstractWorkflow.__init__(self, **keywords)
		
		self.variationSrcPath = self.insertHomePath(self.variationSrcPath, self.home_path)
		self.tabixPath =  self.insertHomePath(self.tabixPath, self.home_path)
		
		self.chr_pattern = re.compile(r'(\w+\d+).*')
		self.contig_id_pattern = re.compile(r'Contig(\d+).*')
		
		self.connectDB()
		
		if not self.data_dir:
			self.data_dir = self.db.data_dir
	
	def connectDB(self):
		"""
		2012.6.5
			overwrite the parent class
		"""
		self.db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema,\
				port=self.port)
		self.db_250k.setup(create_tables=False)
	
	
	def registerExecutables(self, workflow=None):
		"""
		2012.2.15
		"""
		AbstractWorkflow.registerExecutables(self)
		
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		variationSrcPath = self.variationSrcPath
		vervetSrcPath = self.vervetSrcPath
		
		#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		executableClusterSizeMultiplierList = []
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		