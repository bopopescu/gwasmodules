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
	option_default_dict.update(AbstractWorkflow.db_option_dict.copy())
	
	option_default_dict.update({
						('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						
						('schema', 0, ): ['stock_250k', '', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, '', 1, 'database password', ],\
						('genome_dbname', 1, ): ['genome', '', 1, 'genome database name', ],\
						('genome_schema', 1, ): ['genome', '', 1, 'genome database schema name', ],\
						
						
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
		self.db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, db_user=self.db_user,
				db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname, schema=self.schema,\
				port=self.port)
		self.db_250k.setup(create_tables=False)
	
	def addStock_250kDBJob(self, executable=None, outputFile=None, run_type=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, sshDBTunnel=None, **keywords):
		"""
		2012.11.13 expand all short-arguments to become long ones
		2012.6.5
		"""
		extraArgumentList = []
		if run_type:
			extraArgumentList.append('--run_type %s'%(run_type))
		
		job = self.addGenericDBJob(executable=executable, inputFile=None, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, extraArguments=extraArguments,\
						job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel,\
						objectWithDBArguments=self,\
						**keywords)
		return job
	
	
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
		
		Stock_250kDB = Executable(namespace=namespace, name="Stock_250kDB", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		Stock_250kDB.addPFN(PFN("file://" + os.path.join(self.variationSrcPath, "db/Stock_250kDB.py"), site_handler))
		executableClusterSizeMultiplierList.append((Stock_250kDB, 0))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		