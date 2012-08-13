#!/usr/bin/env python
"""
2012.2.15
	an abstract class for pegasus workflows that work on arabidopsis variation data
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractWorkflow import AbstractWorkflow


class AbstractVariationWorkflow(AbstractWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractWorkflow.option_default_dict.copy()
	
	option_default_dict.update({
						('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['banyan', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'stock_250k database name', ],\
						('genome_dbname', 1, ): ['genome', 'g', 1, 'genome database name', ],\
						('schema', 0, ): ['', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('port', 0, ):[None, '', 1, 'database port number'],\
						
						('results_directory', 0, ):["", 't', 1, 'The results directory. Default is None. use the one given by db.'],\
						
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
		
		import re
		self.chr_pattern = re.compile(r'(\w+\d+).*')
		self.contig_id_pattern = re.compile(r'Contig(\d+).*')
	
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
		
		executableList = []
		
		for executable in executableList:
			executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			self.addExecutable(executable)
			setattr(self, executable.name, executable)