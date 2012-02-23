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
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow


class AbstractVariationWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = {
						('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'stock_250k database name', ],\
						('genome_dbname', 1, ): ['genome', 'g', 1, 'genome database name', ],\
						('schema', 0, ): ['', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('port', 0, ):[None, '', 1, 'database port number'],\
						
						('results_directory', 0, ):["", 't', 1, 'The results directory. Default is None. use the one given by db.'],\
						
						("pymodulePath", 1, ): ["%s/script/pymodule", '', 1, 'path to the pymodule folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						('tabixPath', 1, ): ["%s/bin/tabix", '', 1, 'path to the tabix binary', ],\
						("javaPath", 1, ): ["/usr/bin/java", 'J', 1, 'java interpreter binary'],\
						
						("variationSrcPath", 1, ): ["%s/script/variation/src", 'S', 1, 'variation source code folder'],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("site_handler", 1, ): ["condorpool", 'B', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'D', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('clusters_size', 1, int):[15, 'C', 1, 'For short jobs that will be clustered, how many of them should be clustered int one'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
						}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.2.15
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.pymodulePath = self.insertHomePath(self.pymodulePath, self.home_path)
		self.variationSrcPath = self.insertHomePath(self.variationSrcPath, self.home_path)
		self.tabixPath =  self.insertHomePath(self.tabixPath, self.home_path)
		self.vervetSrcPath =  self.insertHomePath(self.vervetSrcPath, self.home_path)
		
		import re
		self.chr_pattern = re.compile(r'(\w+\d+).*')
		self.contig_id_pattern = re.compile(r'Contig(\d+).*')
	
	def registerExecutables(self, workflow):
		"""
		2012.2.15
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		variationSrcPath = self.variationSrcPath
		vervetSrcPath = self.vervetSrcPath
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "shell/mkdirWrap.sh"), site_handler))
		mkdirWrap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%workflow.clusters_size))
		workflow.addExecutable(mkdirWrap)
		workflow.mkdirWrap = mkdirWrap
		
		#mv to rename files and move them
		mv = Executable(namespace=namespace, name="mv", version=version, os=operatingSystem, arch=architecture, installed=True)
		mv.addPFN(PFN("file://" + "/bin/mv", site_handler))
		mv.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mv)
		workflow.mv = mv
		
		#2011-11-28
		mergeSameHeaderTablesIntoOne = Executable(namespace=namespace, name="MergeSameHeaderTablesIntoOne", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		mergeSameHeaderTablesIntoOne.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/MergeSameHeaderTablesIntoOne.py"), site_handler))
		#long arguments will happen, so no clustering.
		#2012.12.8 bug fixed now.
		mergeSameHeaderTablesIntoOne.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mergeSameHeaderTablesIntoOne)
		workflow.mergeSameHeaderTablesIntoOne = mergeSameHeaderTablesIntoOne