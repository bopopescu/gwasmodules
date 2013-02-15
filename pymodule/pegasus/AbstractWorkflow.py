#!/usr/bin/env python
"""
2012.5.23
	a common class for pegasus workflows
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData
from pymodule.yhio import NextGenSeq
from pymodule.pegasus import yh_pegasus
from Pegasus.DAX3 import *
from Pegasus import DAX3

class AbstractWorkflow(ADAG):
	__doc__ = __doc__
	db_option_dict = {
					('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
					('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
					('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
					('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
					('db_user', 1, ): [None, 'u', 1, 'database username', ],\
					('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
					('port', 0, ):[None, '', 1, 'database port number'],\
					('commit', 0, ):[None, '', 0, 'commit database transaction if there is db transaction'],\
					("data_dir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
					("local_data_dir", 0, ): ["", 'D', 1, 'this one should contain same files as data_dir but accessible locally.\
							If not given, use the default stored in db (db.data_dir). This argument is used to find all input files available.\n\
							It should be different from data_dir only when you generate a workflow on one computer and execute it on another which has different data_dir.'],\
						
					}
	option_default_dict = {
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("pymodulePath", 1, ): ["%s/script/pymodule", '', 1, 'path to the pymodule folder'],\
						("variationSrcPath", 1, ): ["%s/script/variation/src", '', 1, 'variation source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("javaPath", 1, ): ["/usr/bin/java", 'J', 1, 'path to java interpreter binary'],\
						("plinkPath", 1, ): ["%s/bin/plink", '', 1, 'path to the plink binary, http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('clusters_size', 1, int):[30, 'C', 1, 'For short jobs that will be clustered, how many of them should be clustered int one'],\
						('pegasusFolderName', 0, ): ['folder', 'F', 1, 'the folder relative to pegasus workflow root to contain input & output.\
								It will be created during the pegasus staging process. It is useful to separate multiple workflows.\
								If empty, everything is in the pegasus root.', ],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('needSSHDBTunnel', 0, int):[0, 'H', 0, 'DB-interacting jobs need a ssh tunnel (running on cluster behind firewall).'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
						}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	
	pathToInsertHomePathList = ['javaPath', 'pymodulePath', 'vervetSrcPath', 'plinkPath', 'variationSrcPath']

	def __init__(self,  **keywords):
		"""
		2012.5.23
		"""
		# Create a abstract dag
		ADAG.__init__(self, "myworkflow")
		"""
		2012.8.29 methods of ADAG
		>>> dir(a)
		['__doc__', '__init__', '__module__', '__str__', '__unicode__', 'addDAG', 'addDAX', 'addDependency', 
		'addExecutable', 'addFile', 'addInvoke', 'addJob', 'addTransformation', 'clearDependencies', 
		'clearExecutables', 'clearFiles', 'clearInvokes', 'clearJobs', 'clearTransformations', 'count', 
		'dependencies', 'depends', 'executables', 'files', 'getJob', 'hasDependency', 'hasExecutable', 
		'hasFile', 'hasInvoke', 'hasJob', 'hasTransformation', 'index', 'invocations', 'invoke', 'jobs',
		 'name', 'nextJobID', 'removeDependency', 'removeExecutable', 'removeFile', 'removeInvoke', 
		 'removeJob', 'removeTransformation', 'sequence', 'toXML', 'transformations', 'writeXML', 'writeXMLFile']
		
		"""
		
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		#change the workflow name to reflect the output filename
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		self.name = workflowName
		
		for pathName in self.pathToInsertHomePathList:
			absPath = self.insertHomePath(getattr(self, pathName, None), self.home_path)
			if absPath:
				setattr(self, pathName, absPath)
			else:
				sys.stderr.write("Warning: %s has empty absolute path. Skip.\n"%(pathName))
			
		#self.pymodulePath = self.insertHomePath(self.pymodulePath, self.home_path)
		#self.vervetSrcPath =  self.insertHomePath(self.vervetSrcPath, self.home_path)
		
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		self.architecture = "x86_64"
		self.operatingSystem = "linux"
		self.namespace = "workflow"
		self.version="1.0"
		
		self.commandline = ' '.join(sys.argv)
		
		#2012.9.25 global counter
		self.no_of_jobs = 0
		
		if hasattr(self, 'contigMaxRankBySize') and hasattr(self, 'contigMinRankBySize'):
			#2013.2.6 non-public schema dbs should be connected before the main vervetdb or other db (schema=public) is connected.
			self.chr2size = self.getTopNumberOfContigs(contigMaxRankBySize=self.contigMaxRankBySize, \
													contigMinRankBySize=self.contigMinRankBySize)
		else:
			self.chr2size = {}
	
		self.connectDB()
		self.extra__init__()
	
	def extra__init__(self):
		"""
		2013.2.14
		"""
		pass
	
	def writeXML(self, out):
		"""
		2012.8.29
			copied from /usr/lib/pegasus/python/Pegasus/DAX3.py because i want to output comment.
			overwrite its parent. ADAG.writeXML()
			Write the ADAG as XML to a stream
		"""
		import datetime, pwd, os, sys
		
		# Preamble
		out.write('<?xml version="1.0" encoding="UTF-8"?>\n')
		
		# Metadata
		out.write('<!-- generated: %s -->\n' % datetime.datetime.now())
		out.write('<!-- generated by: %s -->\n' % pwd.getpwuid(os.getuid())[0])
		out.write('<!-- generator: python -->\n')
		out.write('<!-- commandline: %s -->\n'%(self.commandline.replace("--", "~~")))	#2012.9.4 -- is not allowed in xml comment.
		
		# Open tag
		out.write('<adag xmlns="%s" ' % DAX3.SCHEMA_NAMESPACE)
		out.write('xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ')
		out.write('xsi:schemaLocation="%s %s" ' % (DAX3.SCHEMA_NAMESPACE, DAX3.SCHEMA_LOCATION))
		out.write('version="%s" ' % DAX3.SCHEMA_VERSION)
		out.write('name="%s"' % self.name)
		if self.count: out.write(' count="%d"' % self.count)
		if self.index: out.write(' index="%d"' % self.index)
		out.write('>\n')
		
		# Invocations
		for i in self.invocations:
			out.write('\t')
			i.toXML().write(stream=out, level=1)
			out.write('\n')
		
		# Files
		for f in self.files:
			out.write('\t')
			f.toXML().write(stream=out, level=1)
			out.write('\n')
		
		# Executables
		for e in self.executables:
			out.write('\t')
			e.toXML().write(stream=out, level=1)
			out.write('\n')
		
		# Transformations
		for t in self.transformations:
			out.write('\t')
			t.toXML().write(stream=out, level=1)
			out.write('\n')
		
		# Jobs
		keys = self.jobs.keys()
		keys.sort()
		for job_id in keys:
			job = self.jobs[job_id]
			out.write('\t')
			job.toXML().write(stream=out, level=1)
			out.write('\n')
		
		# Dependencies
		# Since we store dependencies as tuples, but we need to print them as nested elements
		# we first build a map of all the children that maps child -> [(parent,label),...]
		children = {}
		for dep in self.dependencies:
			if not dep.child in children:
				children[dep.child] = []
			children[dep.child].append((dep.parent, dep.edge_label))
		
		# Now output all the xml in sorted order by child, then parent
		keys = children.keys()
		keys.sort()
		for child in keys:
			out.write('\t')
			c = DAX3.Element("child",[("ref",child)])
			parents = children[child]
			parents.sort()
			for parent, edge_label in parents:
				p = DAX3.Element("parent",[
					("ref", parent),
					("edge-label", edge_label)
				])
				c.element(p)
			c.write(stream=out, level=1)
			out.write('\n')
		
		# Close tag
		out.write('</adag>\n')
	
	def constructOneExecutableObject(self, path=None, name=None):
		"""
		2013.2.7
		"""
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		site_handler = self.site_handler
		
		executable = Executable(namespace=namespace, name=name, version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		executable.addPFN(PFN("file://" + os.path.expanduser(path), site_handler))
		return executable
	
	def getTopNumberOfContigs(self, **keywords):
		"""
		2013.2.6 placeholder
		"""
		pass
	
	def connectDB(self):
		"""
		2013.1.25 placeholder, to establish db connection
		"""
		self.db = None
	
	def processListArguments(self, listArgumentName_data_type_ls=None, emptyContent=[]):
		"""
		2012.10.5
			moved to ProcessOptions
		2012.8.15
		"""
		return ProcessOptions.processListArguments(listArgumentName_data_type_ls=listArgumentName_data_type_ls,\
									emptyContent=emptyContent, class_to_have_attr=self)
	
	def initiateWorkflow(self, workflowName=None):
		"""
		2012.5.23
			AbstractWorkflow is now a derivative of ADAG.
		2011-11-22
		"""
		"""
		# Create a abstract dag
		workflow = ADAG(workflowName)
		workflow.site_handler = self.site_handler
		workflow.input_site_handler = self.input_site_handler
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		workflow.architecture = "x86_64"
		workflow.operatingSystem = "linux"
		workflow.namespace = "workflow"
		workflow.version="1.0"
		#clusters_size controls how many jobs will be aggregated as a single job.
		workflow.clusters_size = self.clusters_size
		"""
		return self
	
	def insertHomePath(self, inputPath, home_path):
		"""
		2013.1.4 inputPath could be None
		2012.5.23 copied from AbstractNGSWorkflow
		2012.1.9
		"""
		if inputPath:
			if inputPath.find('%s')!=-1:
				inputPath = inputPath%home_path
		else:
			inputPath = None
		return inputPath
	
	def registerJars(self, workflow=None):
		"""
		2012.5.23
			register jars to be used in the worflow
		"""
		pass
	
	def registerCustomJars(self, workflow=None):
		"""
		2012.5.23
		"""
		pass
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012-8.15
		"""
		pass
	
	def registerCommonExecutables(self, workflow=None):
		"""
		2012.7.25
			add noDefaultClustersSizeExecutableList
		2011-11-22
		"""
		self.registerExecutables(workflow=workflow)

	def registerExecutables(self, workflow=None):
		"""
		2012.7.4
			added cp
		2012.1.9 a symlink to registerCommonExecutables()
		"""
		if not workflow:
			workflow = self
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#executableList = []
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		#noClusteringExecutableSet = set()	#2012.8.2 you don't want to cluster for some jobs.
		
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "shell/mkdirWrap.sh"), site_handler))
		executableClusterSizeMultiplierList.append((mkdirWrap, 1))
		
		#mv to rename files and move them
		mv = Executable(namespace=namespace, name="mv", version=version, os=operatingSystem, arch=architecture, installed=True)
		mv.addPFN(PFN("file://" + "/bin/mv", site_handler))
		executableClusterSizeMultiplierList.append((mv, 1))
		
		#the copy command
		cp = Executable(namespace=namespace, name="cp", version=version, os=operatingSystem, arch=architecture, installed=True)
		cp.addPFN(PFN("file://" + "/bin/cp", site_handler))
		executableClusterSizeMultiplierList.append((cp, 1))
		
		gzip = Executable(namespace=namespace, name="gzip", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		gzip.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "shell/gzip.sh"), site_handler))
		executableClusterSizeMultiplierList.append((gzip, 1))
		
		SelectLineBlockFromFile = Executable(namespace=namespace, name="SelectLineBlockFromFile", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		SelectLineBlockFromFile.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/filter/SelectLineBlockFromFile.py"), \
										site_handler))
		executableClusterSizeMultiplierList.append((SelectLineBlockFromFile, 1))
		
		AbstractPlot =  Executable(namespace=namespace, name="AbstractPlot", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		AbstractPlot.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "plot/AbstractPlot.py"), site_handler))
		executableClusterSizeMultiplierList.append((AbstractPlot, 0))
		
		PlotLD = Executable(namespace=namespace, name="PlotLD", version=version, os=operatingSystem, arch=architecture, installed=True)
		PlotLD.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "plot/PlotLD.py"), site_handler))
		executableClusterSizeMultiplierList.append((PlotLD, 0))
		
		PlotYAsBar = Executable(namespace=namespace, name="PlotYAsBar", version=version, os=operatingSystem, arch=architecture, installed=True)
		PlotYAsBar.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "plot/PlotYAsBar.py"), site_handler))
		executableClusterSizeMultiplierList.append((PlotYAsBar, 0))
		
		DrawHistogram = Executable(namespace=namespace, name="DrawHistogram", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		DrawHistogram.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "plot/DrawHistogram.py"), site_handler))
		executableClusterSizeMultiplierList.append((DrawHistogram, 0))
		
		#2012.8.13 SelectRowsFromMatrix is a derivative of AbstractMatrixFileWalker, so use addAbstractMatrixFileWalkerJob()
		SelectRowsFromMatrix = Executable(namespace=namespace, name="SelectRowsFromMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		SelectRowsFromMatrix.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/filter/SelectRowsFromMatrix.py"), site_handler))
		executableClusterSizeMultiplierList.append((SelectRowsFromMatrix, 1))
		
		#2012.8.15 ancestor of SelectRowsFromMatrix, 
		AbstractMatrixFileWalker  = Executable(namespace=namespace, name="AbstractMatrixFileWalker", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		AbstractMatrixFileWalker.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "yhio/AbstractMatrixFileWalker.py"), site_handler))
		executableClusterSizeMultiplierList.append((AbstractMatrixFileWalker, 1))
		
		#2012.8.13
		OutputVCFSiteGap = Executable(namespace=namespace, name="OutputVCFSiteGap", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		OutputVCFSiteGap.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/computer/OutputVCFSiteGap.py"), site_handler))
		executableClusterSizeMultiplierList.append((OutputVCFSiteGap, 1))
		
		java = Executable(namespace=namespace, name="java", version=version, os=operatingSystem, arch=architecture, installed=True)
		java.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((java, 1))
		
		plink =  Executable(namespace=namespace, name="plink", \
						version=version, os=operatingSystem, arch=architecture, installed=True)
		plink.addPFN(PFN("file://" + self.plinkPath, site_handler))
		executableClusterSizeMultiplierList.append((plink, 1))
		
		#2012.8.10 different plinks so that you can differentiate between different types of plink jobs
		plinkMerge =  Executable(namespace=namespace, name="plinkMerge", \
						version=version, os=operatingSystem, arch=architecture, installed=True)
		plinkMerge.addPFN(PFN("file://" + self.plinkPath, site_handler))
		executableClusterSizeMultiplierList.append((plinkMerge, 0))
		
		plinkIBD =  Executable(namespace=namespace, name="plinkIBD", \
						version=version, os=operatingSystem, arch=architecture, installed=True)
		plinkIBD.addPFN(PFN("file://" + self.plinkPath, site_handler))
		executableClusterSizeMultiplierList.append((plinkIBD, 0))
		
		plinkConvert =  Executable(namespace=namespace, name="plinkConvert", \
						version=version, os=operatingSystem, arch=architecture, installed=True)
		plinkConvert.addPFN(PFN("file://" + self.plinkPath, site_handler))
		executableClusterSizeMultiplierList.append((plinkConvert, 1))
		
		plinkLDPrune =  Executable(namespace=namespace, name="plinkLDPrune", \
						version=version, os=operatingSystem, arch=architecture, installed=True)
		plinkLDPrune.addPFN(PFN("file://" + self.plinkPath, site_handler))
		executableClusterSizeMultiplierList.append((plinkLDPrune, 1))
		
		plinkExtract =  Executable(namespace=namespace, name="plinkExtract", \
						version=version, os=operatingSystem, arch=architecture, installed=True)
		plinkExtract.addPFN(PFN("file://" + self.plinkPath, site_handler))
		executableClusterSizeMultiplierList.append((plinkExtract, 1))
		
		plinkNoClustering =  Executable(namespace=namespace, name="plinkNoClustering", \
						version=version, os=operatingSystem, arch=architecture, installed=True)
		plinkNoClustering.addPFN(PFN("file://" + self.plinkPath, site_handler))
		executableClusterSizeMultiplierList.append((plinkNoClustering, 0))
		
		ConvertBjarniSNPFormat2Yu = Executable(namespace=namespace, name="ConvertBjarniSNPFormat2Yu", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertBjarniSNPFormat2Yu.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "yhio/ConvertBjarniSNPFormat2Yu.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertBjarniSNPFormat2Yu, 1))
		
		ConvertVCF2BjarniFormat = Executable(namespace=namespace, name="ConvertVCF2BjarniFormat", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertVCF2BjarniFormat.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "pegasus/mapper/converter/ConvertVCF2BjarniFormat.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertVCF2BjarniFormat, 1))
		
		
		ConvertYuSNPFormat2Bjarni = Executable(namespace=namespace, name="ConvertYuSNPFormat2Bjarni", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertYuSNPFormat2Bjarni.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "yhio/ConvertYuSNPFormat2Bjarni.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2Bjarni, 1))
		
		ConvertYuSNPFormat2EigenStrat = Executable(namespace=namespace, name="ConvertYuSNPFormat2EigenStrat", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertYuSNPFormat2EigenStrat.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "yhio/ConvertYuSNPFormat2EigenStrat.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2EigenStrat, 1))
		
		ConvertYuSNPFormat2TPED_TFAM = Executable(namespace=namespace, name="ConvertYuSNPFormat2TPED_TFAM", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertYuSNPFormat2TPED_TFAM.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "yhio/ConvertYuSNPFormat2TPED_TFAM.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2TPED_TFAM, 1))
		
		
		DrawMatrix = Executable(namespace=namespace, name="DrawMatrix", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		DrawMatrix.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "plot/DrawMatrix.py"), site_handler))
		executableClusterSizeMultiplierList.append((DrawMatrix, 1))
		
		#2012.10.7
		Draw2DHistogramOfMatrix = Executable(namespace=namespace, name="Draw2DHistogramOfMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		Draw2DHistogramOfMatrix.addPFN(PFN("file://" + \
						os.path.join(self.pymodulePath, "plot/Draw2DHistogramOfMatrix.py"), site_handler))
		executableClusterSizeMultiplierList.append((Draw2DHistogramOfMatrix, 0))
		
		calcula = Executable(namespace=namespace, name="calcula", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		calcula.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py"), \
						site_handler))
		executableClusterSizeMultiplierList.append((calcula, 0.02))
		
		CalculatePairwiseDistanceOutOfSNPXStrainMatrix = Executable(namespace=namespace, name="CalculatePairwiseDistanceOutOfSNPXStrainMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		CalculatePairwiseDistanceOutOfSNPXStrainMatrix.addPFN(PFN("file://" + \
						os.path.join(self.pymodulePath, "pegasus/mapper/CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py"), \
						site_handler))
		executableClusterSizeMultiplierList.append((CalculatePairwiseDistanceOutOfSNPXStrainMatrix, 0.5))
		
		CalculateMedianMeanOfInputColumn = Executable(namespace=namespace, name="CalculateMedianMeanOfInputColumn", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		CalculateMedianMeanOfInputColumn.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/CalculateMedianMeanOfInputColumn"), site_handler))
		executableClusterSizeMultiplierList.append((CalculateMedianMeanOfInputColumn, 1))
		
		SampleRows = Executable(namespace=namespace, name="SampleRows", version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		SampleRows.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "statistics/SampleRows.py"), site_handler))
		executableClusterSizeMultiplierList.append((SampleRows, 0.5))
		
		#this is a reducer
		EstimateOutliersIn2DData = Executable(namespace=namespace, name="EstimateOutliersIn2DData", version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		EstimateOutliersIn2DData.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "statistics/EstimateOutliersIn2DData.py"), site_handler))
		executableClusterSizeMultiplierList.append((EstimateOutliersIn2DData, 0))
		
		#2013.2.3 use samtools to extract consensus from bam files
		ExtractConsensusSequenceFromAlignment = Executable(namespace=namespace, name="ExtractConsensusSequenceFromAlignment", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		ExtractConsensusSequenceFromAlignment.addPFN(PFN("file://" + \
			os.path.join(self.pymodulePath, "pegasus/mapper/alignment/ExtractConsensusSequenceFromAlignment.sh"), site_handler))
		executableClusterSizeMultiplierList.append((ExtractConsensusSequenceFromAlignment, 1))
		
		#2013.2.4, wrapper around psmc's splitfa, a program that splits fasta files
		splitfa = Executable(namespace=namespace, name="splitfa", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		splitfa.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/splitter/splitfa.sh"), site_handler))
		executableClusterSizeMultiplierList.append((splitfa, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
		#2013.2.7 convert, an image swissknife program, part of imagemagick 
		self.addOneExecutableFromPathAndAssignProperClusterSize(path="file:///usr/bin/convert", \
												name='convertImage', clusterSizeMultipler=1)
		#2013.2.10
		self.addOneExecutableFromPathAndAssignProperClusterSize(path="file://" + os.path.join(self.pymodulePath, "pegasus/mapper/runShellCommand.sh"), \
												name='runShellCommand', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'pegasus/mapper/converter/ConvertMSOutput2FASTQ.py'), \
										name='ConvertMSOutput2FASTQ', clusterSizeMultipler=1)
		
		#2013.2.11 moved from vervet/src/reduce to pymodule/pegasus/reducer
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'pegasus/reducer/MergeGenotypeMatrix.py'), \
										name='MergeGenotypeMatrix', clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'pegasus/reducer/MergeSameHeaderTablesIntoOne.py'), \
										name='mergeSameHeaderTablesIntoOne', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'pegasus/reducer/MergeSameHeaderTablesIntoOne.py'), \
										name='MergeSameHeaderTablesIntoOne', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'pegasus/reducer/ReduceMatrixByAverageColumnsWithSameKey.py'), \
										name='ReduceMatrixByAverageColumnsWithSameKey', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'pegasus/reducer/ReduceMatrixByChosenColumn.py'), \
										name='ReduceMatrixByChosenColumn', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'pegasus/reducer/ReduceMatrixByMergeColumnsWithSameKey.py'), \
										name='ReduceMatrixByMergeColumnsWithSameKey', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'pegasus/reducer/ReduceMatrixBySumSameKeyColsAndThenDivide.py'), \
										name='ReduceMatrixBySumSameKeyColsAndThenDivide', clusterSizeMultipler=0)
		
		
		
	def addExecutableAndAssignProperClusterSize(self, executableClusterSizeMultiplierList=[], defaultClustersSize=None):
		"""
		2012.8.31
			make sure the profile of clusters.size is not added already.
		2012.8.9
			
		"""
		if defaultClustersSize is None:
			defaultClustersSize = self.clusters_size
		for executableClusterSizeMultiplierTuple in executableClusterSizeMultiplierList:
			executable = executableClusterSizeMultiplierTuple[0]
			if len(executableClusterSizeMultiplierTuple)==1:
				clusterSizeMultipler = 1
			else:
				clusterSizeMultipler = executableClusterSizeMultiplierTuple[1]
			self.addOneExecutableAndAssignProperClusterSize(executable=executable, \
										clusterSizeMultipler=clusterSizeMultipler, defaultClustersSize=defaultClustersSize)
		
		
	def addOneExecutableAndAssignProperClusterSize(self, executable=None, clusterSizeMultipler=1, defaultClustersSize=None):
		"""
		2013.2.7, split out of addExecutableAndAssignProperClusterSize()
			clusterSizeMultipler could be any real value >=0. 0 means no clustering, 1=default clustering size.
			
			i.e.
			self.addOneExecutableAndAssignProperClusterSize(executable=CompareTwoGWAssociationLocusByPhenotypeVector, clusterSizeMultipler=0)
		"""
		executable = self.setOrChangeExecutableClusterSize(executable=executable, \
											clusterSizeMultipler=clusterSizeMultipler, defaultClustersSize=defaultClustersSize)
		if not self.hasExecutable(executable):
			self.addExecutable(executable)	#removeExecutable() is its counterpart
			setattr(self, executable.name, executable)
	
	def setOrChangeExecutableClusterSize(self, executable=None, clusterSizeMultipler=1, defaultClustersSize=None):
		"""
		2013.2.10
			split out of addOneExecutableAndAssignProperClusterSize()
			it'll remove the clustering profile if the new clusterSize is <1 
		"""
		if defaultClustersSize is None:
			defaultClustersSize = self.clusters_size
		clusterSize = int(defaultClustersSize*clusterSizeMultipler)
		clusteringProf = Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusterSize)
		if executable.hasProfile(clusteringProf):	#2012.8.26 check this first
			executable.removeProfile(clusteringProf)
		if clusterSize>1:
			executable.addProfile(clusteringProf)
		return executable
		
	
	def addOneExecutableFromPathAndAssignProperClusterSize(self, path=None, name=None, clusterSizeMultipler=1):
		"""
		2013.2.7
			combination of constructOneExecutableObject() & addOneExecutableAndAssignProperClusterSize()
		"""
		if clusterSizeMultipler is None:
			clusterSizeMultipler = 1
		executable = self.constructOneExecutableObject(path=path, name=name)
		self.addOneExecutableAndAssignProperClusterSize(executable=executable, clusterSizeMultipler=clusterSizeMultipler)
	
	
	def getFilesWithProperSuffixFromFolder(self, inputFolder=None, suffix='.h5'):
		"""
		2012.3.21
			moved from variation/src/FindGenomeWideLDPatternBetweenSNPsAndPeakWorkflow.py
		"""
		sys.stderr.write("Getting files with %s as suffix from %s ..."%(suffix, inputFolder))
		inputFnameLs = []
		counter = 0
		for filename in os.listdir(inputFolder):
			prefix, file_suffix = os.path.splitext(filename)
			counter += 1
			if file_suffix==suffix:
				inputFnameLs.append(os.path.join(inputFolder, filename))
		sys.stderr.write("%s files out of %s total.\n"%(len(inputFnameLs), counter))
		return inputFnameLs
	
	def getFilesWithSuffixFromFolderRecursive(self, inputFolder=None, suffixSet=set(['.h5']), fakeSuffix='.gz', inputFnameLs=[]):
		"""
		2012.4.30
			similar to getFilesWithProperSuffixFromFolder() but recursively go through all sub-folders
				and it uses utils.getRealPrefixSuffixOfFilenameWithVariableSuffix() to get the suffix.
		"""
		sys.stderr.write("Getting files with %s as suffix (%s as fake suffix) from %s ...\n"%(repr(suffixSet), fakeSuffix, inputFolder))
		counter = 0
		from pymodule import utils
		for filename in os.listdir(inputFolder):
			inputFname = os.path.join(inputFolder, filename)
			counter += 1
			if os.path.isfile(inputFname):
				prefix, file_suffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(filename, fakeSuffix=fakeSuffix)
				if file_suffix in suffixSet:
					inputFnameLs.append(inputFname)
			elif os.path.isdir(inputFname):
				self.getFilesWithSuffixFromFolderRecursive(inputFname, suffixSet=suffixSet, fakeSuffix=fakeSuffix, inputFnameLs=inputFnameLs)
		sys.stderr.write("%s files out of %s total.\n"%(len(inputFnameLs), counter))
		#return inputFnameLs
	
	def registerAllInputFiles(self, workflow=None, inputFnameLs=[], input_site_handler=None, pegasusFolderName='', **keywords):
		"""
		2012.3.9
			copied from variation.src.LDBetweenTwoSNPDataWorkflow.py
		2012.3.3
		"""
		sys.stderr.write("Registering %s input file ..."%(len(inputFnameLs)))
		returnData = PassingData(jobDataLs = [])
		counter = 0
		for inputFname in inputFnameLs:
			counter += 1
			inputF = File(os.path.join(pegasusFolderName, os.path.basename(inputFname)))
			inputF.addPFN(PFN("file://" + inputFname, input_site_handler))
			inputF.abspath = inputFname
			self.addFile(inputF)
			returnData.jobDataLs.append(PassingData(output=inputF, jobLs=[]))
		sys.stderr.write(" %s files registered.\n"%(len(returnData.jobDataLs)))
		return returnData
	
	def registerFilesAsInputToJob(self, job, inputFileList):
		"""
		2011-11-25
		"""
		for inputFile in inputFileList:
			job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
	
	def registerOneInputFile(self, workflow=None, inputFname=None, input_site_handler=None, folderName="", useAbsolutePathAsPegasusFileName=False,\
							pegasusFileName=None):
		"""
		2013.2.14 added argument useAbsolutePathAsPegasusFileName
			This would render the file to be referred as the absolute path on the running computer.
			And pegasus will not seek to symlink or copy/transfer the file.
			set it to True only when you dont want to add the file to the job as INPUT dependency (as it's accessed through abs path).
		2013.1.10 make sure the file is not registed with the workflow already
		2012.3.22
			add abspath attribute to file.
		2012.3.1
			add argument folderName, which will put the file in specific pegasus workflow folder
		2011.12.21
		"""
		if workflow is None:
			workflow = self
		if input_site_handler is None:
			input_site_handler = self.input_site_handler
		if not pegasusFileName:
			if useAbsolutePathAsPegasusFileName:
				pegasusFileName = os.path.abspath(inputFname)	#this will stop symlinking/transferring , and also no need to indicate them as file dependency for jobs.
			else:
				pegasusFileName = os.path.join(folderName, os.path.basename(inputFname))
		pegasusFile = File(pegasusFileName)
		pegasusFile.abspath = os.path.abspath(inputFname)
		pegasusFile.absPath = pegasusFile.abspath
		pegasusFile.addPFN(PFN("file://" + pegasusFile.abspath, input_site_handler))
		if not workflow.hasFile(pegasusFile):	#2013.1.10
			workflow.addFile(pegasusFile)
		return pegasusFile
	
	def registerOneJar(self, name=None, path=None, site_handler=None, workflow=None, folderName="", useAbsolutePathAsPegasusFileName=False):
		"""
		2013.2.14
			useAbsolutePathAsPegasusFileName=True if you do not plan to add a jar file as INPUT dependency for jobs
		"""
		if workflow is None:
			workflow = self
		if site_handler is None:
			site_handler = self.site_handler	#usually they are same
		if not folderName:
			folderName = "jar"
		pegasusFile = self.registerOneInputFile(workflow=workflow, inputFname=path, input_site_handler=site_handler, \
											folderName=folderName, useAbsolutePathAsPegasusFileName=useAbsolutePathAsPegasusFileName)
		setattr(workflow, name, pegasusFile)
	
	def registerBlastNucleotideDatabaseFile(self, ntDatabaseFname=None, input_site_handler=None, folderName=""):
		"""
		2012.10.8
			moved from BlastWorkflow.py
		2012.5.23
		"""
		if input_site_handler is None:
			input_site_handler = self.input_site_handler
		return yh_pegasus.registerRefFastaFile(self, ntDatabaseFname, registerAffiliateFiles=True, \
									input_site_handler=input_site_handler,\
									checkAffiliateFileExistence=True, addPicardDictFile=False, \
									affiliateFilenameSuffixLs=['nin', 'nhr', 'nsq'],\
									folderName=folderName)
	
	
	def addStatMergeJob(self, workflow=None, statMergeProgram=None, outputF=None, \
					parentJobLs=None, extraOutputLs=None,\
					extraDependentInputLs=None, transferOutput=True, extraArguments=None, \
					namespace=None, version=None, job_max_memory=1000, **keywords):
		"""
		2012.8.10
			use addGenericJob()
		2012.4.3
			make argument namespace, version optional
		2011-11-28
			moved from CalculateVCFStatPipeline.py
		2011-11-17
			add argument extraArguments
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		extraArgumentList = []
		key2ObjectForJob = {}
		
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericJob(executable=statMergeProgram, inputFile=None, outputFile=outputF, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
		return job
	
	def addInputToStatMergeJob(self, workflow=None, statMergeJob=None, inputF=None, \
							parentJobLs=None, \
							namespace=None, version=None, extraDependentInputLs=None):
		"""
		i.e. :	self.addInputToStatMergeJob(statMergeJob=associationLocusJob, parentJobLs=[associationPeakJob])
		2013.2.6 make sure parentJob is of instance Job
		2012.10.8
			inputF is optional, if not given, parentJobLs must be given, and parentJobLs[0].output is inputF. 
		2012.9.12 parentJobLs and extraDependentInputLs could be None
		2011-11-28
			moved from CalculateVCFStatPipeline.py
		"""
		if workflow is None:
			workflow = self
		if inputF is None and parentJobLs is not None:	#2012.10.8
			parentJob = parentJobLs[0]
			if hasattr(parentJob, 'output'):
				inputF = parentJob.output
		if inputF:
			isAdded = self.addJobUse(statMergeJob, file=inputF, transfer=True, register=True, link=Link.INPUT)
			if isAdded:
				statMergeJob.addArguments(inputF)
		
		if extraDependentInputLs:
			for inputFile in extraDependentInputLs:
				if inputFile:
					isAdded = self.addJobUse(statMergeJob, file=inputFile, transfer=True, register=True, link=Link.INPUT)
		
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					self.addJobDependency(parentJob=parentJob, childJob=statMergeJob)
	
	
	def addConvertImageJob(self, inputFile=None, inputArgumentOption=None, \
					outputFile=None, outputArgumentOption=None, density=None, \
					resizeDimension=None, \
					parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, job_max_memory=2000,\
					key2ObjectForJob=None, **keywords):
		"""
		2013.2.7 use imagemagick's convert to convert images. examples:
			plotOutputFile = File('%s.eps'%(plotOutputFnamePrefix))
			plotPNGOutputFile = File('%s.png'%(plotOutputFnamePrefix))
			#change dpi to 300
			self.addConvertImageJob(inputFile=plotOutputFile, \
					outputFile=plotPNGOutputFile, density=300, \
					resizeDimension=None, \
					parentJobLs=[psmc_plotJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
					extraArguments=None, frontArgumentList=None, job_max_memory=500)
			
			#resize by demanding the width = 1800, height will scale accordingly
			self.addConvertImageJob(inputFile=plotOutputFile, \
					outputFile=plotPNGOutputFile, density=None, \
					resizeDimension=1800, \
					parentJobLs=[psmc_plotJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
					extraArguments=None, frontArgumentList=None, job_max_memory=500)
			#resize by demanding the dimension=1800X900
			self.addConvertImageJob(inputFile=plotOutputFile, \
					outputFile=plotPNGOutputFile, density=None, \
					resizeDimension='1800X900', \
					parentJobLs=[psmc_plotJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
					extraArguments=None, frontArgumentList=None, job_max_memory=500)
		"""
		if frontArgumentList is None:
			frontArgumentList = []
		if extraOutputLs is None:
			extraOutputLs = []
		if density is not None:
			frontArgumentList.append("-density %s"%(density))
		if resizeDimension is not None:
			frontArgumentList.append("-resize %s"%(resizeDimension))
		#do not pass the inputFileList to addGenericJob() because db arguments need to be added before them. 
		job = self.addGenericJob(workflow=self, executable=self.convertImage, inputFile=inputFile, \
						inputArgumentOption=inputArgumentOption, outputFile=outputFile, \
						outputArgumentOption=outputArgumentOption, inputFileList=None, parentJobLs=parentJobLs, \
						extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
						transferOutput=transferOutput, \
						frontArgumentList=frontArgumentList, extraArguments=extraArguments, extraArgumentList=extraArgumentList,\
						job_max_memory=job_max_memory, key2ObjectForJob=key2ObjectForJob,\
						**keywords)
		return job
	
	def addGenericDBJob(self, workflow=None, executable=None, inputFile=None, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", inputFileList=None, \
					parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					extraArguments=None, extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
					key2ObjectForJob=None, objectWithDBArguments=None, **keywords):
		"""
		2012.10.6
			similar to addGenericJob but these are generic jobs that need db-interacting arguments.
			
			Argument inputFileList is a list of input files to be added to commandline as the last arguments.
				they would also be added as the job's dependent input.
				Difference from extraDependentInputLs: the latter is purely for dependency purpose, not added as input arguments.
					So if files have been put in inputFileList, then they shouldn't be in extraDependentInputLs.
		"""
		#do not pass the inputFileList to addGenericJob() because db arguments need to be added before them. 
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=inputFile, \
						inputArgumentOption=inputArgumentOption, outputFile=outputFile, \
						outputArgumentOption=outputArgumentOption, inputFileList=None, parentJobLs=parentJobLs, \
						extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
						transferOutput=transferOutput, extraArguments=extraArguments, extraArgumentList=extraArgumentList,\
						job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel, key2ObjectForJob=key2ObjectForJob,\
						**keywords)
		if objectWithDBArguments is None:
			objectWithDBArguments = self
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=objectWithDBArguments)
				
		#2012.10.6 add all input files to the last (after db arguments,) otherwise, it'll mask others (cuz these don't have options).
		if inputFileList:
			for inputFile in inputFileList:
				if inputFile:
					isAdded = self.addJobUse(job, file=inputFile, transfer=True, register=True, link=Link.INPUT)
					if isAdded:	#2012.11.21
						job.addArguments(inputFile)
						job.inputLs.append(inputFile)
		
		#2012.10.6 set the job.input
		if getattr(job, 'input', None) is None and job.inputLs:
			job.input = job.inputLs[0]
		return job
		
	def addGenericFile2DBJob(self, workflow=None, executable=None, inputFile=None, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", inputFileList=None, \
					data_dir=None, logFile=None, commit=False,\
					parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					extraArguments=None, extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
					key2ObjectForJob=None, objectWithDBArguments=None, **keywords):
		"""
		2012.12.21 a generic wrapper for jobs that "inserting" data (from file) into database
		"""
		if extraArgumentList is None:
			extraArgumentList = []
		if extraOutputLs is None:
			extraOutputLs = []
		
		if data_dir:
			extraArgumentList.append('--data_dir %s'%(data_dir))
		if commit:
			extraArgumentList.append('--commit')
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
			extraOutputLs.append(logFile)
		#do not pass the inputFileList to addGenericJob() because db arguments need to be added before them. 
		job = self.addGenericDBJob(workflow=workflow, executable=executable, inputFile=inputFile, \
						inputArgumentOption=inputArgumentOption, outputFile=outputFile, \
						outputArgumentOption=outputArgumentOption, inputFileList=inputFileList, parentJobLs=parentJobLs, \
						extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
						transferOutput=transferOutput, extraArguments=extraArguments, extraArgumentList=extraArgumentList,\
						job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel, key2ObjectForJob=key2ObjectForJob,\
						objectWithDBArguments=objectWithDBArguments, **keywords)
		return job
	
	def addJobUse(self, job=None, file=None, transfer=True, register=True, link=None):
		"""
		2012.10.18 check whether that file is a use of job already.
		"""
		use = Use(file.name, link=link, register=register, transfer=transfer, optional=None, \
								namespace=job.namespace,\
								version=job.version, executable=None)	#, size=None
		if job.hasUse(use):
			return False
		else:
			job.addUse(use)
			return True
	
	def addJobDependency(self, workflow=None, parentJob=None, childJob=None):
		"""
		2013.2.6 make sure parentJob is of instance Job, sometimes, it could be a fake job, like PassingData(output=...)
		2012.10.18 check whether that the dependency exists already
		"""
		if workflow is None:
			workflow = self
		addedOrNot = True
		if isinstance(parentJob, Job):
			dep = Dependency(parent=parentJob, child=childJob)
			if not workflow.hasDependency(dep):
				workflow.addDependency(dep)
				addedOrNot = True
			else:
				addedOrNot = False
		else:	#2013.2.6
			#sys.stderr.write("Warning: parent job %s is not a Job-instance.\n"%(repr(parentJob)))
			addedOrNot = False
		return addedOrNot
		
		
	def addGenericJob(self, workflow=None, executable=None, inputFile=None, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", inputFileList=None, \
					parentJob=None, parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=None, **keywords):
		"""
		2013.2.7 add argument frontArgumentList, a list of arguments to be put in front of anything else
		2012.10.18
			add job.parentJobLs to store its parent job(s).
		2012.10.15
			add argument parentJob, same as parentJobLs, but just one job
		2012.10.6
			add argument inputFileList, which is a list of input files to be added to commandline as the last arguments
				they would also be added as the job's dependent input.
				Difference from extraDependentInputLs: the latter is purely for dependency purpose, not added as input arguments.
					So if files have been put in inputFileList, then they shouldn't be in extraDependentInputLs.
		2012.8.17 if transferOutput is None, do not register output files as OUTPUT with transfer flag
		2012.8.2
			add argument inputArgumentOption, outputArgumentOption so that user 
		2012.7.31 add argument key2ObjectForJob, which is a dictionary with strings as key, to set key:object for each job
		#2012.7.28 if job.output is not set, set it to the 1st entry of job.outputLs
		2012.6.27 add job.outputLs to hold more output files.
		2012.6.1
			add argument extraOutputLs
		2012.5.24
			generic job addition function for other functions to use
		"""
		if workflow is None:
			workflow =self
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.outputLs = []	#2012.6.27 to hold more output files
		job.inputLs = []
		if frontArgumentList:	#2013.2.7
			job.addArguments(*frontArgumentList)
		if inputFile:
			if inputArgumentOption:
				job.addArguments(inputArgumentOption)
			job.addArguments(inputFile)
			isAdded = self.addJobUse(job, file=inputFile, transfer=True, register=True, link=Link.INPUT)
			#job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
			job.input = inputFile
			job.inputLs.append(inputFile)
		if outputFile:
			if outputArgumentOption:
				job.addArguments(outputArgumentOption)
			job.addArguments(outputFile)
			if transferOutput is not None:	#2012.8.17
				self.addJobUse(job, file=outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
				#job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			job.output = outputFile
			job.outputLs.append(outputFile)
		if extraArgumentList:
			job.addArguments(*extraArgumentList)
		
		if extraArguments:
			job.addArguments(extraArguments)
		
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel,\
									no_of_cpus=no_of_cpus, max_walltime=max_walltime)
		workflow.addJob(job)
		job.parentJobLs = []	#2012.10.18
		if parentJob:	#2012.10.15
			isAdded = self.addJobDependency(workflow=workflow, parentJob=parentJob, childJob=job)
			if isAdded:
				job.parentJobLs.append(parentJob)
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					isAdded = self.addJobDependency(workflow=workflow, parentJob=parentJob, childJob=job)
					if isAdded:
						job.parentJobLs.append(parentJob)
		if extraDependentInputLs:
			for input in extraDependentInputLs:
				if input:
					isAdded = self.addJobUse(job, file=input, transfer=True, register=True, link=Link.INPUT)
					#job.uses(input, transfer=True, register=True, link=Link.INPUT)
					if isAdded:
						job.inputLs.append(input)
		if extraOutputLs:
			for output in extraOutputLs:
				if output:
					job.outputLs.append(output)
					if transferOutput is not None:	#2012.8.17
						self.addJobUse(job, file=output, transfer=transferOutput, register=True, link=Link.OUTPUT)
						#job.uses(output, transfer=transferOutput, register=True, link=Link.OUTPUT)
		if key2ObjectForJob:
			for key, object in key2ObjectForJob.iteritems():
				setattr(job, key, object)	#key should be a string.
		
		#2012.10.6 add all input files to the last (after db arguments,) otherwise, it'll mask others (cuz these don't have options).
		if inputFileList:
			for inputFile in inputFileList:
				if inputFile:
					job.addArguments(inputFile)
					isAdded = self.addJobUse(job, file=inputFile, transfer=True, register=True, link=Link.INPUT)
					#job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
					if isAdded:
						job.inputLs.append(inputFile)
		#2012.10.9 make sure outputList
		job.outputList = job.outputLs
		#2012.7.28 if job.output is not set, set it to the 1st entry of job.outputLs
		if getattr(job, 'output', None) is None and job.outputLs:
			job.output = job.outputLs[0]
		if getattr(job, 'input', None) is None and job.inputLs:
			job.input = job.inputLs[0]
		self.no_of_jobs += 1
		return job
	
	def setJobOutputFileTransferFlag(self, job=None, transferOutput=False, outputLs=None):
		"""
		2012.8.17
			assume all output files in job.outputLs
		"""
		if outputLs is None and getattr(job, 'outputLs', None):
			outputLs = job.outputLs
		
		for output in outputLs:
			job.uses(output, transfer=transferOutput, register=True, link=Link.OUTPUT)
	
	def addCalculateDepthMeanMedianModeJob(self, workflow=None, executable=None, \
							inputFile=None, outputFile=None, alignmentID=None, fractionToSample=0.001, \
							noOfLinesInHeader=0, whichColumn=2, maxNumberOfSamplings=1E7, inputStatName=None,\
							parentJobLs=None, job_max_memory = 500, extraArguments=None, \
							transferOutput=False, **keywords):
		"""
		2013.1.8 moved from vervet.src.alignment.InspectAlignmentPipeline and use addGenericJob()
		2012.6.15 turn maxNumberOfSamplings into integer when passing it to the job
		2012.5.7
			a job to take input of samtoolsDepth
		"""
		extraArgumentList = []
		if alignmentID is not None:
			extraArgumentList.append("--alignmentID %s"%(alignmentID))
		if noOfLinesInHeader is not None:
			extraArgumentList.append("--noOfLinesInHeader %s"%(noOfLinesInHeader))
		if fractionToSample is not None:
			extraArgumentList.append("--fractionToSample %s"%(fractionToSample))
		if whichColumn is not None:
			extraArgumentList.append("--whichColumn %s"%(whichColumn))
		if maxNumberOfSamplings is not None:
			extraArgumentList.append('--maxNumberOfSamplings %d'%(maxNumberOfSamplings))
		if inputStatName is not None:
			extraArgumentList.append("--inputStatName %s"%(inputStatName))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericJob(workflow=workflow, executable=executable, inputFile=inputFile, outputFile=outputFile, \
				parentJobLs=parentJobLs, extraDependentInputLs=None, \
				extraOutputLs=None,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=None, \
				sshDBTunnel=None, job_max_memory=job_max_memory, **keywords)
		return job
	
	def addDBArgumentsToOneJob(self, job=None, objectWithDBArguments=None):
		"""
		2012.8.17
			use long arguments , rather than short ones
		2012.6.5
			tired of adding all these arguments to db-interacting jobs
		"""
		if objectWithDBArguments is None:
			objectWithDBArguments = self
		job.addArguments("--drivername", objectWithDBArguments.drivername, "--hostname", objectWithDBArguments.hostname, \
						"--dbname", objectWithDBArguments.dbname, \
						"--db_user", objectWithDBArguments.db_user, "--db_passwd %s"%objectWithDBArguments.db_passwd)
		if objectWithDBArguments.schema:
			job.addArguments("--schema", objectWithDBArguments.schema)
		if getattr(objectWithDBArguments, 'port', None):
			job.addArguments("--port=%s"%(objectWithDBArguments.port))
		return job
	
	def addGzipSubWorkflow(self, workflow=None, inputData=None, transferOutput=True,\
						outputDirPrefix="", topOutputDirJob=None, report=True, **keywords):
		"""
		2012.8.2 bugfix.
		2012.7.19
		"""
		if workflow is None:
			workflow = self
		if report:
			sys.stderr.write("Adding gzip jobs for %s input job data ... "%(len(inputData.jobDataLs)))
		
		returnData = PassingData(topOutputDirJob=None)
		returnData.jobDataLs = []
		if len(inputData.jobDataLs)>0:
			if topOutputDirJob is None:
				topOutputDir = "%sGzip"%(outputDirPrefix)
				topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
			
			returnData.topOutputDirJob = topOutputDirJob
			for jobData in inputData.jobDataLs:
				for inputF in jobData.fileList:
					inputFBaseName = os.path.basename(inputF.name)
					outputF = File(os.path.join(topOutputDirJob.output, '%s.gz'%(inputFBaseName)))
					key2ObjectForJob = {}
					extraArgumentList = []
					#make sure set inputArgumentOption&outputArgumentOption to None, \
					# otherwise addGenericJob will add "-i" and "-o" in front of it
					job= self.addGenericJob(workflow=workflow, executable=workflow.gzip, inputFile=inputF,
								inputArgumentOption=None, outputArgumentOption=None,  outputFile=outputF, \
								parentJobLs=[topOutputDirJob]+jobData.jobLs, extraDependentInputLs=None, \
								extraOutputLs=[],\
								transferOutput=transferOutput, \
								extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, \
								job_max_memory=200, **keywords)
					"""	
					# 2012.8.2 wrong, because -i and -o will be added in front.
					abstractMapperJob = self.addAbstractMapperLikeJob(workflow, executable=workflow.gzip, \
							inputF=None, outputF=outputF, \
							parentJobLs=[topOutputDirJob]+jobData.jobLs, transferOutput=transferOutput, job_max_memory=200,\
							extraArguments=None, extraDependentInputLs=[inputF], )
					"""
					returnData.jobDataLs.append(PassingData(jobLs=[job], vcfFile=None, file=outputF,\
											fileList=[outputF]))
		if report:
			sys.stderr.write("no_of_jobs = %s.\n"%(self.no_of_jobs))
		return returnData
	
	def addAbstractMapperLikeJob(self, workflow=None, executable=None, \
					inputVCF=None, inputF=None, outputF=None, extraOutputLs=None,\
					parentJobLs=None, transferOutput=True, job_max_memory=2000,\
					extraArguments=None, extraArgumentList=None, extraDependentInputLs=None, \
					sshDBTunnel=None, **keywords):
		"""
		2012.10.8 call addGenericJob() instead
		2012.7.19
			moved from AbstractNGSWorkflow to here.
			add argument inputF. inputVCF is not generic enough.
		2012.5.11
		"""
		if inputF is None:	#2012.7.19
			inputF = inputVCF
		job= self.addGenericJob(workflow=workflow, executable=executable, inputFile=inputF, outputFile=outputF, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArguments=extraArguments,\
				extraArgumentList=extraArgumentList, key2ObjectForJob=None, \
				sshDBTunnel=sshDBTunnel, job_max_memory=job_max_memory, **keywords)
		return job
	
	def addSelectLineBlockFromFileJob(self, executable=None, inputFile=None, outputFile=None,\
					startLineNumber=None, stopLineNumber=None, parentJobLs=None, extraDependentInputLs=None, \
					transferOutput=False, \
					extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.7.30
		"""
		extraArgumentList = ['-s %s'%(startLineNumber), '-t %s'%(stopLineNumber)]
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		return job
	
	def getJVMMemRequirment(self, job_max_memory=4000, minMemory=2000):
		"""
		2012.8.2
			job_max_memory could be set by user to lower than minMemory.
			but minMemory makes sure it's never too low.
		"""
		MaxPermSize_user = job_max_memory*2/5
		mxMemory_user = job_max_memory*3/5
		MaxPermSize= min(35000, max(minMemory, MaxPermSize_user))
		PermSize=MaxPermSize*3/4
		mxMemory = max(minMemory, job_max_memory)
		msMemory = mxMemory*3/4
		memRequirementInStr = "-Xms%sm -Xmx%sm -XX:PermSize=%sm -XX:MaxPermSize=%sm"%\
				(msMemory, mxMemory, PermSize, MaxPermSize)
		memRequirement = MaxPermSize + mxMemory
		
		return PassingData(memRequirementInStr=memRequirementInStr, memRequirement=memRequirement)
	
	def addPlotLDJob(self, workflow=None, executable=None, inputFile=None, inputFileList=None, outputFile=None, \
					outputFnamePrefix=None,
					whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, title=None, \
					logY=None, valueForNonPositiveYValue=-1, \
					missingDataNotation='-nan',\
					xColumnPlotLabel=None, chrLengthColumnHeader=None, chrColumnHeader=None, \
					minChrLength=1000000, xColumnHeader=None, pos2ColumnHeader=None, minNoOfTotal=100, maxNoOfTotal=None,\
					figureDPI=300, formatString='.', ylim_type=2, samplingRate=0.0001,  need_svg=False, logCount=False, \
					minDist=None, maxDist=None, movingAverageType=2,\
					parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
		"""
		2012.10.25
			expose argument missingDataNotation, minDist, maxDist
		2012.8.18
			use addAbstractPlotJob()
		2012.8.2 moved from vervet/src/CalculateVCFStatPipeline.py
		2012.8.1
		('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
			('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
			('title', 1, ): [None, 't', 1, 'title for the figure.'],\
			('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
			('formatString', 1, ): ['.', 'a', 1, 'formatString passed to matplotlib plot'],\
			('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
			('samplingRate', 1, float): [0.001, 's', 1, 'how often you include the data'],\
			('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
			('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnHeader', 0, ): ["", 'W', 1, 'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
			('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, 'plot label for data of the whichColumn', ],\
			('chrLengthColumnHeader', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
			('chrColumnHeader', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
			('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
			('pos1ColumnLabel', 1, ): ['POS1', 'l', 1, 'label of the 1st position column', ],\
			('pos2ColumnLabel', 1, ): ['POS2', 'p', 1, 'label of the 2nd position column', ],\
			('posColumnPlotLabel', 1, ): ['distance', 'x', 1, 'x-axis label in  plot', ],\
			
		"""
		if extraArguments is None:
			extraArguments = ""
		if extraArgumentList is None:
			extraArgumentList = []
		if logCount:
			extraArguments += " --logCount "
		if minChrLength is not None:
			extraArguments += " --minChrLength %s "%(minChrLength)
		if chrLengthColumnHeader:
			extraArgumentList.append("--chrLengthColumnHeader %s"%(chrLengthColumnHeader))
		if chrColumnHeader:
			extraArgumentList.append("--chrColumnHeader %s"%(chrColumnHeader))
		if pos2ColumnHeader:
			extraArgumentList.append(' --pos2ColumnHeader %s '%(pos2ColumnHeader))
		if minDist:
			extraArgumentList.append('--minDist %s'%(minDist))
		if maxDist:
			extraArgumentList.append('--maxDist %s'%(maxDist))
		if maxNoOfTotal:
			extraArgumentList.append("--maxNoOfTotal %s"%(maxNoOfTotal))
		if movingAverageType:
			extraArgumentList.append("--movingAverageType %s"%(movingAverageType))
		
		return self.addAbstractPlotJob(workflow=workflow, executable=executable, inputFileList=inputFileList, \
							inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
							whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=whichColumnPlotLabel, \
							logY=logY, valueForNonPositiveYValue=valueForNonPositiveYValue, \
							missingDataNotation=missingDataNotation,\
							xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel, title=title, \
							minNoOfTotal=minNoOfTotal, \
							figureDPI=figureDPI, formatString=formatString, ylim_type=ylim_type, samplingRate=samplingRate, need_svg=need_svg, \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraArgumentList=extraArgumentList, extraArguments=extraArguments, transferOutput=transferOutput, \
							job_max_memory=job_max_memory, \
							**keywords)
		
		
	
	def addAbstractPlotJob(self, workflow=None, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
					logX=None, logY=None, valueForNonPositiveYValue=-1, \
					xScaleLog=0, yScaleLog=0, \
					missingDataNotation='NA',\
					xColumnHeader=None, xColumnPlotLabel=None, title=None, \
					minNoOfTotal=100, maxNoOfTotal=None,\
					figureDPI=300, formatString='.', ylim_type=2, samplingRate=0.001, need_svg=False, \
					inputFileFormat=None, outputFileFormat=None,\
					parentJob=None, parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, \
					sshDBTunnel=False, \
					objectWithDBArguments=None, **keywords):
		"""
		2012.12.3 added argument title, logX, logY
		2012.10.16 added argument sshDBTunnel, objectWithDBArguments
		2012.8.31 add argument missingDataNotation
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		2012.8.2 (check AbstractMatrixFileWalker.py or AbstractPlot.py for updated arguments)
			('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
			('minNoOfTotal', 1, int): [100, 'M', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
			('title', 0, ): [None, 't', 1, 'title for the figure.'],\
			('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
			('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
			('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
			('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnHeader', 0, ): ["", 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
			('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
			('logY', 0, int): [0, '', 1, 'value 0: nothing; 1: log(), 2: -log(). replacing self.logWhichColumn.'],\
			('valueForNonPositiveYValue', 1, float): [50, '', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
					what yValue should be.'],\
			('xColumnHeader', 1, ): ['', 'l', 1, 'header of the x-axis data column, ' ],\
			('xColumnPlotLabel', 0, ): ['', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
			('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
			
			inputFileFormat   1: csv-like plain text file; 2: YHPyTables.YHFile; 3: HDF5MatrixFile; . "1"(default)
			
		"""
		extraOutputLs = []
		key2ObjectForJob = {}
		if executable is None:
			executable = self.AbstractPlot
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if inputFileList:
			extraDependentInputLs.extend(inputFileList)
		if extraArgumentList is None:
			extraArgumentList = []
		
		if outputFnamePrefix:
			extraArgumentList.append('--outputFnamePrefix %s'%(outputFnamePrefix))
			if outputFile is None:
				extraOutputLs.append(File('%s.png'%(outputFnamePrefix)))
				if need_svg:
					extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
		if minNoOfTotal is not None:
			extraArgumentList.append('--minNoOfTotal %s'%(minNoOfTotal))
		if maxNoOfTotal:
			extraArgumentList.append("--maxNoOfTotal %s"%(maxNoOfTotal))
		if figureDPI:
			extraArgumentList.append('--figureDPI %s'%(figureDPI))
		if formatString:
			extraArgumentList.append('--formatString %s'%(formatString))
		if ylim_type:
			extraArgumentList.append('--ylim_type %s'%(ylim_type))
		if samplingRate is not None:
			extraArgumentList.append('--samplingRate %s'%(samplingRate))
		
		if xColumnHeader:
			extraArgumentList.append('--xColumnHeader %s'%(xColumnHeader))
		if xColumnPlotLabel:
			extraArgumentList.append("--xColumnPlotLabel %s"%(xColumnPlotLabel))
		if whichColumnHeader:
			extraArgumentList.append("--whichColumnHeader %s"%(whichColumnHeader))
		if whichColumn:
			extraArgumentList.append("--whichColumn %s"%(whichColumn))
		if whichColumnPlotLabel:
			extraArgumentList.append("--whichColumnPlotLabel %s"%(whichColumnPlotLabel))
		if title:
			extraArgumentList.append("--title %s"%(title))
		if logX:
			extraArgumentList.append("--logX %s"%(logX))
		if logY:
			extraArgumentList.append('--logY %s'%(logY))
		if xScaleLog:
			extraArgumentList.append("--xScaleLog %s"%(xScaleLog))
		if yScaleLog:
			extraArgumentList.append("--yScaleLog %s"%(yScaleLog))
		
		if valueForNonPositiveYValue:
			extraArgumentList.append("--valueForNonPositiveYValue %s"%(valueForNonPositiveYValue))
		if inputFileFormat:
			extraArgumentList.append("--inputFileFormat %s"%(inputFileFormat))
		if outputFileFormat:
			extraArgumentList.append("--outputFileFormat %s"%(outputFileFormat))
		if need_svg:
			extraArgumentList.append('--need_svg')
			if not outputFnamePrefix:
				outputFnamePrefix = os.path.splitext(outputFile.name)[0]	#2012.8.20 bugfix.
			extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		if objectWithDBArguments:
			#if this is given, then it's a db interacting job, and inputFileList is not added into extraArgumentList
			job = self.addGenericDBJob(workflow=workflow, executable=executable, inputFile=inputFile, outputFile=outputFile, \
				inputFileList = inputFileList,\
				parentJob=parentJob, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, objectWithDBArguments=objectWithDBArguments, **keywords)
		else:
			#add all input files to the last
			#only for non-db interacting jobs
			if inputFileList:
				extraArgumentList.extend(inputFileList)
			job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=inputFile, outputFile=outputFile, \
				parentJob=parentJob, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		return job
	
	def addAbstractMatrixFileWalkerJob(self, workflow=None, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, \
					logY=None, valueForNonPositiveYValue=-1, \
					minNoOfTotal=10,\
					samplingRate=1, \
					inputFileFormat=None, outputFileFormat=None,\
					parentJob=None, parentJobLs=None, \
					extraDependentInputLs=None, extraArgumentList=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000, sshDBTunnel=False, \
					objectWithDBArguments=None, **keywords):
		"""
		2012.11.25 more arguments, logY, inputFileFormat, outputFileFormat
		2012.10.16 added argument sshDBTunnel, objectWithDBArguments
		2012.10.15 added extraArgumentList, parentJob
		2012.8.15
			('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
			('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
			('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnHeader', 0, ): ["", 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
			('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
			('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
				only effective when logWhichColumn is toggled. '],\
			('valueForNonPositiveYValue', 1, float): [50, '', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
					what yValue should be.'],\
			
		"""
		return self.addAbstractPlotJob(workflow=workflow, executable=executable, inputFileList=inputFileList, \
							inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
							whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=None, \
							logY=logY, \
							valueForNonPositiveYValue=valueForNonPositiveYValue, \
							missingDataNotation=None,\
							xColumnHeader=None, xColumnPlotLabel=None, \
							minNoOfTotal=minNoOfTotal, \
							figureDPI=None, formatString=None, ylim_type=None, samplingRate=samplingRate, need_svg=False, \
							parentJob=parentJob, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraArgumentList=extraArgumentList,\
							extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
							sshDBTunnel=sshDBTunnel, objectWithDBArguments=objectWithDBArguments, **keywords)
	
	def addDrawHistogramJob(self, workflow=None, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
					logY=None, valueForNonPositiveYValue=-1, missingDataNotation='NA', title=None, \
					minNoOfTotal=10,\
					figureDPI=100, formatString='.', ylim_type=2, samplingRate=0.001, need_svg=False, \
					logCount=False, inputFileFormat=None, \
					parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
		"""
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		2012.8.2
			('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
			('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
			('title', 0, ): [None, 't', 1, 'title for the figure.'],\
			('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
			('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
			('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
			('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnHeader', 0, ): ["", 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
			('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
			('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
				only effective when logWhichColumn is toggled. '],\
			('valueForNonPositiveYValue', 1, float): [50, '', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
					what yValue should be.'],\
			('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
			
		"""
		if extraArguments is None:
			extraArguments = ""
		if logCount:
			extraArguments += " --logCount "
		return self.addAbstractPlotJob(workflow=workflow, executable=executable, inputFileList=inputFileList, \
							inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
							whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=whichColumnPlotLabel, \
							logY=logY, valueForNonPositiveYValue=valueForNonPositiveYValue, \
							missingDataNotation=missingDataNotation,\
							xColumnHeader=None, xColumnPlotLabel=None, title=title, \
							minNoOfTotal=minNoOfTotal, \
							figureDPI=figureDPI, formatString=formatString, ylim_type=ylim_type, \
							samplingRate=samplingRate, need_svg=need_svg, \
							inputFileFormat=inputFileFormat,\
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
							**keywords)
	
	def addDraw2DHistogramOfMatrixJob(self, workflow=None, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
				outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
				logX=False, logY=False, logZ=False, valueForNonPositiveYValue=-1, \
				missingDataNotation='NA',\
				xColumnHeader=None, xColumnPlotLabel=None, \
				minNoOfTotal=100,\
				figureDPI=300, formatString='.', samplingRate=0.001, need_svg=False, \
				inputFileFormat=None, outputFileFormat=None,\
				zColumnHeader=None, \
				parentJobLs=None, \
				extraDependentInputLs=None, \
				extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
		"""
		2013.2.8 added argument inputFileFormat
		2013.2.7 executable could be None, default is self.Draw2DHistogramOfMatrix
		2012.11.28 change logX, logY, logZ
		2012.10.7
			
		"""
		if extraArgumentList is None:
			extraArgumentList = []
		if zColumnHeader:
			extraArgumentList.append("--zColumnHeader %s"%(zColumnHeader))
		if logZ:
			extraArgumentList.append("--logZ %s"%(logZ))
		if executable is None:
			executable = self.Draw2DHistogramOfMatrix
		return self.addAbstractPlotJob(workflow=workflow, executable=executable, inputFileList=inputFileList, \
							inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
							whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=whichColumnPlotLabel, \
							logX=logX, logY=logY, valueForNonPositiveYValue=valueForNonPositiveYValue, \
							missingDataNotation=missingDataNotation,\
							xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel, \
							minNoOfTotal=minNoOfTotal, \
							figureDPI=figureDPI, formatString=formatString, ylim_type=None, \
							samplingRate=samplingRate, need_svg=need_svg, \
							inputFileFormat=inputFileFormat, outputFileFormat=outputFileFormat,\
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
							**keywords)
	
	def addMkDirJob(self, workflow=None, executable=None, outputDir=None, namespace=None, version=None,\
			parentJobLs=None, extraDependentInputLs=None):
		"""
		2013.2.11, wrapper around yh_pegasus.addMkDirJob
			i.e. 
			simulateOutputDirJob = self.addMkDirJob(outputDir=simulateOutputDir)
		"""
		from pymodule.pegasus import yh_pegasus
		if workflow is None:
			workflow=self
		if namespace is None:
			namespace = workflow.namespace
		if version is None:
			version = workflow.version
		if executable is None:
			executable = workflow.mkdirWrap
		
		return yh_pegasus.addMkDirJob(workflow=workflow, mkdir=executable, outputDir=outputDir, namespace=namespace, \
							version=version,\
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs)
	
	def addPlinkJob(self, workflow=None, executable=None, inputFileList=None, tpedFile=None, tfamFile=None,\
				pedFile=None, famFile=None, mapFile=None, bedFile=None, bimFile=None,\
				inputFnamePrefix=None, inputOption='--file', \
				outputFnamePrefix=None, outputOption='--out',\
				makeBED=False, calculateMendelError=False, checkSex=False, \
				LDPruneWindowSize=100, LDPruneWindowShiftSize=5, LDPruneByPairwiseR2=False, LDPruneMinR2=0.1,\
				LDPruneByRegression=False, LDPruneMinVarianceInflationFactor=2,\
				estimatePairwiseGenomeWideIBD=False, estimatePairwiseGenomeWideIBDFreqFile=None, \
				extractSNPFile=None, recodeOutput=False, estimateAlleFrequency=False, \
				mergeListFile=None,\
				parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
				extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.8.28
			add argument
				estimateAlleFrequency, estimate frequency of input file. "--nonfounders" could be added as well.
				estimatePairwiseGenomeWideIBDFreqFile, is the file from which IBD check could draw frequency (rather than estimate from founders)
		2012.8.9
			inputFileList is a list of pegasus Files (.ped, .fam, or .tped, .tfam, etc.) or could be supplied individually.
			
			inputOption could be, "--file" for .ped .map ; "--tfile" for .tped, .tfam; or '--bfile' for .bed, .fam, .bim
		
			if extractSNPFile or mergeListFile is given, either recodeOutput or makeBED have to be on. otherwise, no output.
			http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml
			
			i.e.
			
				bedFnamePrefix = os.path.join(topOutputDir, '%s_bed'%(commonPrefix))
				convertSingleTPED2BEDJob = self.addPlinkJob(executable=self.plink, inputFileList=[], 
									tpedFile=modifyTPEDJob.output, tfamFile=tfamJob.tfamFile,\
					outputFnamePrefix=bedFnamePrefix, outputOption='--out',\
					makeBED=True, \
					extraDependentInputLs=None, transferOutput=transferOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs = convertSingleTPED2BEDParentJobLs)
					
				
				convertMergedTPED2BEDJob = self.addPlinkJob(executable=self.plink, inputFileList=[tpedFileMergeJob.output, tfamJob.tfamFile], \
								inputFnamePrefix=mergedPlinkFnamePrefix, inputOption='--tfile', \
				outputFnamePrefix=mergedPlinkBEDFnamePrefix, outputOption='--out',\
				makeBED=True, \
				extraDependentInputLs=None, transferOutput=transferOutput, \
				extraArguments=None, job_max_memory=2000, parentJobLs=[mergedOutputDirJob, tpedFileMergeJob, tfamJob])
			
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if inputFileList:
			extraDependentInputLs.extend(inputFileList)
		
		extraArgumentList = []
		extraOutputLs = []
		key2ObjectForJob = {}
		if inputOption and inputFnamePrefix:
			extraArgumentList.extend([inputOption, inputFnamePrefix])
		if tpedFile:
			extraDependentInputLs.append(tpedFile)
			extraArgumentList.extend(["--tped", tpedFile])
		if tfamFile:
			extraDependentInputLs.append(tfamFile)
			extraArgumentList.extend(["--tfam", tfamFile])
		if pedFile:
			extraDependentInputLs.append(pedFile)
			extraArgumentList.extend(["--ped", pedFile])
		if famFile:
			extraDependentInputLs.append(famFile)
			extraArgumentList.extend(["--fam", famFile])
		if mapFile:
			extraDependentInputLs.append(mapFile)
			extraArgumentList.extend(["--map", mapFile])
		if bedFile:
			extraDependentInputLs.append(bedFile)
			extraArgumentList.extend(["--bed", bedFile])
		if bimFile:
			extraDependentInputLs.append(bimFile)
			extraArgumentList.extend(["--bim", bimFile])
		
		if outputFnamePrefix and outputOption:
			extraArgumentList.extend([outputOption, outputFnamePrefix])
		else:
			outputFnamePrefix = 'plink'
		
		logFile = File('%s.log'%(outputFnamePrefix))	#2012.8.10 left in the folder dying
		extraOutputLs.append(logFile)
		
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 
		if makeBED:
			extraArgumentList.append('--make-bed')
			suffixAndNameTupleList.extend([['.bed',], ('.fam',), ['.bim',]])		#, binary map file, is excluded for now
		if calculateMendelError:
			extraArgumentList.append('--mendel')
			suffixAndNameTupleList.extend([('.mendel',), ('.imendel',), ('.fmendel',), ('.lmendel',)])
			#its output is not tab-delimited. rather it's space (multi) delimited.
		if checkSex:
			extraArgumentList.append('--check-sex')
			suffixAndNameTupleList.extend([('.sexcheck',), ('.hh', )])	#.sexcheck file is accessible as job.sexcheckFile.
				#.hh is heterozygous haplotype genotypes 
		if LDPruneByPairwiseR2:
			extraArgumentList.append('--indep-pairwise %s %s %s'%(LDPruneWindowSize, LDPruneWindowShiftSize, LDPruneMinR2))
			suffixAndNameTupleList.extend([('.prune.in',), ('.prune.out',)])	#".prune.in" is accessible as job.prune_inFile
		if LDPruneByRegression:
			extraArgumentList.append('--indep %s %s %s'%(LDPruneWindowSize, LDPruneWindowShiftSize, LDPruneMinVarianceInflationFactor))
			suffixAndNameTupleList.extend([('.prune.in',), ('.prune.out',)])	#".prune.in" is accessible as job.prune_inFile
		if estimatePairwiseGenomeWideIBD:
			extraArgumentList.append('--genome')
			suffixAndNameTupleList.extend([('.genome',)])	#.genome is accessible as job.genomeFile
			if estimatePairwiseGenomeWideIBDFreqFile:	#2012.8.28
				extraArgumentList.extend(['--read-freq', estimatePairwiseGenomeWideIBDFreqFile])
				extraDependentInputLs.append(estimatePairwiseGenomeWideIBDFreqFile)
		if extractSNPFile:
			extraArgumentList.extend(['--extract', extractSNPFile])
			extraDependentInputLs.append(extractSNPFile)
		if recodeOutput:
			extraArgumentList.extend(['--recode',])
			suffixAndNameTupleList.extend([('.ped',), ('.map',)])
		if estimateAlleFrequency:	#2012.8.28
			extraArgumentList.append('--freq')
			suffixAndNameTupleList.extend([('.frq',)])
		
		if mergeListFile:
			extraArgumentList.extend(['--merge-list', mergeListFile])
			extraDependentInputLs.append(mergeListFile)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, suffixAndNameTupleList=suffixAndNameTupleList, \
													extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
		
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
		return job
	
	def setupMoreOutputAccordingToSuffixAndNameTupleList(self, outputFnamePrefix=None, suffixAndNameTupleList=None, extraOutputLs=None, key2ObjectForJob=None):
		"""
		2012.8.16
			split from addPlinkJob()
		"""
		for suffixNameTuple in suffixAndNameTupleList:
			if len(suffixNameTuple)==1:
				suffix = suffixNameTuple[0]
				name = suffix[1:].replace('.', '_')	#replace dot with underscore. as dot is used to access method/attribute of python object
				# i.e. ".prune.in" is accessible as job.prune_inFile
			elif len(suffixNameTuple)>=2:
				suffix, name = suffixNameTuple[:2]
			outputFile = File('%s%s'%(outputFnamePrefix, suffix))
			extraOutputLs.append(outputFile)
			key2ObjectForJob['%sFile'%(name)] = outputFile