#!/usr/bin/env python
"""
2012.5.23
	a common class for pegasus workflows
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from Pegasus.DAX3 import *


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
		
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		#change the workflow name to reflect the output filename
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		self.name = workflowName
		
		for pathName in self.pathToInsertHomePathList:
			setattr(self, pathName, self.insertHomePath(getattr(self, pathName, None), self.home_path))
		#self.pymodulePath = self.insertHomePath(self.pymodulePath, self.home_path)
		#self.vervetSrcPath =  self.insertHomePath(self.vervetSrcPath, self.home_path)
		
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		self.architecture = "x86_64"
		self.operatingSystem = "linux"
		self.namespace = "workflow"
		self.version="1.0"
	
	def processListArguments(self, listArgumentName_data_type_ls=None, emptyContent=[]):
		"""
		2012.8.15
		"""
		listArgumentName2hasContent = {}
		for listArgumentName, data_type in listArgumentName_data_type_ls:
			listArgumentValue = getattr(self, listArgumentName, None)
			if listArgumentValue:
				setattr(self, listArgumentName, getListOutOfStr(listArgumentValue, data_type=data_type))
				listArgumentName2hasContent[listArgumentName]=True
			else:
				setattr(self, listArgumentName, emptyContent)
				listArgumentName2hasContent[listArgumentName]=False
		return listArgumentName2hasContent
	
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
		2012.5.23 copied from AbstractNGSWorkflow
		2012.1.9
		"""
		if inputPath.find('%s')!=-1:
			inputPath = inputPath%home_path
		return inputPath
	
	def registerJars(self):
		"""
		2012.5.23
			register jars to be used in the worflow
		"""
		pass
	
	def registerCustomJars(self):
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
		
		mergeSameHeaderTablesIntoOne = Executable(namespace=namespace, name="mergeSameHeaderTablesIntoOne", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		mergeSameHeaderTablesIntoOne.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/MergeSameHeaderTablesIntoOne.py"), site_handler))
		executableClusterSizeMultiplierList.append((mergeSameHeaderTablesIntoOne,0))
		
		ReduceMatrixByChosenColumn = Executable(namespace=namespace, name="ReduceMatrixByChosenColumn", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		ReduceMatrixByChosenColumn.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/ReduceMatrixByChosenColumn.py"), site_handler))
		executableClusterSizeMultiplierList.append((ReduceMatrixByChosenColumn,0))
		
		ReduceMatrixBySumSameKeyColsAndThenDivide = Executable(namespace=namespace, name="ReduceMatrixBySumSameKeyColsAndThenDivide", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		ReduceMatrixBySumSameKeyColsAndThenDivide.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/ReduceMatrixBySumSameKeyColsAndThenDivide.py"), \
											site_handler))
		executableClusterSizeMultiplierList.append((ReduceMatrixBySumSameKeyColsAndThenDivide,0))
		
		
		ReduceMatrixByAverageColumnsWithSameKey = Executable(namespace=namespace, name="ReduceMatrixByAverageColumnsWithSameKey", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		ReduceMatrixByAverageColumnsWithSameKey.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/ReduceMatrixByAverageColumnsWithSameKey.py"), \
												site_handler))
		executableClusterSizeMultiplierList.append((ReduceMatrixByAverageColumnsWithSameKey,0))
		
		
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
		SelectLineBlockFromFile.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/SelectLineBlockFromFile.py"), \
										site_handler))
		executableClusterSizeMultiplierList.append((SelectLineBlockFromFile, 1))
		
		AbstractPlot =  Executable(namespace=namespace, name="AbstractPlot", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		AbstractPlot.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "plot/AbstractPlot.py"), site_handler))
		executableClusterSizeMultiplierList.append((AbstractPlot, 0))
		
		PlotLD = Executable(namespace=namespace, name="PlotLD", version=version, os=operatingSystem, arch=architecture, installed=True)
		PlotLD.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "plot/PlotLD.py"), site_handler))
		executableClusterSizeMultiplierList.append((PlotLD, 0))
		
		PlotXYAsBarChart = Executable(namespace=namespace, name="PlotXYAsBarChart", version=version, os=operatingSystem, arch=architecture, installed=True)
		PlotXYAsBarChart.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "plot/PlotXYAsBarChart.py"), site_handler))
		executableClusterSizeMultiplierList.append((PlotXYAsBarChart, 0))
		
		DrawHistogram = Executable(namespace=namespace, name="DrawHistogram", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		DrawHistogram.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "plot/DrawHistogram.py"), site_handler))
		executableClusterSizeMultiplierList.append((DrawHistogram, 0))
		
		#2012.8.13 SelectRowsFromMatrix is a derivative of AbstractMatrixFileWalker, so use addAbstractMatrixFileWalkerJob()
		SelectRowsFromMatrix = Executable(namespace=namespace, name="SelectRowsFromMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		SelectRowsFromMatrix.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/SelectRowsFromMatrix.py"), site_handler))
		executableClusterSizeMultiplierList.append((SelectRowsFromMatrix, 1))
		
		#2012.8.15 ancestor of SelectRowsFromMatrix, 
		AbstractMatrixFileWalker  = Executable(namespace=namespace, name="AbstractMatrixFileWalker", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		AbstractMatrixFileWalker.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "AbstractMatrixFileWalker.py"), site_handler))
		executableClusterSizeMultiplierList.append((AbstractMatrixFileWalker, 1))
		
		#2012.8.13
		OutputVCFSiteGap = Executable(namespace=namespace, name="OutputVCFSiteGap", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		OutputVCFSiteGap.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/OutputVCFSiteGap.py"), site_handler))
		executableClusterSizeMultiplierList.append((OutputVCFSiteGap, 1))
		
		java = Executable(namespace=namespace, name="java", version=version, os=operatingSystem, arch=architecture, installed=True)
		java.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((java, 1))
		
		plink =  Executable(namespace=namespace, name="plink", \
						version=version, os=operatingSystem, arch=architecture, installed=True)
		plink.addPFN(PFN("file://" + self.plinkPath, site_handler))
		executableClusterSizeMultiplierList.append((plink, 1))
		
		ConvertBjarniSNPFormat2Yu = Executable(namespace=namespace, name="ConvertBjarniSNPFormat2Yu", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertBjarniSNPFormat2Yu.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "ConvertBjarniSNPFormat2Yu.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertBjarniSNPFormat2Yu, 1))
		
		ConvertVCF2BjarniFormat = Executable(namespace=namespace, name="ConvertVCF2BjarniFormat", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertVCF2BjarniFormat.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "pegasus/mapper/ConvertVCF2BjarniFormat.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertVCF2BjarniFormat, 1))
		
		
		ConvertYuSNPFormat2Bjarni = Executable(namespace=namespace, name="ConvertYuSNPFormat2Bjarni", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertYuSNPFormat2Bjarni.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "ConvertYuSNPFormat2Bjarni.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2Bjarni, 1))
		
		ConvertYuSNPFormat2EigenStrat = Executable(namespace=namespace, name="ConvertYuSNPFormat2EigenStrat", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertYuSNPFormat2EigenStrat.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "ConvertYuSNPFormat2EigenStrat.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2EigenStrat, 1))
		
		ConvertYuSNPFormat2TPED_TFAM = Executable(namespace=namespace, name="ConvertYuSNPFormat2TPED_TFAM", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		ConvertYuSNPFormat2TPED_TFAM.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "ConvertYuSNPFormat2TPED_TFAM.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2TPED_TFAM, 1))
		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def addExecutableAndAssignProperClusterSize(self, executableClusterSizeMultiplierList=[], defaultClustersSize=None):
		"""
		2012.8.9
			
		"""
		if defaultClustersSize is None:
			defaultClustersSize = self.clusters_size
		for executableAndclusterSizeMultipler in executableClusterSizeMultiplierList:
			executable = executableAndclusterSizeMultipler[0]
			if len(executableAndclusterSizeMultipler)==1:
				clusterSizeMultipler = 1
			elif len(executableAndclusterSizeMultipler)>1:
				clusterSizeMultipler = executableAndclusterSizeMultipler[1]
			clusterSize = int(defaultClustersSize*clusterSizeMultipler)
			if clusterSize>1:
				executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusterSize))
			self.addExecutable(executable)
			setattr(self, executable.name, executable)
		
		
	
	
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
	
	def registerAllInputFiles(self, workflow=None, inputFnameLs=[], input_site_handler=None, pegasusFolderName=''):
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
	
	def registerOneInputFile(self, workflow=None, inputFname=None, folderName=""):
		"""
		2012.3.22
			add abspath attribute to file.
		2012.3.1
			add argument folderName, which will put the file in specific pegasus workflow folder
		2011.12.21
		"""
		file = File(os.path.join(folderName, os.path.basename(inputFname)))
		file.abspath = os.path.abspath(inputFname)
		file.absPath = os.path.abspath(inputFname)
		file.addPFN(PFN("file://" + file.abspath, self.input_site_handler))
		self.addFile(file)
		return file
	
	def addStatMergeJob(self, workflow=None, statMergeProgram=None, outputF=None, \
					parentJobLs=None, \
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
		extraOutputLs = []
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
							parentJobLs=[], \
							namespace=None, version=None, extraDependentInputLs=[]):
		"""
		2011-11-28
			moved from CalculateVCFStatPipeline.py
		"""
		statMergeJob.addArguments(inputF)
		statMergeJob.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		for input in extraDependentInputLs:
			if input:
				statMergeJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=statMergeJob)
	
	def addGenericJob(self, workflow=None, executable=None, inputFile=None, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", \
					parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					extraArguments=None, extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
					key2ObjectForJob=None, **keywords):
		"""
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
		if inputFile:
			if inputArgumentOption:
				job.addArguments(inputArgumentOption)
			job.addArguments(inputFile)
			job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
			job.input = inputFile
			job.inputLs.append(inputFile)
		if outputFile:
			if outputArgumentOption:
				job.addArguments(outputArgumentOption)
			job.addArguments(outputFile)
			if transferOutput is not None:	#2012.8.17
				job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			job.output = outputFile
			job.outputLs.append(outputFile)
		if extraArgumentList:
			job.addArguments(*extraArgumentList)
		
		if extraArguments:
			job.addArguments(extraArguments)
		
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel)
		workflow.addJob(job)
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					workflow.depends(parent=parentJob, child=job)
		if extraDependentInputLs:
			for input in extraDependentInputLs:
				if input:
					job.uses(input, transfer=True, register=True, link=Link.INPUT)
					job.inputLs.append(input)
		if extraOutputLs:
			for output in extraOutputLs:
				if output:
					job.outputLs.append(output)
					if transferOutput is not None:	#2012.8.17
						job.uses(output, transfer=transferOutput, register=True, link=Link.OUTPUT)
		if key2ObjectForJob:
			for key, object in key2ObjectForJob.iteritems():
				setattr(job, key, object)	#key should be a string.
		#2012.7.28 if job.output is not set, set it to the 1st entry of job.outputLs
		if getattr(job, 'output', None) is None and job.outputLs:
			job.output = job.outputLs[0]
		if getattr(job, 'input', None) is None and job.inputLs:
			job.input = job.inputLs[0]
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
						outputDirPrefix="", **keywords):
		"""
		2012.8.2 bugfix.
		2012.7.19
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding gzip jobs for %s input job data ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sGzip"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for jobData in inputData.jobDataLs:
			for inputF in jobData.fileList:
				inputFBaseName = os.path.basename(inputF.name)
				outputF = File(os.path.join(topOutputDir, '%s.gz'%(inputFBaseName)))
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
				returnData.jobDataLs.append(PassingData(jobLs=[job], vcfFile=None, \
										fileList=[outputF]))
				no_of_jobs += 1
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		return returnData
	
	def addAbstractMapperLikeJob(self, workflow=None, executable=None, \
					inputVCF=None, inputF=None, outputF=None, \
					parentJobLs=[], namespace=None, version=None, transferOutput=True, job_max_memory=200,\
					extraArguments=None, extraDependentInputLs=[]):
		"""
		2012.7.19
			moved from AbstractNGSWorkflow to here.
			add argument inputF. inputVCF is not generic enough.
		2012.5.11
		"""
		if inputF is None:	#2012.7.19
			inputF = inputVCF
		#2011-9-22 union of all samtools intervals for one contig
		job = Job(namespace=getattr(self, 'namespace', namespace), name=executable.name, \
						version=getattr(self, 'version', version))
		job.addArguments("-i", inputF, "-o", outputF)
		
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		for input in extraDependentInputLs:
			if input:
				job.uses(input, transfer=True, register=True, link=Link.INPUT)
		job.output = outputF
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=job)
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
					whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
					logWhichColumn=False, positiveLog=False, valueForNonPositiveYValue=50, \
					xColumnPlotLabel=None, chrLengthColumnHeader=None, chrColumnHeader=None, \
					minChrLength=1000000, xColumnHeader=None, pos2ColumnHeader=None, minNoOfTotal=100,\
					figureDPI=300, formatString='.', ylim_type=2, samplingRate=0.0001,  need_svg=False, logCount=False, \
					parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
		"""
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
		return self.addAbstractPlotJob(workflow=workflow, executable=executable, inputFileList=inputFileList, \
							inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
							whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=whichColumnPlotLabel, logWhichColumn=logWhichColumn, \
							positiveLog=positiveLog, valueForNonPositiveYValue=valueForNonPositiveYValue, \
							xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel, \
							minNoOfTotal=minNoOfTotal, \
							figureDPI=figureDPI, formatString=formatString, ylim_type=ylim_type, samplingRate=samplingRate, need_svg=need_svg, \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraArgumentList=extraArgumentList, extraArguments=extraArguments, transferOutput=transferOutput, \
							job_max_memory=job_max_memory, \
							**keywords)
		
		
	
	def addAbstractPlotJob(self, workflow=None, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
					logWhichColumn=False, positiveLog=False, valueForNonPositiveYValue=50, \
					xColumnHeader=None, xColumnPlotLabel=None, \
					minNoOfTotal=100,\
					figureDPI=300, formatString='.', ylim_type=2, samplingRate=0.001, need_svg=False, \
					parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
		"""
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		2012.8.2
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
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
			('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
				only effective when logWhichColumn is toggled. '],\
			('valueForNonPositiveYValue', 1, float): [50, '', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
					what yValue should be.'],\
			('xColumnHeader', 1, ): ['', 'l', 1, 'header of the x-axis data column, ' ],\
			('xColumnPlotLabel', 0, ): ['', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
			('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
			
		"""
		if executable is None:
			executable = self.AbstractPlot
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if inputFileList:
			extraDependentInputLs.extend(inputFileList)
		if extraArgumentList is None:
			extraArgumentList = []
		if xColumnHeader:
			extraArgumentList.append('--xColumnHeader %s'%(xColumnHeader))
		extraOutputLs = []
		key2ObjectForJob = {}
		if outputFnamePrefix:
			extraArgumentList.append('--outputFnamePrefix %s'%(outputFnamePrefix))
			if outputFile is None:
				extraOutputLs.append(File('%s.png'%(outputFnamePrefix)))
				if need_svg:
					extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
		if minNoOfTotal:
			extraArgumentList.append('--minNoOfTotal %s'%(minNoOfTotal))
		if figureDPI:
			extraArgumentList.append('--figureDPI %s'%(figureDPI))
		if formatString:
			extraArgumentList.append('--formatString %s'%(formatString))
		if ylim_type:
			extraArgumentList.append('--ylim_type %s'%(ylim_type))
		if samplingRate is not None:
			extraArgumentList.append('--samplingRate %s'%(samplingRate))
		if whichColumnHeader:
			extraArgumentList.append("--whichColumnHeader %s"%(whichColumnHeader))
		if whichColumn:
			extraArgumentList.append("--whichColumn %s"%(whichColumn))
		if logWhichColumn:
			extraArgumentList.append('--logWhichColumn')
			if positiveLog:
				extraArgumentList.append('--positiveLog')
		if whichColumnPlotLabel:
			extraArgumentList.append("--whichColumnPlotLabel %s"%(whichColumnPlotLabel))
		if xColumnPlotLabel:
			extraArgumentList.append("--xColumnPlotLabel %s"%(xColumnPlotLabel))
		if valueForNonPositiveYValue:
			extraArgumentList.append("--valueForNonPositiveYValue %s"%(valueForNonPositiveYValue))
		if need_svg:
			extraArgumentList.append('--need_svg')
			if not outputFnamePrefix:
				outputFnamePrefix = os.path.splitext(outputFile.name)[0]	#2012.8.20 bugfix.
			extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		#add all input files to the last
		if inputFileList:
			extraArgumentList.extend(inputFileList)
		job= self.addGenericJob(workflow=workflow, executable=executable, inputFile=inputFile, outputFile=outputFile, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
		return job
	
	def addAbstractMatrixFileWalkerJob(self, workflow=None, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, \
					logWhichColumn=False, positiveLog=False, valueForNonPositiveYValue=50, \
					minNoOfTotal=10,\
					samplingRate=0.001, \
					parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
		"""
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
							whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=None, logWhichColumn=logWhichColumn, \
							positiveLog=positiveLog, valueForNonPositiveYValue=valueForNonPositiveYValue, \
							minNoOfTotal=minNoOfTotal, \
							figureDPI=None, formatString=None, ylim_type=None, samplingRate=samplingRate, need_svg=False, \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
							**keywords)
	
	def addDrawHistogramJob(self, workflow=None, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
					logWhichColumn=False, positiveLog=False, valueForNonPositiveYValue=50, \
					minNoOfTotal=10,\
					figureDPI=100, formatString='.', ylim_type=2, samplingRate=0.001, need_svg=False, \
					logCount=False,\
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
							logWhichColumn=logWhichColumn, \
							positiveLog=positiveLog, valueForNonPositiveYValue=valueForNonPositiveYValue, \
							minNoOfTotal=minNoOfTotal, \
							figureDPI=figureDPI, formatString=formatString, ylim_type=ylim_type, \
							samplingRate=samplingRate, need_svg=need_svg, \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
							**keywords)
	
	def addPlinkJob(self, workflow=None, executable=None, inputFileList=None, tpedFile=None, tfamFile=None,\
				pedFile=None, famFile=None, mapFile=None, bedFile=None, bimFile=None,\
				inputFnamePrefix=None, inputOption='--file', \
				outputFnamePrefix=None, outputOption='--out',\
				makeBED=False, calculateMendelError=False, checkSex=False, \
				LDPruneWindowSize=100, LDPruneWindowShiftSize=5, LDPruneByPairwiseR2=False, LDPruneMinR2=0.1,\
				LDPruneByRegression=False, LDPruneMinVarianceInflationFactor=2,\
				estimatePairwiseGenomeWideIBD=False,\
				extractSNPFile=None, recodeOutput=False,\
				mergeListFile=None,\
				parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
				extraArguments=None, job_max_memory=2000, **keywords):
		"""
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
		if extractSNPFile:
			extraArgumentList.extend(['--extract', extractSNPFile])
			extraDependentInputLs.append(extractSNPFile)
		if recodeOutput:
			extraArgumentList.extend(['--recode',])
			suffixAndNameTupleList.extend([('.ped',), ('.map',)])
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