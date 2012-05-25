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
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("javaPath", 1, ): ["/usr/bin/java", 'J', 1, 'java interpreter binary'],\
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
		
		self.javaPath =  self.insertHomePath(self.javaPath, self.home_path)
		self.pymodulePath = self.insertHomePath(self.pymodulePath, self.home_path)
		self.vervetSrcPath =  self.insertHomePath(self.vervetSrcPath, self.home_path)
		
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		self.architecture = "x86_64"
		self.operatingSystem = "linux"
		self.namespace = "workflow"
		self.version="1.0"
		
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
	
	def registerExecutables(self):
		"""
		2012.1.9 a symlink to registerCommonExecutables()
		"""
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableList = []
		
		mergeSameHeaderTablesIntoOne = Executable(namespace=namespace, name="mergeSameHeaderTablesIntoOne", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		mergeSameHeaderTablesIntoOne.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/MergeSameHeaderTablesIntoOne.py"), site_handler))
		executableList.append(mergeSameHeaderTablesIntoOne)
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "shell/mkdirWrap.sh"), site_handler))
		executableList.append(mkdirWrap)
		
		#mv to rename files and move them
		mv = Executable(namespace=namespace, name="mv", version=version, os=operatingSystem, arch=architecture, installed=True)
		mv.addPFN(PFN("file://" + "/bin/mv", site_handler))
		executableList.append(mv)
		
		for executable in executableList:
			executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			self.addExecutable(executable)
			setattr(self, executable.name, executable)
		
	
	def registerCommonExecutables(self):
		"""
		2011-11-22
		"""
		pass
	
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
		file.addPFN(PFN("file://" + file.abspath, self.input_site_handler))
		self.addFile(file)
		return file
	
	def addStatMergeJob(self, workflow=None, statMergeProgram=None, outputF=None, \
					parentJobLs=[], \
					extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					namespace=None, version=None, **keywords):
		"""
		2012.4.3
			make argument namespace, version optional
		2011-11-28
			moved from CalculateVCFStatPipeline.py
		2011-11-17
			add argument extraArguments
		"""
		statMergeJob = Job(namespace=getattr(self, 'namespace', namespace), name=statMergeProgram.name, \
						version=getattr(self, 'version', version))
		statMergeJob.addArguments('-o', outputF)
		if extraArguments:
			statMergeJob.addArguments(extraArguments)
		statMergeJob.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		statMergeJob.output = outputF
		self.addJob(statMergeJob)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=statMergeJob)
		for input in extraDependentInputLs:
			if input:
				statMergeJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		return statMergeJob
	
	def addInputToStatMergeJob(self, workflow=None, statMergeJob=None, inputF=None, \
							parentJobLs=[], \
							namespace=None, version=None, extraDependentInputLs=[]):
		"""
		2011-11-28
			moved from CalculateVCFStatPipeline.py
		"""
		statMergeJob.addArguments(inputF)
		statMergeJob.uses(inputF, transfer=False, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=statMergeJob)
	
	def addGenericJob(self, executable=None, inputFile=None, outputFile=None, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, extraArgumentList=[], job_max_memory=2000,  sshDBTunnel=None, \
						**keywords):
		"""
		2012.5.24
			generic job addition function for other functions to use
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		
		if inputFile:
			job.addArguments("-i", inputFile)
		if outputFile:
			job.addArguments("-o", outputFile)
		if extraArgumentList:
			job.addArguments(*extraArgumentList)
		
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		if outputFile:
			job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			job.output = outputFile
			
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel)
		self.addJob(job)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			if input:
				job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job