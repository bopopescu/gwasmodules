#!/usr/bin/env python
"""
Examples:

	%s -d /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta 
		-i ~/script/vervet/data/OphoffMethylation/DMR330K_ProbeSeq.fasta  -a 2 -l condorpool -j condorpool 
		-C 1 -o workflow/BlastDMR330K_ProbeSeqAgainst524.xml

2011-11-22
	a common class for pegasus workflows that work on NGS (next-gen sequencing) data
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0],)

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, utils
from Pegasus.DAX3 import *
from AbstractWorkflow import AbstractWorkflow

class BlastWorkflow(AbstractWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractWorkflow.option_default_dict.copy()
	option_default_dict.update({
						("inputFname", 1, ): ["", 'i', 1, 'the input fasta file'],\
						("formatdbPath", 1, ): ["/usr/bin/formatdb", 'f', 1, 'path to formatdb, index fasta database file'],\
						("blastallPath", 1, ): ["/usr/bin/blastall", 's', 1, 'path to blastall'],\
						("blockSize", 1, int): [1000, 'c', 1, 'how many sequences each blast job handles'],\
						('databaseFname', 1, ): ['', 'd', 1, 'filename of the database to blast against, must be indexed', ],\
						('minNoOfIdentities', 0, int): [None, 'm', 1, 'minimum number of identities between a query and target', ],\
						('maxNoOfMismatches', 0, int): [None, 'a', 1, 'minimum number of mismatches between a query and target', ],\
						('minIdentityPercentage', 0, float): [None, 'n', 1, 'minimum percentage of identities between a query and target', ],\
						
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.5.23
		"""
		AbstractWorkflow.__init__(self, **keywords)
	
	def getNoOfSequencesFromFasta(self, inputFastaFname=None):
		"""
		2012.5.24
		"""
		sys.stderr.write("Getting number of sequences from %s ..."%(inputFastaFname))
		inf = utils.openGzipFile(inputFastaFname)
		no_of_sequences = 0
		for line in inf:
			if line[0]=='>':
				no_of_sequences += 1
		del inf
		sys.stderr.write("%s sequences.\n"%(no_of_sequences))
		return no_of_sequences
	
	def addSplitFastaFileJob(self, executable=None, inputFile=None, outputFnamePrefix=None, \
						noOfSequencesPerSplitFile=1000, filenameSuffix="", noOfTotalSequences=1000000,\
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=500, **keywords):
		"""
		2012.5.24
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		
		noOfSplitFiles = int(math.ceil(noOfTotalSequences/float(noOfSequencesPerSplitFile)))
		suffixLength = len(repr(noOfSplitFiles))
		
		job.addArguments("-i", inputFile, "-l %s"%(noOfSequencesPerSplitFile), \
						"-O", outputFnamePrefix, '-f %s'%(filenameSuffix), '-a %s'%(suffixLength))
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		job.outputList = []
		for i in xrange(noOfSplitFiles):	#start from 0
			splitFname = utils.comeUpSplitFilename(outputFnamePrefix=outputFnamePrefix, suffixLength=suffixLength, fileOrder=i,\
											filenameSuffix=filenameSuffix)
			splitFile = File(splitFname)
			
			job.outputList.append(splitFile)
			job.uses(splitFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			if input:
				job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addBlastWrapperJob(self, executable=None, inputFile=None, outputFile=None, outputFnamePrefix=None, databaseFile=None,\
						maxNoOfMismatches=None, minNoOfIdentities=None, minIdentityPercentage=None, blastallPath=None, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.5.24
		"""
		extraArgumentList = ['-l %s'%blastallPath, '--databaseFname', databaseFile]
		if outputFnamePrefix:
			extraArgumentList.append("--outputFnamePrefix %s "%(outputFnamePrefix))
		if maxNoOfMismatches:
			extraArgumentList.append('--maxNoOfMismatches %s'%maxNoOfMismatches)
		if minNoOfIdentities:
			extraArgumentList.append("--minNoOfIdentities %s"%(minNoOfIdentities))
		if minIdentityPercentage:
			extraArgumentList.append("--minIdentityPercentage %s"%(minIdentityPercentage))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		return self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory)
		
	def addJobs(self, inputData=None, outputDirPrefix="", ntDatabaseFileList=None, noOfTotalSequences=None, \
			transferOutput=True, makeBlastDBJob=None):
		"""
		2012.5.24
		"""
		
		sys.stderr.write("Adding blast jobs for %s input ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		topOutputDir = "%sBlast"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(self, mkdir=self.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		allBlastResultFile = File(os.path.join(topOutputDir, 'blast.tsv'))
		allBlastMergeJob = self.addStatMergeJob(statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
							outputF=allBlastResultFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob])
		no_of_jobs += 1
		
		ntDatabaseFile = ntDatabaseFileList[0]
		returnData = PassingData()
		returnData.jobDataLs = []
		
		for jobData in inputData.jobDataLs:
			inputF = jobData.output
			outputFnamePrefix = os.path.join(topOutputDir, os.path.splitext(os.path.basename(inputF.name))[0])
			
			splitFastaJob = self.addSplitFastaFileJob(executable=self.SplitFastaFile, inputFile=inputF, outputFnamePrefix=outputFnamePrefix, \
						noOfSequencesPerSplitFile=self.blockSize, filenameSuffix=".fasta", noOfTotalSequences=noOfTotalSequences,\
						parentJobLs=jobData.jobLs + [topOutputDirJob], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=500)
			no_of_jobs += 1
			for splitFastaOutput in splitFastaJob.outputList:
				outputFile = File('%s.tsv'%(splitFastaOutput.name))
				blastJob = self.addBlastWrapperJob(executable=self.BlastWrapper, inputFile=splitFastaOutput, outputFile=outputFile, \
								outputFnamePrefix=splitFastaOutput.name , databaseFile=ntDatabaseFile,\
								maxNoOfMismatches=self.maxNoOfMismatches, minNoOfIdentities=self.minNoOfIdentities, \
								minIdentityPercentage=self.minIdentityPercentage, blastallPath=self.blastallPath, \
								parentJobLs=[splitFastaJob, makeBlastDBJob], extraDependentInputLs=ntDatabaseFileList, transferOutput=False, \
								extraArguments=None, job_max_memory=1000)
				
				#add output to some reduce job
				self.addInputToStatMergeJob(statMergeJob=allBlastMergeJob, \
								inputF=blastJob.output, parentJobLs=[blastJob])
				no_of_jobs += 1
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		#include the tfam (outputList[1]) into the fileList
		returnData.jobDataLs.append(PassingData(jobLs=[allBlastMergeJob], file=allBlastResultFile, \
											fileList=[allBlastResultFile]))
		return returnData
	
				
		
	def registerCustomExecutables(self, workflow=None, **keywords):
		"""
		2012.5.23
		"""
		if workflow is None:
			workflow = self
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)		
		blastall = Executable(namespace=namespace, name="blastall", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		blastall.addPFN(PFN("file://" + os.path.join(self.blastallPath), site_handler))
		executableClusterSizeMultiplierList.append((blastall, 0.2))
		
		formatdb = Executable(namespace=namespace, name="formatdb", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		formatdb.addPFN(PFN("file://" + os.path.join(self.formatdbPath), site_handler))
		executableClusterSizeMultiplierList.append((formatdb, 0))
		
		BlastWrapper = Executable(namespace=namespace, name="BlastWrapper", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		BlastWrapper.addPFN(PFN("file://" + os.path.join(self.pymodulePath, 'pegasus/mapper/alignment/BlastWrapper.py'), site_handler))
		executableClusterSizeMultiplierList.append((BlastWrapper, 0.1))
		
		SplitFastaFile = Executable(namespace=namespace, name="SplitFastaFile", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		SplitFastaFile.addPFN(PFN("file://" + os.path.join(self.pymodulePath, 'pegasus/mapper/splitter/SplitFastaFile.py'), site_handler))
		executableClusterSizeMultiplierList.append((SplitFastaFile, 0.1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)

	
	def addMakeBlastDBJob(self, executable=None, inputFile=None, \
						parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=500, **keywords):
		"""
		2012.10.9 use addGenericJob() instead
		2012.5.24
			untested
		"""
		extraOutputLs = []
		for suffix in ['.nin', '.nhr', '.nsq']:	#start from 0
			dbIndexFile = File('%s%s'%(inputFile.name, suffix))
			extraOutputLs.append(dbIndexFile)
		
		extraArgumentList = ["-p F"]
		job = self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
						key2ObjectForJob=None,\
						job_max_memory=job_max_memory)
		return job
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables()
		
		
		inputData = PassingData(jobDataLs = [])
		inputFile = self.registerOneInputFile(inputFname=self.inputFname, folderName=self.pegasusFolderName)
		inputData.jobDataLs.append(PassingData(output=inputFile, jobLs=[]))
		noOfTotalSequences= self.getNoOfSequencesFromFasta(inputFastaFname=self.inputFname)
		
		ntDatabaseFileList = self.registerBlastNucleotideDatabaseFile(ntDatabaseFname=self.databaseFname, \
																	input_site_handler=self.input_site_handler)
		ntDatabaseFile = ntDatabaseFileList[0]

		if len(ntDatabaseFileList)<4:	#some nt-database index file is missing
			sys.stderr.write("Adding blast-db-making job...")
			makeBlastDBJob = self.addMakeBlastDBJob(executable=self.formatdb,\
												inputFile=ntDatabaseFile, transferOutput=True)
			#add the index files to the ntDatabaseFileList
			ntDatabaseFileList = [ntDatabaseFile] + makeBlastDBJob.outputList
			sys.stderr.write(".\n")
		else:
			makeBlastDBJob = None
		
		self.addJobs(inputData=inputData, outputDirPrefix=self.pegasusFolderName, ntDatabaseFileList=ntDatabaseFileList, \
					noOfTotalSequences=noOfTotalSequences, \
					transferOutput=True, makeBlastDBJob=makeBlastDBJob)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)

if __name__ == '__main__':
	main_class = BlastWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
