#!/usr/bin/env python
"""
2011-11-22
	a common class for pegasus workflows that work on NGS (next-gen sequencing) data
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *


class AbstractNGSWorkflow(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/gatk/dist", '', 1, 'GATK folder containing its jar binaries'],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("javaPath", 1, ): ["/usr/bin/java", 'J', 1, 'java interpreter binary'],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("localDataDir", 0, ): ["", 'D', 1, 'localDataDir should contain same files as dataDir but accessible locally.\
									If not given, use the default stored in db. This argument is used to find all input files available.'],\
						
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('clusters_size', 1, int):[5, 'C', 1, 'For short jobs that will be clustered, how many of them should be clustered int one'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.samtools_path = self.samtools_path%self.home_path
		self.picard_path = self.picard_path%self.home_path
		self.gatk_path = self.gatk_path%self.home_path
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
	def registerJars(self, workflow, ):
		"""
		2011-11-22
			register jars to be used in the worflow
		"""
		site_handler = self.site_handler
		#add the MergeSamFiles.jar file into workflow
		abs_path = os.path.join(self.picard_path, 'MergeSamFiles.jar')
		mergeSamFilesJar = File(abs_path)	#using abs_path avoids add this jar to every job as Link.INPUT
		mergeSamFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(mergeSamFilesJar)
		workflow.mergeSamFilesJar = mergeSamFilesJar
		
		abs_path = os.path.join(self.picard_path, 'BuildBamIndex.jar')
		BuildBamIndexFilesJar = File(abs_path)
		BuildBamIndexFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(BuildBamIndexFilesJar)
		workflow.BuildBamIndexFilesJar = BuildBamIndexFilesJar
		
		abs_path = os.path.join(self.gatk_path, 'GenomeAnalysisTK.jar')
		genomeAnalysisTKJar = File(abs_path)
		genomeAnalysisTKJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(genomeAnalysisTKJar)
		workflow.genomeAnalysisTKJar = genomeAnalysisTKJar
		
		abs_path = os.path.join(self.picard_path, 'CreateSequenceDictionary.jar')
		createSequenceDictionaryJar = File(abs_path)
		createSequenceDictionaryJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(createSequenceDictionaryJar)
		workflow.createSequenceDictionaryJar = createSequenceDictionaryJar
		
		abs_path = os.path.join(self.picard_path, 'AddOrReplaceReadGroupsAndCleanSQHeader.jar')
		addOrReplaceReadGroupsAndCleanSQHeaderJar = File(abs_path)
		addOrReplaceReadGroupsAndCleanSQHeaderJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(addOrReplaceReadGroupsAndCleanSQHeaderJar)
		workflow.addOrReplaceReadGroupsAndCleanSQHeaderJar = addOrReplaceReadGroupsAndCleanSQHeaderJar
		
		abs_path = os.path.join(self.picard_path, 'AddOrReplaceReadGroups.jar')
		addOrReplaceReadGroupsJar = File(abs_path)
		addOrReplaceReadGroupsJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(addOrReplaceReadGroupsJar)
		workflow.addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar
		
		abs_path = os.path.join(self.picard_path, 'MarkDuplicates.jar')
		MarkDuplicatesJar = File(abs_path)
		MarkDuplicatesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(MarkDuplicatesJar)
		workflow.MarkDuplicatesJar = MarkDuplicatesJar
	
	def registerCommonExecutables(self, workflow):
		"""
		2011-11-22
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/mkdirWrap.sh"), site_handler))
		mkdirWrap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%workflow.clusters_size))
		workflow.addExecutable(mkdirWrap)
		workflow.mkdirWrap = mkdirWrap
		
		#mv to rename files and move them
		mv = Executable(namespace=namespace, name="mv", version=version, os=operatingSystem, arch=architecture, installed=True)
		mv.addPFN(PFN("file://" + "/bin/mv", site_handler))
		mv.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mv)
		workflow.mv = mv
		
		selectAndSplit = Executable(namespace=namespace, name="SelectAndSplitAlignment", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		selectAndSplit.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "SelectAndSplitAlignment.py"), site_handler))
		workflow.addExecutable(selectAndSplit)
		workflow.selectAndSplit = selectAndSplit
		
		selectAndSplitFasta = Executable(namespace=namespace, name="SelectAndSplitFastaRecords", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		selectAndSplitFasta.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "SelectAndSplitFastaRecords.py"), site_handler))
		workflow.addExecutable(selectAndSplitFasta)
		workflow.selectAndSplitFasta = selectAndSplitFasta
		
		samtools = Executable(namespace=namespace, name="samtools", version=version, os=operatingSystem, arch=architecture, installed=True)
		samtools.addPFN(PFN("file://" + self.samtools_path, site_handler))
		workflow.addExecutable(samtools)
		workflow.samtools = samtools
		
		java = Executable(namespace=namespace, name="java", version=version, os=operatingSystem, arch=architecture, installed=True)
		java.addPFN(PFN("file://" + self.javaPath, site_handler))
		workflow.addExecutable(java)
		workflow.java = java
		
		addOrReplaceReadGroupsJava = Executable(namespace=namespace, name="addOrReplaceReadGroupsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		addOrReplaceReadGroupsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		workflow.addExecutable(addOrReplaceReadGroupsJava)
		workflow.addOrReplaceReadGroupsJava = addOrReplaceReadGroupsJava
		
		genotyperJava = Executable(namespace=namespace, name="genotyperJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		genotyperJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		workflow.addExecutable(genotyperJava)
		workflow.genotyperJava = genotyperJava
		
		BuildBamIndexFilesJava = Executable(namespace=namespace, name="BuildBamIndexFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		BuildBamIndexFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		workflow.addExecutable(BuildBamIndexFilesJava)
		workflow.BuildBamIndexFilesJava = BuildBamIndexFilesJava
		
		createSequenceDictionaryJava = Executable(namespace=namespace, name="createSequenceDictionaryJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		createSequenceDictionaryJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		workflow.addExecutable(createSequenceDictionaryJava)
		workflow.createSequenceDictionaryJava = createSequenceDictionaryJava
		
		DOCWalkerJava = Executable(namespace=namespace, name="DepthOfCoverageWalker", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		#no clusters_size for this because it could run on a whole bam for hours
		#DOCWalkerJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		DOCWalkerJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		workflow.addExecutable(DOCWalkerJava)
		workflow.DOCWalkerJava = DOCWalkerJava
		
		VariousReadCountJava = Executable(namespace=namespace, name="VariousReadCountJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		VariousReadCountJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#no clusters_size for this because it could run on a whole bam for hours
		#VariousReadCountJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(VariousReadCountJava)
		workflow.VariousReadCountJava = VariousReadCountJava
		
		
		MarkDuplicatesJava = Executable(namespace=namespace, name="MarkDuplicatesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		MarkDuplicatesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		workflow.addExecutable(MarkDuplicatesJava)
		workflow.MarkDuplicatesJava = MarkDuplicatesJava
		
		CallVariantBySamtools = Executable(namespace=namespace, name="CallVariantBySamtools", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		CallVariantBySamtools.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/CallVariantBySamtools.sh"), site_handler))
		workflow.addExecutable(CallVariantBySamtools)
		workflow.CallVariantBySamtools = CallVariantBySamtools
		
		genotypeCallByCoverage = Executable(namespace=namespace, name="GenotypeCallByCoverage", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		genotypeCallByCoverage.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "GenotypeCallByCoverage.py"), site_handler))
		workflow.addExecutable(genotypeCallByCoverage)
		workflow.genotypeCallByCoverage = genotypeCallByCoverage
		
		mergeGenotypeMatrix = Executable(namespace=namespace, name="MergeGenotypeMatrix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		mergeGenotypeMatrix.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "MergeGenotypeMatrix.py"), site_handler))
		workflow.addExecutable(mergeGenotypeMatrix)
		workflow.mergeGenotypeMatrix = mergeGenotypeMatrix
		
		bgzip_tabix = Executable(namespace=namespace, name="bgzip_tabix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		bgzip_tabix.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/bgzip_tabix.sh"), site_handler))
		bgzip_tabix.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(bgzip_tabix)
		workflow.bgzip_tabix = bgzip_tabix
		
		vcf_convert = Executable(namespace=namespace, name="vcf_convert", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		vcf_convert.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_convert.sh"), site_handler))
		vcf_convert.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(vcf_convert)
		workflow.vcf_convert = vcf_convert
		
		vcf_isec = Executable(namespace=namespace, name="vcf_isec", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_isec.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_isec.sh"), site_handler))
		vcf_isec.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(vcf_isec)
		workflow.vcf_isec = vcf_isec
		
		vcf_concat = Executable(namespace=namespace, name="vcf_concat", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_concat.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		#vcf_concat might involve a very long argument, which would be converted to *.arg and then pegasus bug
		#vcf_concat.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(vcf_concat)
		workflow.vcf_concat = vcf_concat
		
		concatGATK = Executable(namespace=namespace, name="concatGATK", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		concatGATK.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		concatGATK.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(concatGATK)
		workflow.concatGATK = concatGATK
		
		concatSamtools = Executable(namespace=namespace, name="concatSamtools", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		concatSamtools.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		concatSamtools.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(concatSamtools)
		workflow.concatSamtools = concatSamtools
		
		calcula = Executable(namespace=namespace, name="CalculatePairwiseDistanceOutOfSNPXStrainMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		calcula.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py"), site_handler))
		workflow.addExecutable(calcula)
		workflow.calcula = calcula
		
		#2011-11-28
		mergeSameHeaderTablesIntoOne = Executable(namespace=namespace, name="MergeSameHeaderTablesIntoOne", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		mergeSameHeaderTablesIntoOne.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/MergeSameHeaderTablesIntoOne.py"), site_handler))
		#long arguments will happen, so no clustering
		#mergeSameHeaderTablesIntoOne.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mergeSameHeaderTablesIntoOne)
		workflow.mergeSameHeaderTablesIntoOne = mergeSameHeaderTablesIntoOne
	
	def initiateWorkflow(self, workflowName):
		"""
		2011-11-22
		"""
		# Create a abstract dag
		workflow = ADAG(workflowName)
		workflow.site_handler = self.site_handler
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		workflow.architecture = "x86_64"
		workflow.operatingSystem = "linux"
		workflow.namespace = "workflow"
		workflow.version="1.0"
		#clusters_size controls how many jobs will be aggregated as a single job.
		workflow.clusters_size = self.clusters_size
		return workflow
	
	@classmethod
	def addBAMIndexJob(cls, workflow, BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					inputBamF=None,\
					parentJob=None, namespace='workflow', version='1.0',\
					stageOutFinalOutput=True, javaMaxMemory=2500,\
					**keywords):
		"""
		2011-11-20
		"""
		index_sam_job = Job(namespace=namespace, name=BuildBamIndexFilesJava.name, version=version)
		baiFile = File('%s.bai'%inputBamF.name)
		index_sam_job.addArguments("-Xms128m", "-Xmx%sm"%(javaMaxMemory), "-jar", BuildBamIndexFilesJar, "VALIDATION_STRINGENCY=LENIENT", \
						"INPUT=", inputBamF, "OUTPUT=", baiFile)
		yh_pegasus.setJobProperRequirement(index_sam_job, job_max_memory=javaMaxMemory)
		index_sam_job.bamFile = inputBamF
		index_sam_job.baiFile = baiFile
		if stageOutFinalOutput:
			index_sam_job.uses(inputBamF, transfer=True, register=True, link=Link.INPUT)
			index_sam_job.uses(baiFile, transfer=True, register=True, link=Link.OUTPUT)
		else:
			pass	#don't register the files so leave them there
		workflow.addJob(index_sam_job)
		if parentJob:
			workflow.depends(parent=parentJob, child=index_sam_job)
		return index_sam_job
	
	def registerFilesAsInputToJob(self, job, inputFileList):
		"""
		2011-11-25
		"""
		for inputFile in inputFileList:
			job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
	
	def addRefFastaFaiIndexJob(self, workflow, samtools=None, refFastaF=None, ):
		"""
		2011-11-25
			the *.fai file of refFastaF is required for GATK
		"""
		refFastaIndexFname = '%s.fai'%(refFastaF.name)	# the .fai file is required for GATK
		refFastaIndexF = File(refFastaIndexFname)
		#not os.path.isfile(refFastaIndexFname)
		fastaIndexJob = Job(namespace=workflow.namespace, name=samtools.name, version=workflow.version)
		fastaIndexJob.addArguments("faidx", refFastaF)
		fastaIndexJob.uses(refFastaF,  transfer=True, register=False, link=Link.INPUT)
		fastaIndexJob.uses(refFastaIndexFname, transfer=True, register=False, link=Link.OUTPUT)
		fastaIndexJob.refFastaIndexF = refFastaIndexF
		yh_pegasus.setJobProperRequirement(fastaIndexJob, job_max_memory=1000)
		workflow.addJob(fastaIndexJob)
		return fastaIndexJob
	
	def addRefFastaDictJob(self, workflow, createSequenceDictionaryJava=None, refFastaF=None):
		"""
		2011-11-25
			# the .dict file is required for GATK
		"""
		refFastaDictFname = '%s.dict'%(os.path.splitext(refFastaF.name)[0])
		refFastaDictF = File(refFastaDictFname)
		#not os.path.isfile(refFastaDictFname) or 
		fastaDictJob = Job(namespace=workflow.namespace, name=createSequenceDictionaryJava.name, version=workflow.version)
		fastaDictJob.addArguments('-jar', createSequenceDictionaryJar, \
				'REFERENCE=', refFastaF, 'OUTPUT=', refFastaDictF)
		fastaDictJob.uses(refFastaF,  transfer=True, register=False, link=Link.INPUT)
		fastaDictJob.uses(refFastaDictF, transfer=True, register=False, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(fastaDictJob, job_max_memory=1000)
		fastaDictJob.refFastaDictF = refFastaDictF
		workflow.addJob(fastaDictJob)
		return fastaDictJob
	
	def addStatMergeJob(self, workflow, statMergeProgram=None, outputF=None, \
					parentJobLs=[], \
					namespace=None, version=None, extraDependentInputLs=[], transferOutput=True, extraArguments=None):
		"""
		2011-11-28
			moved from CalculateVCFStatPipeline.py
		2011-11-17
			add argument extraArguments
		"""
		statMergeJob = Job(namespace=namespace, name=statMergeProgram.name, version=version)
		statMergeJob.addArguments('-o', outputF)
		if extraArguments:
			statMergeJob.addArguments(extraArguments)
		statMergeJob.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		statMergeJob.output = outputF
		workflow.addJob(statMergeJob)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=statMergeJob)
		for input in extraDependentInputLs:
			statMergeJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		return statMergeJob
	
	def addInputToStatMergeJob(self, workflow, statMergeJob=None, inputF=None, \
							parentJobLs=[], \
							namespace=None, version=None, extraDependentInputLs=[]):
		"""
		2011-11-28
			moved from CalculateVCFStatPipeline.py
		"""
		statMergeJob.addArguments(inputF)
		statMergeJob.uses(inputF, transfer=False, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=statMergeJob)