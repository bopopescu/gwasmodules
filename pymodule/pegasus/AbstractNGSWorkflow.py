#!/usr/bin/env python
"""
2011-11-22
	a common class for pegasus workflows that work on NGS (next-gen sequencing) data
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, Genome
from Pegasus.DAX3 import *
from AbstractWorkflow import AbstractWorkflow

class AbstractNGSWorkflow(AbstractWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractWorkflow.option_default_dict.copy()
	option_default_dict.update(AbstractWorkflow.db_option_dict)
	
	option_default_dict.update({
						('ref_ind_seq_id', 1, int): [524, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/gatk/dist", '', 1, 'GATK folder containing its jar binaries'],\
						('tabixPath', 1, ): ["%s/bin/tabix", '', 1, 'path to the tabix binary', ],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("localDataDir", 0, ): ["", 'D', 1, 'localDataDir should contain same files as dataDir but accessible locally.\
									If not given, use the default stored in db. This argument is used to find all input files available.'],\
						
						('checkEmptyVCFByReading', 0, int):[0, 'E', 0, 'toggle to check if a vcf file is empty by reading its content'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractWorkflow.__init__(self, **keywords)
		#from pymodule import ProcessOptions
		#self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
		#												class_to_have_attr=self)
		self.samtools_path = self.insertHomePath(self.samtools_path, self.home_path)
		self.picard_path =  self.insertHomePath(self.picard_path, self.home_path)
		self.gatk_path =  self.insertHomePath(self.gatk_path, self.home_path)
		self.tabixPath =  self.insertHomePath(self.tabixPath, self.home_path)
		
		import re
		self.chr_pattern = re.compile(r'(\w+\d+).*')
		self.contig_id_pattern = re.compile(r'Contig(\d+).*')
	
	def registerJars(self, workflow=None, ):
		"""
		2011-11-22
			register jars to be used in the worflow
		"""
		AbstractWorkflow.registerJars(self)
		
		site_handler = self.site_handler
		#add the MergeSamFiles.jar file into workflow
		abs_path = os.path.join(self.picard_path, 'MergeSamFiles.jar')
		mergeSamFilesJar = File(abs_path)	#using abs_path avoids add this jar to every job as Link.INPUT
		mergeSamFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		self.addFile(mergeSamFilesJar)
		self.mergeSamFilesJar = mergeSamFilesJar
		
		abs_path = os.path.join(self.picard_path, 'BuildBamIndex.jar')
		BuildBamIndexFilesJar = File(abs_path)
		BuildBamIndexFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		self.addFile(BuildBamIndexFilesJar)
		self.BuildBamIndexFilesJar = BuildBamIndexFilesJar
		
		abs_path = os.path.join(self.gatk_path, 'GenomeAnalysisTK.jar')
		genomeAnalysisTKJar = File(abs_path)
		genomeAnalysisTKJar.addPFN(PFN("file://" + abs_path, site_handler))
		self.addFile(genomeAnalysisTKJar)
		self.genomeAnalysisTKJar = genomeAnalysisTKJar
		
		abs_path = os.path.join(self.picard_path, 'CreateSequenceDictionary.jar')
		createSequenceDictionaryJar = File(abs_path)
		createSequenceDictionaryJar.addPFN(PFN("file://" + abs_path, site_handler))
		self.addFile(createSequenceDictionaryJar)
		self.createSequenceDictionaryJar = createSequenceDictionaryJar
		
		abs_path = os.path.join(self.picard_path, 'AddOrReplaceReadGroupsAndCleanSQHeader.jar')
		addOrReplaceReadGroupsAndCleanSQHeaderJar = File(abs_path)
		addOrReplaceReadGroupsAndCleanSQHeaderJar.addPFN(PFN("file://" + abs_path, site_handler))
		self.addFile(addOrReplaceReadGroupsAndCleanSQHeaderJar)
		self.addOrReplaceReadGroupsAndCleanSQHeaderJar = addOrReplaceReadGroupsAndCleanSQHeaderJar
		
		abs_path = os.path.join(self.picard_path, 'AddOrReplaceReadGroups.jar')
		addOrReplaceReadGroupsJar = File(abs_path)
		addOrReplaceReadGroupsJar.addPFN(PFN("file://" + abs_path, site_handler))
		self.addFile(addOrReplaceReadGroupsJar)
		self.addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar
		
		abs_path = os.path.join(self.picard_path, 'MarkDuplicates.jar')
		MarkDuplicatesJar = File(abs_path)
		MarkDuplicatesJar.addPFN(PFN("file://" + abs_path, site_handler))
		self.addFile(MarkDuplicatesJar)
		self.MarkDuplicatesJar = MarkDuplicatesJar
		
		abs_path = os.path.join(self.picard_path, 'SplitReadFile.jar')
		SplitReadFileJar = File(abs_path)
		SplitReadFileJar.addPFN(PFN("file://" + abs_path, site_handler))
		self.addFile(SplitReadFileJar)
		self.SplitReadFileJar = SplitReadFileJar
	
	def registerCustomJars(self, workflow=None, ):
		"""
		2012.1.9
		"""
		pass
	
	def registerExecutables(self, workflow=None):
		"""
		2012.1.9 a symlink to registerCommonExecutables()
		"""
		AbstractWorkflow.registerExecutables(self)
		self.registerCommonExecutables(workflow=workflow)
	
	def registerCommonExecutables(self, workflow=None):
		"""
		2012.7.25
			add noDefaultClustersSizeExecutableList
		2011-11-22
		"""
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableList = []
		noDefaultClustersSizeExecutableList = []
		
		
		selectAndSplit = Executable(namespace=namespace, name="SelectAndSplitAlignment", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		selectAndSplit.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "SelectAndSplitAlignment.py"), site_handler))
		self.addExecutable(selectAndSplit)
		self.selectAndSplit = selectAndSplit
		
		selectAndSplitFasta = Executable(namespace=namespace, name="SelectAndSplitFastaRecords", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		selectAndSplitFasta.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "SelectAndSplitFastaRecords.py"), site_handler))
		self.addExecutable(selectAndSplitFasta)
		self.selectAndSplitFasta = selectAndSplitFasta
		
		samtools = Executable(namespace=namespace, name="samtools", version=version, os=operatingSystem, arch=architecture, installed=True)
		samtools.addPFN(PFN("file://" + self.samtools_path, site_handler))
		self.addExecutable(samtools)
		self.samtools = samtools
		
		java = Executable(namespace=namespace, name="java", version=version, os=operatingSystem, arch=architecture, installed=True)
		java.addPFN(PFN("file://" + self.javaPath, site_handler))
		self.addExecutable(java)
		self.java = java
		
		addOrReplaceReadGroupsJava = Executable(namespace=namespace, name="addOrReplaceReadGroupsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		addOrReplaceReadGroupsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		self.addExecutable(addOrReplaceReadGroupsJava)
		self.addOrReplaceReadGroupsJava = addOrReplaceReadGroupsJava
		
		genotyperJava = Executable(namespace=namespace, name="genotyperJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		genotyperJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#clustering is controlled by a separate parameter
		#genotyperJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(genotyperJava)
		self.genotyperJava = genotyperJava
		
		BuildBamIndexFilesJava = Executable(namespace=namespace, name="BuildBamIndexFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		BuildBamIndexFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		self.addExecutable(BuildBamIndexFilesJava)
		self.BuildBamIndexFilesJava = BuildBamIndexFilesJava
		
		createSequenceDictionaryJava = Executable(namespace=namespace, name="createSequenceDictionaryJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		createSequenceDictionaryJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		self.addExecutable(createSequenceDictionaryJava)
		self.createSequenceDictionaryJava = createSequenceDictionaryJava
		
		DOCWalkerJava = Executable(namespace=namespace, name="DepthOfCoverageWalker", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		#no clusters_size for this because it could run on a whole bam for hours
		#DOCWalkerJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		DOCWalkerJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		self.addExecutable(DOCWalkerJava)
		self.DOCWalkerJava = DOCWalkerJava
		
		VariousReadCountJava = Executable(namespace=namespace, name="VariousReadCountJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		VariousReadCountJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#no clusters_size for this because it could run on a whole bam for hours
		#VariousReadCountJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(VariousReadCountJava)
		self.VariousReadCountJava = VariousReadCountJava
		
		
		MarkDuplicatesJava = Executable(namespace=namespace, name="MarkDuplicatesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		MarkDuplicatesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#MarkDuplicatesJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(MarkDuplicatesJava)
		self.MarkDuplicatesJava = MarkDuplicatesJava
		
		SelectVariantsJava = Executable(namespace=namespace, name="SelectVariantsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		SelectVariantsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableList.append(SelectVariantsJava)
		
		CombineVariantsJava = Executable(namespace=namespace, name="CombineVariantsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		CombineVariantsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		CombineVariantsJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(CombineVariantsJava)
		self.CombineVariantsJava = CombineVariantsJava
		
		CallVariantBySamtools = Executable(namespace=namespace, name="CallVariantBySamtools", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		CallVariantBySamtools.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/CallVariantBySamtools.sh"), site_handler))
		#clustering is controlled by a separate parameter
		#CallVariantBySamtools.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(CallVariantBySamtools)
		self.CallVariantBySamtools = CallVariantBySamtools
		
		genotypeCallByCoverage = Executable(namespace=namespace, name="GenotypeCallByCoverage", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		genotypeCallByCoverage.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "GenotypeCallByCoverage.py"), site_handler))
		genotypeCallByCoverage.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(genotypeCallByCoverage)
		self.genotypeCallByCoverage = genotypeCallByCoverage
		
		mergeGenotypeMatrix = Executable(namespace=namespace, name="MergeGenotypeMatrix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		mergeGenotypeMatrix.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "MergeGenotypeMatrix.py"), site_handler))
		mergeGenotypeMatrix.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(mergeGenotypeMatrix)
		self.mergeGenotypeMatrix = mergeGenotypeMatrix
		
		bgzip_tabix = Executable(namespace=namespace, name="bgzip_tabix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		bgzip_tabix.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/bgzip_tabix.sh"), site_handler))
		bgzip_tabix.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(bgzip_tabix)
		self.bgzip_tabix = bgzip_tabix
		
		vcf_convert = Executable(namespace=namespace, name="vcf_convert", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		vcf_convert.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_convert.sh"), site_handler))
		vcf_convert.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(vcf_convert)
		self.vcf_convert = vcf_convert
		
		vcf_isec = Executable(namespace=namespace, name="vcf_isec", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_isec.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_isec.sh"), site_handler))
		vcf_isec.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(vcf_isec)
		self.vcf_isec = vcf_isec
		
		vcf_concat = Executable(namespace=namespace, name="vcf_concat", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_concat.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		vcf_concat.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(vcf_concat)
		self.vcf_concat = vcf_concat
		
		vcfSubsetPath = os.path.join(self.home_path, "bin/vcftools/vcf-subset")
		self.vcfSubsetPath = vcfSubsetPath	#vcfSubsetPath is first argument to vcfSubsetPath
		vcfSubset = Executable(namespace=namespace, name="vcfSubset", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcfSubset.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcfSubset.sh"), site_handler))
		vcfSubset.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(vcfSubset)
		self.vcfSubset = vcfSubset
		
		concatGATK = Executable(namespace=namespace, name="concatGATK", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		concatGATK.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		concatGATK.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(concatGATK)
		self.concatGATK = concatGATK
		
		concatSamtools = Executable(namespace=namespace, name="concatSamtools", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		concatSamtools.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		concatSamtools.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(concatSamtools)
		self.concatSamtools = concatSamtools
		
		calcula = Executable(namespace=namespace, name="CalculatePairwiseDistanceOutOfSNPXStrainMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		calcula.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py"), site_handler))
		calcula.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(calcula)
		self.calcula = calcula
		
		
		
		ReduceMatrixByChosenColumn = Executable(namespace=namespace, name="ReduceMatrixByChosenColumn", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		ReduceMatrixByChosenColumn.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/ReduceMatrixByChosenColumn.py"), site_handler))
		ReduceMatrixByChosenColumn.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(ReduceMatrixByChosenColumn)
		self.ReduceMatrixByChosenColumn = ReduceMatrixByChosenColumn
		
		ReduceMatrixBySumSameKeyColsAndThenDivide = Executable(namespace=namespace, name="ReduceMatrixBySumSameKeyColsAndThenDivide", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		ReduceMatrixBySumSameKeyColsAndThenDivide.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/ReduceMatrixBySumSameKeyColsAndThenDivide.py"), \
															site_handler))
		ReduceMatrixBySumSameKeyColsAndThenDivide.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(ReduceMatrixBySumSameKeyColsAndThenDivide)
		self.ReduceMatrixBySumSameKeyColsAndThenDivide = ReduceMatrixBySumSameKeyColsAndThenDivide
		
		ReduceMatrixByAverageColumnsWithSameKey = Executable(namespace=namespace, name="ReduceMatrixByAverageColumnsWithSameKey", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		ReduceMatrixByAverageColumnsWithSameKey.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/ReduceMatrixByAverageColumnsWithSameKey.py"), site_handler))
		ReduceMatrixByAverageColumnsWithSameKey.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(ReduceMatrixByAverageColumnsWithSameKey)
		self.ReduceMatrixByAverageColumnsWithSameKey = ReduceMatrixByAverageColumnsWithSameKey
		
		vcftoolsPath = os.path.join(self.home_path, "bin/vcftools/vcftools")
		self.vcftoolsPath = vcftoolsPath	#vcftoolsPath is first argument to vcftoolsWrapper
		vcftoolsWrapper = Executable(namespace=namespace, name="vcftoolsWrapper", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcftoolsWrapper.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcftoolsWrapper.sh"), site_handler))
		vcftoolsWrapper.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		vcftoolsWrapper.vcftoolsPath = vcftoolsPath
		self.addExecutable(vcftoolsWrapper)
		self.vcftoolsWrapper = vcftoolsWrapper
		
		tabix = Executable(namespace=namespace, name="tabix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		tabix.addPFN(PFN("file://" + self.tabixPath, site_handler))
		tabix.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(tabix)
		self.tabix = tabix
		
		#2011.12.21 moved from FilterVCFPipeline.py
		FilterVCFByDepthJava = Executable(namespace=namespace, name="FilterVCFByDepth", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		FilterVCFByDepthJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		FilterVCFByDepthJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(FilterVCFByDepthJava)
		self.FilterVCFByDepthJava = FilterVCFByDepthJava
		
		#2011.12.21	for OutputVCFSiteStat.py
		tabixRetrieve = Executable(namespace=namespace, name="tabixRetrieve", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		tabixRetrieve.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/tabixRetrieve.sh"), site_handler))
		tabixRetrieve.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(tabixRetrieve)
		self.tabixRetrieve = tabixRetrieve
		
		#2012.3.1
		MergeFiles = Executable(namespace=namespace, name="MergeFiles", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		MergeFiles.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "shell/MergeFiles.sh"), site_handler))
		executableList.append(MergeFiles)
		
		#2012.7.25
		MergeVCFReplicateHaplotypesJava= Executable(namespace=namespace, name="MergeVCFReplicateHaplotypesJava", \
											version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		MergeVCFReplicateHaplotypesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		noDefaultClustersSizeExecutableList.append(MergeVCFReplicateHaplotypesJava)
		
		for executable in executableList:
			executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			self.addExecutable(executable)
			setattr(self, executable.name, executable)
		
		for executable in noDefaultClustersSizeExecutableList:
			self.addExecutable(executable)
			setattr(self, executable.name, executable)
	
	def addBAMIndexJob(self, workflow=None, BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					inputBamF=None,\
					parentJobLs=[], namespace='workflow', version='1.0',\
					stageOutFinalOutput=True, javaMaxMemory=2500,\
					**keywords):
		"""
		 2012.4.12
			remove argument parentJob and stop adding it to parentJobLs, which causes an insidious bug
				that accumulates parent jobs from multiple calls of addBAMIndexJob() into parentJobLs
				(they all become parents of this bam index job.)
		2012.3.22
			bugfix, change argument parentJobLs's default value from [] to None. [] would make every run have the same parentJobLs 
			proper transfer/register setup
		2011-11-20
		"""
		index_sam_job = Job(namespace=getattr(self, 'namespace', namespace), name=BuildBamIndexFilesJava.name, \
						version=getattr(self, 'version', version))
		baiFile = File('%s.bai'%inputBamF.name)
		index_sam_job.addArguments("-Xms128m", "-Xmx%sm"%(javaMaxMemory), "-jar", BuildBamIndexFilesJar, "VALIDATION_STRINGENCY=LENIENT", \
						"INPUT=", inputBamF, "OUTPUT=", baiFile)
		index_sam_job.bamFile = inputBamF
		index_sam_job.baiFile = baiFile
		index_sam_job.output = baiFile	#2012.3.20
		index_sam_job.uses(inputBamF, transfer=True, register=True, link=Link.INPUT)
		#index_sam_job.uses(inputBamF, transfer=True, register=True, link=Link.OUTPUT)
		index_sam_job.uses(baiFile, transfer=stageOutFinalOutput, register=True, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(index_sam_job, job_max_memory=javaMaxMemory)
		self.addJob(index_sam_job)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=index_sam_job)
		return index_sam_job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.1.9
			abstract function
		"""
		pass
	
	
	def addRefFastaFaiIndexJob(self, workflow=None, samtools=None, refFastaF=None, ):
		"""
		2011-11-25
			the *.fai file of refFastaF is required for GATK
		"""
		refFastaIndexFname = '%s.fai'%(refFastaF.name)	# the .fai file is required for GATK
		refFastaIndexF = File(refFastaIndexFname)
		#not os.path.isfile(refFastaIndexFname)
		fastaIndexJob = Job(namespace=self.namespace, name=samtools.name, version=self.version)
		fastaIndexJob.addArguments("faidx", refFastaF)
		fastaIndexJob.uses(refFastaF,  transfer=True, register=True, link=Link.INPUT)
		fastaIndexJob.uses(refFastaIndexFname, transfer=True, register=True, link=Link.OUTPUT)
		fastaIndexJob.refFastaIndexF = refFastaIndexF
		yh_pegasus.setJobProperRequirement(fastaIndexJob, job_max_memory=1000)
		self.addJob(fastaIndexJob)
		return fastaIndexJob
	
	def addRefFastaDictJob(self, workflow=None, createSequenceDictionaryJava=None, refFastaF=None):
		"""
		2011-11-25
			# the .dict file is required for GATK
		"""
		refFastaDictFname = '%s.dict'%(os.path.splitext(refFastaF.name)[0])
		refFastaDictF = File(refFastaDictFname)
		#not os.path.isfile(refFastaDictFname) or 
		fastaDictJob = Job(namespace=self.namespace, name=createSequenceDictionaryJava.name, version=self.version)
		fastaDictJob.addArguments('-jar', createSequenceDictionaryJar, \
				'REFERENCE=', refFastaF, 'OUTPUT=', refFastaDictF)
		fastaDictJob.uses(refFastaF,  transfer=True, register=True, link=Link.INPUT)
		fastaDictJob.uses(refFastaDictF, transfer=True, register=True, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(fastaDictJob, job_max_memory=1000)
		fastaDictJob.refFastaDictF = refFastaDictF
		self.addJob(fastaDictJob)
		return fastaDictJob
	
	def addRefFastaJobDependency(self, workflow=None, job=None, refFastaF=None, fastaDictJob=None, refFastaDictF=None, fastaIndexJob = None, refFastaIndexF = None):
		"""
		2011-9-14
		"""
		if fastaIndexJob:	#2011-7-22 if job doesn't exist, don't add it. means this job isn't necessary to run.
			self.depends(parent=fastaIndexJob, child=job)
			job.uses(refFastaIndexF, transfer=True, register=True, link=Link.INPUT)
		if fastaDictJob:
			self.depends(parent=fastaDictJob, child=job)
			job.uses(refFastaDictF, transfer=True, register=True, link=Link.INPUT)
		if fastaIndexJob or fastaDictJob:
			job.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
	
	def addVCFFormatConvertJob(self, workflow=None, vcf_convert=None, parentJob=None, inputF=None, outputF=None, \
							namespace=None, version=None, transferOutput=False):
		"""
		2011-11-4
		"""
		vcf_convert_job = Job(namespace=getattr(self, 'namespace', namespace), name=vcf_convert.name, \
							version=getattr(self, 'version', version))
		vcf_convert_job.addArguments(inputF, outputF)
		vcf_convert_job.uses(inputF, transfer=False, register=True, link=Link.INPUT)
		vcf_convert_job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		self.addJob(vcf_convert_job)
		self.depends(parent=parentJob, child=vcf_convert_job)
		return vcf_convert_job
	
	def addBGZIP_tabix_Job(self, workflow=None, bgzip_tabix=None, parentJob=None, inputF=None, outputF=None, \
							namespace=None, version=None, transferOutput=False, parentJobLs=[], tabixArguments=""):
		"""
		2011.12.20
			pass additional tabix arguments to bgzip_tabix shell script
		2011-11-4
		
		"""
		bgzip_tabix_job = Job(namespace=getattr(self, 'namespace', namespace), name=bgzip_tabix.name, \
							version=getattr(self, 'version', version))
		tbi_F = File("%s.tbi"%outputF.name)
		bgzip_tabix_job.addArguments(inputF, outputF)
		if tabixArguments:
			bgzip_tabix_job.addArguments(tabixArguments)
		bgzip_tabix_job.uses(inputF, transfer=False, register=True, link=Link.INPUT)
		bgzip_tabix_job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		bgzip_tabix_job.uses(tbi_F, transfer=transferOutput, register=True, link=Link.OUTPUT)
		bgzip_tabix_job.output = outputF
		bgzip_tabix_job.tbi_F = tbi_F
		self.addJob(bgzip_tabix_job)
		if parentJob:
			self.depends(parent=parentJob, child=bgzip_tabix_job)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=bgzip_tabix_job)
		return bgzip_tabix_job
	
	def addVCFConcatJob(self, workflow=None, concatExecutable=None, parentDirJob=None, outputF=None, \
							namespace=None, version=None, transferOutput=True, vcf_job_max_memory=500):
		"""
		2011-11-5
		"""
		#2011-9-22 union of all samtools intervals for one contig
		vcfConcatJob = Job(namespace=getattr(self, 'namespace', namespace), name=concatExecutable.name, \
						version=getattr(self, 'version', version))
		vcfConcatJob.addArguments(outputF)
		vcfConcatJob.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		tbi_F = File("%s.tbi"%outputF.name)
		vcfConcatJob.uses(tbi_F, transfer=transferOutput, register=True, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(vcfConcatJob, job_max_memory=vcf_job_max_memory)
		self.addJob(vcfConcatJob)
		self.depends(parent=parentDirJob, child=vcfConcatJob)
		return vcfConcatJob
	
	def addVCFSubsetJob(self, workflow=None, executable=None, vcfSubsetPath=None, sampleIDFile=None,\
					inputVCF=None, outputF=None, \
					parentJobLs=[], namespace=None, version=None, transferOutput=True, job_max_memory=200,\
					extraArguments=None, extraDependentInputLs=[]):
		"""
		2012.5.10
		"""
		#2011-9-22 union of all samtools intervals for one contig
		job = Job(namespace=getattr(self, 'namespace', namespace), name=executable.name, \
						version=getattr(self, 'version', version))
		job.addArguments(vcfSubsetPath, sampleIDFile, inputVCF, outputF)
		
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputVCF, transfer=True, register=True, link=Link.INPUT)
		job.uses(sampleIDFile, transfer=True, register=True, link=Link.INPUT)
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

	def addVCF2MatrixJob(self, workflow=None, executable=None, inputVCF=None, outputFile=None, \
						refFastaF=None, run_type=3, numberOfReadGroups=10, seqCoverageF=None, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.5.8
			executable is GenotypeCallByCoverage
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		#2012.5.8 "-n 10" means numberOfReadGroups is 10 but it's irrelevant when "-y 3" (run_type =3, read from vcf without filter)
		job.addArguments("-i", inputVCF, "-n 10", "-o", outputFile, '-y 3')
		if refFastaF:
			job.addArguments("-e", refFastaF)
			job.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
		if extraArguments:
			job.addArguments(extraArguments)
		if seqCoverageF:
			job.addArguments("-q", seqCoverageF)
			job.uses(seqCoverageF, transfer=True, register=True, link=Link.INPUT)
		job.uses(inputVCF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputFile
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			if input:
				job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addCalculatePairwiseDistanceFromSNPXStrainMatrixJob(self, workflow=None, executable=None, inputFile=None, outputFile=None, \
						min_MAF=0, max_NA_rate=0.4, convertHetero2NA=0, hetHalfMatchDistance=0.5,\
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.5.11
			executable is CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py
			
			#add the pairwise distance matrix job after filter is done
			calcula_job = Job(namespace=namespace, name=calcula.name, version=version)
			
			calcula_job.addArguments("-i", genotypeCallOutput, "-n", str(self.min_MAF), \
						"-o", calculaOutput, '-m', repr(self.max_NA_rate), '-c', str(self.convertHetero2NA),\
						'-H', repr(self.hetHalfMatchDistance))
			calcula_job.uses(genotypeCallOutput, transfer=False, register=False, link=Link.INPUT)
			calcula_job.uses(calculaOutput, transfer=True, register=False, link=Link.OUTPUT)
			
			workflow.addJob(calcula_job)
			workflow.depends(parent=genotypeCallByCoverage_job, child=calcula_job)
			workflow.depends(parent=matrixDirJob, child=calcula_job)
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		job.addArguments("-i", inputFile,  "-n %s"%(min_MAF), \
					"-o", outputFile, '-m %s'%(max_NA_rate), '-c %s'%(convertHetero2NA),\
					'-H %s'%(hetHalfMatchDistance))
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputFile
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			if input:
				job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	@classmethod
	def findProperVCFDirIdentifier(cls, vcfDir, defaultName='vcf1'):
		"""
		2011.11.28
		"""
		#name to distinguish between vcf1Dir, and vcf2Dir
		folderNameLs = vcfDir.split('/')
		vcf1Name=None
		for i in range(1, len(folderNameLs)+1):
			if folderNameLs[-i]:	#start from the end to find the non-empty folder name
				#the trailing "/" on self.vcf1Dir could mean  an empty string
				vcf1Name = folderNameLs[-i]
				break
		if not vcf1Name:
			vcf1Name = defaultName
		return vcf1Name
	
	def addSelectVariantsJob(self, workflow=None, SelectVariantsJava=None, genomeAnalysisTKJar=None, inputF=None, outputF=None, \
					refFastaFList=[], parentJobLs=None, \
					extraDependentInputLs=[], transferOutput=True, extraArguments=None, job_max_memory=2000, interval=None,\
					**keywords):
		"""
		2011-12.5
		"""
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		refFastaF = refFastaFList[0]
		job = Job(namespace=self.namespace, name=SelectVariantsJava.name, version=self.version)
		job.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T SelectVariants", "-R", refFastaF, \
						"--variant", inputF, "-o ", outputF, "-L %s"%(interval))
		if extraArguments:
			job.addArguments(extraArguments)
		for refFastaFile in refFastaFList:
			job.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		self.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def getVCFFileID2path(self, inputDir):
		"""
		2011-12-1
		"""
		sys.stderr.write("Getting all vcf files from %s "%(inputDir))
		vcfFileID2path = {}
		for inputFname in os.listdir(inputDir):
			inputAbsPath = os.path.join(os.path.abspath(inputDir), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False) and not NextGenSeq.isVCFFileEmpty(inputAbsPath):
				fileID = self.chr_pattern.search(inputFname).group(1)
				if fileID in vcfFileID2path:
					sys.stderr.write("fileID %s already has value (%s) in dictionary. but now a 2nd file %s overwrites previous value.\n"%\
									(fileID, vcfFileID2path.get(fileID), inputFname))
				vcfFileID2path[fileID] = inputAbsPath
		sys.stderr.write("  found %s files.\n"%(len(vcfFileID2path)))
		return vcfFileID2path
	
	def addCheckTwoVCFOverlapJob(self, workflow=None, executable=None, vcf1=None, vcf2=None, chromosome=None, chrLength=None, \
					outputFnamePrefix=None, parentJobLs=[], \
					extraDependentInputLs=[], transferOutput=False, extraArguments=None, job_max_memory=1000, \
					**keywords):
		"""
		2011.12.9
		"""
		overlapSitePosF = File('%s_overlapSitePos.tsv'%(outputFnamePrefix))
		outputF = File('%s_overlap.tsv'%(outputFnamePrefix))
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		job.addArguments("-i", vcf1, "-j", vcf2, "-c", chromosome, "-l %s"%(chrLength), '-o', outputF,\
						'-O', outputFnamePrefix)
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(vcf1, transfer=True, register=True, link=Link.INPUT)
		job.uses(vcf2, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.uses(overlapSitePosF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.overlapSitePosF = overlapSitePosF
		job.output = outputF
		self.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addFilterVCFByDepthJob(self, workflow=None, FilterVCFByDepthJava=None, genomeAnalysisTKJar=None, \
							refFastaFList=None, inputVCFF=None, outputVCFF=None, outputSiteStatF=None, \
							parentJobLs=[], alnStatForFilterF=None, \
							job_max_memory=1000, extraDependentInputLs=[], onlyKeepBiAllelicSNP=False, \
							namespace=None, version=None, transferOutput=False, **keywords):
		"""
		2011.12.20
			moved from FilterVCFPipeline.py
			add argument transferOutput, outputSiteStatF
		"""
		# Add a mkdir job for any directory.
		filterByDepthJob = Job(namespace=getattr(self, 'namespace', namespace), name=FilterVCFByDepthJava.name, \
							version=getattr(self, 'version', version))
		refFastaF = refFastaFList[0]
		filterByDepthJob.addArguments("-Xmx%sm"%(job_max_memory), "-jar", genomeAnalysisTKJar, "-R", refFastaF, "-T FilterVCFByDepth", \
							"--variant", inputVCFF, "-depthFname", alnStatForFilterF)
		if outputVCFF:
			filterByDepthJob.addArguments("-o", outputVCFF)
			filterByDepthJob.uses(outputVCFF, transfer=transferOutput, register=True, link=Link.OUTPUT)
			filterByDepthJob.output = outputVCFF
		
		if outputSiteStatF:
			filterByDepthJob.addArguments("-ssFname", outputSiteStatF)
			filterByDepthJob.uses(outputSiteStatF, transfer=transferOutput, register=True, link=Link.OUTPUT)
			filterByDepthJob.outputSiteStatF = outputSiteStatF
		
		if onlyKeepBiAllelicSNP:
			filterByDepthJob.addArguments("--onlyKeepBiAllelicSNP")
		
		for refFastaFile in refFastaFList:
			filterByDepthJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		filterByDepthJob.uses(alnStatForFilterF, transfer=True, register=True, link=Link.INPUT)
		filterByDepthJob.uses(inputVCFF, transfer=True, register=True, link=Link.INPUT)		
		for input in extraDependentInputLs:
			filterByDepthJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		self.addJob(filterByDepthJob)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=filterByDepthJob)
		yh_pegasus.setJobProperRequirement(filterByDepthJob, job_max_memory=job_max_memory)
		return filterByDepthJob
	
	def addFilterJobByvcftools(self, workflow=None, vcftoolsWrapper=None, inputVCFF=None, outputFnamePrefix=None, \
							parentJobLs=None, snpMisMatchStatFile=None, minMAC=2, minMAF=None, maxSNPMissingRate=0.9,\
							namespace=None, version=None, extraDependentInputLs=[], transferOutput=False, \
							outputFormat='--recode', extraArguments=None):
		"""
		2012.5.9
			moved from FilterVCFPipeline.py
			add argument outputFormat to make "--recode" explicit and allow other output formats
			add argument extraArguments
			
			this could be just used to output vcf in various formats
		2011-11-21
			argument vcftools is replaced with a wrapper, which takes vcftools path as 1st argument
		"""
		# Add a mkdir job for any directory.
		vcftoolsJob = Job(namespace=getattr(self, 'namespace', namespace), name=vcftoolsWrapper.name, \
						version=getattr(self, 'version', version))
		vcftoolsJob.addArguments(self.vcftoolsPath)	#2011-11-21
		if inputVCFF.name[-2:]=='gz':
			vcftoolsJob.addArguments("--gzvcf", inputVCFF)
		else:
			vcftoolsJob.addArguments("--vcf", inputVCFF)
		vcftoolsJob.addArguments("--out", outputFnamePrefix, outputFormat)
		if snpMisMatchStatFile:
			vcftoolsJob.addArguments("--positions", snpMisMatchStatFile)
			vcftoolsJob.uses(snpMisMatchStatFile, transfer=True, register=True, link=Link.INPUT)
		
		if maxSNPMissingRate is not None:
			vcftoolsJob.addArguments("--geno %s"%(1-maxSNPMissingRate))
		if minMAF is not None:
			vcftoolsJob.addArguments("--maf %s"%(minMAF))
		if minMAC is not None:
			vcftoolsJob.addArguments("--mac %s"%(minMAC))
		if extraArguments:
			vcftoolsJob.addArguments(extraArguments)
		vcftoolsJob.uses(inputVCFF, transfer=True, register=True, link=Link.INPUT)
		for input in extraDependentInputLs:
			if input:
				vcftoolsJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		if outputFormat=='--recode':	#2012.5.9
			outputVCFF = File("%s.recode.vcf"%(outputFnamePrefix))
			vcftoolsJob.uses(outputVCFF, transfer=transferOutput, register=True, link=Link.OUTPUT)
			vcftoolsJob.output = outputVCFF
		elif outputFormat=='--plink-tped':
			output1 = File("%s.tped"%(outputFnamePrefix))
			output2 = File("%s.tfam"%(outputFnamePrefix))
			vcftoolsJob.output = output1
			vcftoolsJob.outputList=[output1, output2]
			for output in vcftoolsJob.outputList:
				vcftoolsJob.uses(output, transfer=transferOutput, register=True, link=Link.OUTPUT)
		elif outputFormat=='--plink':
			output1 = File("%s.ped"%(outputFnamePrefix))
			output2 = File("%s.fam"%(outputFnamePrefix))
			vcftoolsJob.output = output1
			vcftoolsJob.outputList=[output1, output2]
			for output in vcftoolsJob.outputList:
				vcftoolsJob.uses(output, transfer=transferOutput, register=True, link=Link.OUTPUT)
		elif outputFormat=='--IMPUTE':
			output1 = File("%s.impute.hap"%(outputFnamePrefix))
			output2 = File("%s.impute.hap.legend"%(outputFnamePrefix))
			output3 = File("%s.impute.hap.indv"%(outputFnamePrefix))
			vcftoolsJob.output = output1
			vcftoolsJob.outputList=[output1, output2, output3]
			for output in vcftoolsJob.outputList:
				vcftoolsJob.uses(output, transfer=transferOutput, register=True, link=Link.OUTPUT)
		elif outputFormat=='--BEAGLE-GL':
			output1 = File("%s.BEAGLE.GL"%(outputFnamePrefix))
			vcftoolsJob.output = output1
			vcftoolsJob.outputList=[output1]
			for output in vcftoolsJob.outputList:
				vcftoolsJob.uses(output, transfer=transferOutput, register=True, link=Link.OUTPUT)
		else:
			vcftoolsJob.output = None
			vcftoolsJob.outputList = []
		self.addJob(vcftoolsJob)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=vcftoolsJob)
		return vcftoolsJob
	
	def addCalculateTwoVCFSNPMismatchRateJob(self, workflow=None, executable=None, \
							vcf1=None, vcf2=None, snpMisMatchStatFile=None, \
							maxSNPMismatchRate=1.0, parentJobLs=[], \
							job_max_memory=1000, extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2011.12.20
			
		"""
		job = Job(namespace=self.namespace, name=executable.name, \
												version=self.version)
		job.addArguments("-i", vcf1, "-j", vcf2, \
						"-m %s"%(maxSNPMismatchRate), '-o', snpMisMatchStatFile)
		job.uses(vcf1, transfer=True, register=True, link=Link.INPUT)
		job.uses(vcf2, transfer=True, register=True, link=Link.INPUT)
		job.uses(snpMisMatchStatFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=job)
		return job

	def addTabixRetrieveJob(self, workflow=None, executable=None, tabixPath=None, \
							inputF=None, outputF=None, regionOfInterest=None, includeHeader=True,\
							parentJobLs=[], job_max_memory=100, extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2011.12.20
			The executable should be tabixRetrieve (a tabix shell wrapper).
			
			http://samtools.sourceforge.net/tabix.shtml
			run something like below to extract data from regionOfInterest out of bgzipped&tabix-indexed file.
				tabix sorted.gff.gz chr1:10,000,000-20,000,000
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		job.addArguments(tabixPath)
		job.addArguments(inputF, outputF, regionOfInterest)
		if includeHeader:
			job.addArguments("-h")
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		
		job.output = outputF
		
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=job)
		return job
	
	
	def getContigIDFromFname(self, filename):
		"""
		2012.7.14 copied to pymodule.Genome
		2011-10-20
			
			If filename is like .../Contig0.filter_by_vcftools.recode.vcf.gz,
				It returns "0", excluding the "Contig".
				If you want "Contig" included, use getChrIDFromFname().
			If search fails, it returns the prefix in the basename of filename.
		"""
		return Genome.getContigIDFromFname(filename)
	
	def getChrFromFname(self, filename):
		"""
		2012.7.14 copied to pymodule.Genome
		2011-10-20
			filename example: Contig0.filter_by_vcftools.recode.vcf.gz
				It returns "Contig0".
				If you want just "0", use getContigIDFromFname().
			If search fails, it returns the prefix in the basename of filename.
		"""
		return Genome.getChrFromFname(filename)
		
	def addPutStuffIntoDBJob(self, workflow=None, executable=None, inputFileLs=[], \
					logFile=None, commit=False, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, sshDBTunnel=0, **keywords):
		"""
		2012.5.8 add sshDBTunnel
		2012.4.3
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd)
		if extraArguments:
			job.addArguments(extraArguments)
		if commit:
			job.addArguments("-c")
		if logFile:
			job.addArguments("--logFilename", logFile)
			job.uses(logFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			job.output = logFile
		for inputFile in inputFileLs:	#2012.4.3 this inputFile addition has to be at last. as it doesn't have options ahead of them. 
			job.addArguments(inputFile)
			job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		self.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addSamtoolsFlagstatJob(self, workflow=None, executable=None, inputF=None, outputF=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.4.3
			samtools (sam_stat.c) has been modified so that it could take one more optional argument to store the
				stats that are usually directed to stdout. 
			inputF is bam file. outputF is to store the output.
			
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		job.addArguments('flagstat', inputF, outputF)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addMergeVCFReplicateGenotypeColumnsJob(self, workflow=None, executable=None, genomeAnalysisTKJar=None, \
						inputF=None, outputF=None, replicateIndividualTag=None, \
						debugHaplotypeDistanceFile=None, \
						debugMajoritySupportFile=None,\
						refFastaFList=[], parentJobLs=[], extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.7.25
			use self.addGenericJob() and moved from AlignmentToTrioCallPipeline.py
		2012.6.1
			change MergeVCFReplicateGenotypeColumns to MergeVCFReplicateHaplotypes
			
		2012.4.2
			java -jar /home/crocea/script/gatk/dist/GenomeAnalysisTK.jar -T MergeVCFReplicateGenotypeColumns 
				-R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
				--variant /tmp/Contig0.vcf -o /tmp/contig0_afterMerge.vcf --onlyKeepBiAllelicSNP --replicateIndividualTag copy
		"""
		#GATK job
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		refFastaF = refFastaFList[0]
		extraArgumentList = [javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "MergeVCFReplicateHaplotypes",\
			"-R", refFastaF, "--variant:VCF", inputF, "--out", outputF,\
			'--onlyKeepBiAllelicSNP', "--replicateIndividualTag %s"%(replicateIndividualTag)]
		if debugHaplotypeDistanceFile:
			extraArgumentList.extend(["--debugHaplotypeDistanceFname", debugHaplotypeDistanceFile])
		if debugMajoritySupportFile:
			extraArgumentList.extend(["--debugMajoritySupportFname", debugMajoritySupportFile])
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			_extraDependentInputLs=[inputF] + refFastaFList
		else:
			import copy
			_extraDependentInputLs = copy.deepcopy(extraDependentInputLs)
			_extraDependentInputLs.append(inputF)
			_extraDependentInputLs.extend(refFastaFList)
		
		# don't pass inputF and outputF to addGenericJob() because it'll add "-i" and "-o" in front of the two respectively
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=_extraDependentInputLs, \
						extraOutputLs=[outputF, debugHaplotypeDistanceFile, debugMajoritySupportFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory)
		
		return job