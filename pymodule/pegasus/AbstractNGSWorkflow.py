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
						('maxContigID', 0, int): [None, 'x', 1, 'if contig/chromosome(non-sex) ID > this number, it will not be included. If None or 0, no restriction.', ],\
						('minContigID', 0, int): [None, 'V', 1, 'if contig/chromosome(non-sex) ID < this number, it will not be included. If None or 0, no restriction.', ],\
						("contigMaxRankBySize", 1, int): [1000, 'N', 1, 'maximum rank number (rank 1=biggest) controls how small a contig to be included in calling'],\
						("contigMinRankBySize", 1, int): [1, 'M', 1, 'minimum rank number (rank 1=biggest contig) controls how big a contig to be included in calling'],\
						
						('checkEmptyVCFByReading', 0, int):[0, 'E', 0, 'toggle to check if a vcf file is empty by reading its content'],\
						('excludeContaminant', 0, int):[0, '', 0, 'toggle this to exclude alignments from contaminated individuals, \n\
		(Individual.is_contaminated=1)'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	
	"""
	2012.9.18
		to silence this kind of error messages:
		
		##### ERROR
		##### ERROR MESSAGE:
				SAM/BAM file SAMFileReader{/u/home/eeskin2/polyacti/NetworkData/scratch/HaplotypeScore/HaplotypeScore_ISQ633_638.2012.Sep.18T110829/individual_alignment/751_634_vs_524_by_2.bam}
				is malformed: read ends with deletion. Cigar: 6M13I5M9D25M51I10D

	"""
	defaultGATKArguments = " --unsafe --validation_strictness SILENT --read_filter BadCigar "
	
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		self.pathToInsertHomePathList.extend(['samtools_path', 'picard_path', 'gatk_path', 'tabixPath'])
		#inserted before AbstractWorkflow.__init__()
		AbstractWorkflow.__init__(self, **keywords)
		#from pymodule import ProcessOptions
		#self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
		#												class_to_have_attr=self)
#		self.samtools_path = self.insertHomePath(self.samtools_path, self.home_path)
#		self.picard_path =  self.insertHomePath(self.picard_path, self.home_path)
#		self.gatk_path =  self.insertHomePath(self.gatk_path, self.home_path)
#		self.tabixPath =  self.insertHomePath(self.tabixPath, self.home_path)
		
		self.chr_pattern = Genome.chr_pattern
		self.contig_id_pattern = Genome.contig_id_pattern
		
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
	
	def registerCommonExecutables(self, workflow=None):
		"""
		2012.7.30
			become a link to registerExecutables()
		2012.7.25
			add noDefaultClustersSizeExecutableList
		2011-11-22
		"""
		self.registerExecutables(workflow=workflow)
	
	def registerExecutables(self, workflow=None):
		"""
		2012.8.7 remove noDefaultClustersSizeExecutableList, use executableClusterSizeMultiplierList instead
		2012.1.9 a symlink to registerCommonExecutables()
		"""
		AbstractWorkflow.registerExecutables(self)
		
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		
		selectAndSplit = Executable(namespace=namespace, name="SelectAndSplitAlignment", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		selectAndSplit.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "SelectAndSplitAlignment.py"), site_handler))
		executableClusterSizeMultiplierList.append((selectAndSplit, 0))
		
		selectAndSplitFasta = Executable(namespace=namespace, name="SelectAndSplitFastaRecords", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		selectAndSplitFasta.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "SelectAndSplitFastaRecords.py"), site_handler))
		executableClusterSizeMultiplierList.append((selectAndSplitFasta, 0))
		
		samtools = Executable(namespace=namespace, name="samtools", version=version, os=operatingSystem, arch=architecture, installed=True)
		samtools.addPFN(PFN("file://" + self.samtools_path, site_handler))
		executableClusterSizeMultiplierList.append((samtools, 0))
		
		
		addOrReplaceReadGroupsJava = Executable(namespace=namespace, name="addOrReplaceReadGroupsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		addOrReplaceReadGroupsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((addOrReplaceReadGroupsJava, 0))
		
		genotyperJava = Executable(namespace=namespace, name="genotyperJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		genotyperJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((genotyperJava, 0))
		#clustering is controlled by a separate parameter
		#genotyperJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		
		BuildBamIndexFilesJava = Executable(namespace=namespace, name="BuildBamIndexFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		BuildBamIndexFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((BuildBamIndexFilesJava, 0.5))
		
		#2012.9.21 same as BuildBamIndexFilesJava, but no clustering
		IndexMergedBamIndexJava =  Executable(namespace=namespace, name="IndexMergedBamIndexJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		IndexMergedBamIndexJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((IndexMergedBamIndexJava, 0))
		
		createSequenceDictionaryJava = Executable(namespace=namespace, name="createSequenceDictionaryJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		createSequenceDictionaryJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((createSequenceDictionaryJava, 0))
		
		DOCWalkerJava = Executable(namespace=namespace, name="DOCWalkerJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		#no clusters_size for this because it could run on a whole bam for hours
		#DOCWalkerJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		DOCWalkerJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((DOCWalkerJava, 0))
		
		VariousReadCountJava = Executable(namespace=namespace, name="VariousReadCountJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		VariousReadCountJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#no clusters_size for this because it could run on a whole bam for hours
		#VariousReadCountJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		executableClusterSizeMultiplierList.append((VariousReadCountJava, 0))
		
		
		MarkDuplicatesJava = Executable(namespace=namespace, name="MarkDuplicatesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		MarkDuplicatesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#MarkDuplicatesJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		executableClusterSizeMultiplierList.append((MarkDuplicatesJava, 0))
		
		SelectVariantsJava = Executable(namespace=namespace, name="SelectVariantsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		SelectVariantsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((SelectVariantsJava, 0.3))
		
		CombineVariantsJava = Executable(namespace=namespace, name="CombineVariantsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		CombineVariantsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((CombineVariantsJava, 0.3))
		
		CallVariantBySamtools = Executable(namespace=namespace, name="CallVariantBySamtools", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		CallVariantBySamtools.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/CallVariantBySamtools.sh"), site_handler))
		#clustering is controlled by a separate parameter
		#CallVariantBySamtools.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		executableClusterSizeMultiplierList.append((CallVariantBySamtools, 0))
		
		GenotypeCallByCoverage = Executable(namespace=namespace, name="GenotypeCallByCoverage", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		GenotypeCallByCoverage.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mapper/GenotypeCallByCoverage.py"), site_handler))
		executableClusterSizeMultiplierList.append((GenotypeCallByCoverage, 1))
		
		MergeGenotypeMatrix = Executable(namespace=namespace, name="MergeGenotypeMatrix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		MergeGenotypeMatrix.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reducer/MergeGenotypeMatrix.py"), site_handler))
		executableClusterSizeMultiplierList.append((MergeGenotypeMatrix, 1))
		
		bgzip_tabix = Executable(namespace=namespace, name="bgzip_tabix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		bgzip_tabix.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/bgzip_tabix.sh"), site_handler))
		#2012.8.8 bgzip_tabix runs really fast. so multiplier set to 4. 
		executableClusterSizeMultiplierList.append((bgzip_tabix, 4))
		
		vcf_convert = Executable(namespace=namespace, name="vcf_convert", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		vcf_convert.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_convert.sh"), site_handler))
		executableClusterSizeMultiplierList.append((vcf_convert, 1))
		
		vcf_isec = Executable(namespace=namespace, name="vcf_isec", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_isec.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_isec.sh"), site_handler))
		executableClusterSizeMultiplierList.append((vcf_isec, 1))
		
		vcf_concat = Executable(namespace=namespace, name="vcf_concat", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_concat.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		executableClusterSizeMultiplierList.append((vcf_concat, 1))
		
		vcfSubsetPath = os.path.join(self.home_path, "bin/vcftools/vcf-subset")
		self.vcfSubsetPath = vcfSubsetPath	#vcfSubsetPath is first argument to vcfSubsetPath
		vcfSubset = Executable(namespace=namespace, name="vcfSubset", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcfSubset.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcfSubset.sh"), site_handler))
		executableClusterSizeMultiplierList.append((vcfSubset, 1))
		
		concatGATK = Executable(namespace=namespace, name="concatGATK", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		concatGATK.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		executableClusterSizeMultiplierList.append((concatGATK, 1))
		
		concatSamtools = Executable(namespace=namespace, name="concatSamtools", version=version, \
									os=operatingSystem, arch=architecture, installed=True)
		concatSamtools.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		executableClusterSizeMultiplierList.append((concatSamtools, 1))
		
		vcftoolsPath = os.path.join(self.home_path, "bin/vcftools/vcftools")
		self.vcftoolsPath = vcftoolsPath	#vcftoolsPath is first argument to vcftoolsWrapper
		vcftoolsWrapper = Executable(namespace=namespace, name="vcftoolsWrapper", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcftoolsWrapper.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcftoolsWrapper.sh"), site_handler))
		vcftoolsWrapper.vcftoolsPath = vcftoolsPath
		executableClusterSizeMultiplierList.append((vcftoolsWrapper, 1))
		#vcftoolsWrapper.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		#self.addExecutable(vcftoolsWrapper)
		#self.vcftoolsWrapper = vcftoolsWrapper
		
		tabix = Executable(namespace=namespace, name="tabix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		tabix.addPFN(PFN("file://" + self.tabixPath, site_handler))
		executableClusterSizeMultiplierList.append((tabix,5))
		
		#2011.12.21 moved from FilterVCFPipeline.py
		FilterVCFByDepthJava = Executable(namespace=namespace, name="FilterVCFByDepthJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		FilterVCFByDepthJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((FilterVCFByDepthJava, 1 ))
		
		#2011.12.21	for OutputVCFSiteStat.py
		tabixRetrieve = Executable(namespace=namespace, name="tabixRetrieve", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		tabixRetrieve.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/tabixRetrieve.sh"), site_handler))
		executableClusterSizeMultiplierList.append((tabixRetrieve, 1 ))
		
		#2012.3.1
		MergeFiles = Executable(namespace=namespace, name="MergeFiles", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		MergeFiles.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "shell/MergeFiles.sh"), site_handler))
		executableClusterSizeMultiplierList.append((MergeFiles, 0 ))
		
		#2012.7.25
		MergeVCFReplicateHaplotypesJava= Executable(namespace=namespace, name="MergeVCFReplicateHaplotypesJava", \
											version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		MergeVCFReplicateHaplotypesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((MergeVCFReplicateHaplotypesJava, 0.5))
		
		#2012.7.29, moved from CheckTwoVCFOverlapPipeline.py
		CheckTwoVCFOverlap = Executable(namespace=namespace, name="CheckTwoVCFOverlap", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		CheckTwoVCFOverlap.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/CheckTwoVCFOverlap.py"), site_handler))
		executableClusterSizeMultiplierList.append((CheckTwoVCFOverlap, 1))
		
		#2012.9.6
		AppendInfo2SmartPCAOutput = Executable(namespace=namespace, name="AppendInfo2SmartPCAOutput", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		AppendInfo2SmartPCAOutput.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/AppendInfo2SmartPCAOutput.py"), site_handler))
		executableClusterSizeMultiplierList.append((AppendInfo2SmartPCAOutput, 0))
		
				
		AddAlignmentFile2DB = Executable(namespace=namespace, name="AddAlignmentFile2DB", version=version, os=operatingSystem,\
										arch=architecture, installed=True)
		AddAlignmentFile2DB.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "db/AddAlignmentFile2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((AddAlignmentFile2DB, 0))
		
		
		mergeSamFilesJava = Executable(namespace=namespace, name="mergeSamFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		mergeSamFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((mergeSamFilesJava, 0))
		
		SortSamFilesJava = Executable(namespace=namespace, name="SortSamFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		SortSamFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((SortSamFilesJava, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	bwaIndexFileSuffixLs = ['amb', 'ann', 'bwt', 'pac', 'sa']
	#, 'nhr', 'nin', 'nsq' are formatdb (blast) output, 2012.10.18 i think
	
	def registerBWAIndexFile(self, refFastaFname=None, input_site_handler=None, folderName=""):
		"""
		2012.10.10
		"""
		if input_site_handler is None:
			input_site_handler = self.input_site_handler
		return yh_pegasus.registerRefFastaFile(self, refFastaFname, registerAffiliateFiles=True, \
									input_site_handler=input_site_handler,\
									checkAffiliateFileExistence=True, addPicardDictFile=False, \
									affiliateFilenameSuffixLs=self.bwaIndexFileSuffixLs,\
									folderName=folderName)
	
	def addBAMIndexJob(self, workflow=None, BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					inputBamF=None,\
					extraArguments=None, parentJobLs=None, extraDependentInputLs=None, \
					transferOutput=True, javaMaxMemory=2500,\
					**keywords):
		"""
		2012.10.18 use addGenericJob() instead
		2012.4.12
			remove argument parentJob and stop adding it to parentJobLs, which causes an insidious bug
				that accumulates parent jobs from multiple calls of addBAMIndexJob() into parentJobLs
				(they all become parents of this bam index job.)
		2012.3.22
			bugfix, change argument parentJobLs's default value from [] to None. [] would make every run have the same parentJobLs 
			proper transfer/register setup
		2011-11-20
		"""
		memRequirementData = self.getJVMMemRequirment(job_max_memory=javaMaxMemory, minMemory=2000)
		job_max_memory = memRequirementData.memRequirement
		javaMemRequirement = memRequirementData.memRequirementInStr
		baiFile = File('%s.bai'%inputBamF.name)		
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', BuildBamIndexFilesJar,\
							"VALIDATION_STRINGENCY=LENIENT", \
							"INPUT=", inputBamF, "OUTPUT=", baiFile]
					#not including 'SORT_ORDER=coordinate'
					#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexFilesJar won't fail.)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputBamF])
		
		job= self.addGenericJob(executable=BuildBamIndexFilesJava, inputFile=None,\
							outputFile=None, outputArgumentOption="-o", \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraOutputLs=[baiFile],\
							transferOutput=transferOutput, \
							extraArgumentList=extraArgumentList, \
							job_max_memory=memRequirementData.memRequirement, **keywords)
		job.bamFile = inputBamF
		job.baiFile = baiFile
		#job.parentJobLs is where the actual alignment job and its bam/sam output are.
		return job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.1.9
			abstract function
		"""
		AbstractWorkflow.registerCustomExecutables(self, workflow=workflow)
	
	
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
		self.no_of_jobs += 1
		return fastaIndexJob
	
	def addRefFastaDictJob(self, workflow=None, createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, refFastaF=None):
		"""
		2012.9.14 bugfix. add argument createSequenceDictionaryJar
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
		self.no_of_jobs += 1
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
							namespace=None, version=None, transferOutput=False, parentJobLs=None, tabixArguments=None,\
							extraDependentInputLs=None, **keywords):
		"""
		2012.8.17 if transferOutput is None, do not register output files as OUTPUT with transfer flag
			use addGenericJob()
		2011.12.20
			pass additional tabix arguments to bgzip_tabix shell script
		2011-11-4
		
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = []
		extraOutputLs = []
		key2ObjectForJob = {}
		
		tbi_F = File("%s.tbi"%outputF.name)
		key2ObjectForJob['tbi_F'] = tbi_F
		extraOutputLs.append(tbi_F)
		# 2012.8.19 add the parentJob to parentJobLs
		if parentJobLs is None:
			parentJobLs = []
		if parentJob:
			parentJobLs.append(parentJob)
		
		job = self.addGenericJob(executable=bgzip_tabix, inputFile=inputF, inputArgumentOption="", \
					outputFile=outputF, outputArgumentOption="", \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=tabixArguments, extraArgumentList=extraArgumentList, job_max_memory=2000,  sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, **keywords)
		return job
	
	def addVCFConcatJob(self, workflow=None, concatExecutable=None, parentDirJob=None, outputF=None, parentJobLs=None, \
					extraDependentInputLs =None, transferOutput=True, vcf_job_max_memory=500, **keywords):
		"""
		2012.8.30
			use addGenericJob() instead
		2011-11-5
		"""
		#2011-9-22 union of all samtools intervals for one contig
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = []
		extraOutputLs = []
		key2ObjectForJob = {}
		
		tbi_F = File("%s.tbi"%outputF.name)
		key2ObjectForJob['tbi_F'] = tbi_F
		extraOutputLs.append(tbi_F)
		# 2012.8.19 add the parentJob to parentJobLs
		if parentJobLs is None:
			parentJobLs = []
		if parentDirJob:
			parentJobLs.append(parentDirJob)
		
		job = self.addGenericJob(executable=concatExecutable, inputFile=None, inputArgumentOption="", \
					outputFile=outputF, outputArgumentOption="", \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=vcf_job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, **keywords)
		return job
	
	def addVCFSubsetJob(self, workflow=None, executable=None, vcfSubsetPath=None, sampleIDFile=None,\
					inputVCF=None, outputF=None, \
					parentJobLs=None, transferOutput=True, job_max_memory=200,\
					extraArguments=None, extraDependentInputLs=None, **keywords):
		"""
		2012.10.10 use addGenericJob
		2012.5.10
		"""
		extraArgumentList = [vcfSubsetPath, sampleIDFile, inputVCF, outputF]
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraDependentInputLs.append(sampleIDFile)
		extraDependentInputLs.append(inputVCF)
		extraOutputLs = [outputF]
		job = self.addGenericJob(executable=executable, inputFile=None, \
					outputFile=None, \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					job_max_memory=job_max_memory, sshDBTunnel=None, \
					key2ObjectForJob=None, **keywords)
		return job

	def addVCF2MatrixJob(self, workflow=None, executable=None, inputVCF=None, outputFile=None, \
						refFastaF=None, run_type=3, numberOfReadGroups=10, seqCoverageF=None, \
						outputDelimiter=None, minDepth=0, \
						parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.9.24 added argument minDepth
		2012.8.20
			add argument outputDelimiter and use addGenericJob()
		2012.5.8
			executable is GenotypeCallByCoverage
		"""
		
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		#2012.5.8 "-n 10" means numberOfReadGroups is 10 but it's irrelevant when "-y 3" (run_type =3, read from vcf without filter)
		extraArgumentList = [ "-n 10", '-y 3']
		extraOutputLs = []
		key2ObjectForJob = {}
		
		if refFastaF:
			extraArgumentList.extend(["-e", refFastaF])
			extraDependentInputLs.append(refFastaF)
		if seqCoverageF:
			extraArgumentList.extend(["-q", seqCoverageF])
			extraDependentInputLs.append(seqCoverageF)
		if outputDelimiter:
			extraArgumentList.append('-u %s'%(outputDelimiter))
		if minDepth is not None and minDepth>=0:
			extraArgumentList.append("--minDepth %s"%(minDepth))
			
		job = self.addGenericJob(executable=executable, inputFile=inputVCF, inputArgumentOption="-i", \
					outputFile=outputFile, outputArgumentOption="-o", \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, job_max_memory=2000,  sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, **keywords)
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
		self.no_of_jobs += 1
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
					refFastaFList=[], sampleIDKeepFile=None, snpIDKeepFile=None, sampleIDExcludeFile=None, \
					parentJobLs=None, \
					extraDependentInputLs=None, transferOutput=True, extraArguments=None, job_max_memory=2000, interval=None,\
					**keywords):
		"""
		2012.10.17 add argument sampleIDKeepFile, snpIDKeepFile, sampleIDExcludeFile
		2012.10.10 use addGenericJob()
		2012.10.5 try add new option "--regenotype" (to extraArguments) to allow re-genotype the selected samples based on their GLs (or PLs)
			does it update the DP INFO field.
		2011-12.5
			http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_SelectVariants.html
			
			option:
				--concordance	RodBinding[VariantContext]	none	Output variants that were also called in this comparison track
				--discordance	RodBinding[VariantContext]	none	Output variants that were not called in this comparison track
				--exclude_sample_file	Set[File]	[]	File containing a list of samples (one per line) to exclude.
						Can be specified multiple times
				--exclude_sample_name	Set[String]	[]	Exclude genotypes from this sample. Can be specified multiple times
				--excludeFiltered	boolean	false	Don't include filtered loci in the analysis
				--excludeNonVariants	boolean	false	Don't include loci found to be non-variant after the subsetting procedure
				--keepIDs	File	NA	Only emit sites whose ID is found in this file (one ID per line)
				--keepOriginalAC	boolean	false	Don't update the AC, AF, or AN values in the INFO field after selecting
				--mendelianViolation	Boolean	false	output mendelian violation sites only
				-mvq	double	0.0	Minimum genotype QUAL score for each trio member required to accept a site as a violation
				--regenotype	Boolean	false	re-genotype the selected samples based on their GLs (or PLs)
				--remove_fraction_genotypes	double	0.0	Selects a fraction (a number between 0 and 1) of the total
						genotypes at random from the variant track and sets them to nocall
				--restrictAllelesTo	NumberAlleleRestriction	ALL	Select only variants of a particular allelicity.
						Valid options are ALL (default), MULTIALLELIC or BIALLELIC
				--sample_expressions	Set[String]	NA	Regular expression to select many samples from the ROD tracks provided.
					Can be specified multiple times
				--sample_file	Set[File]	NA	File containing a list of samples (one per line) to include.
					Can be specified multiple times
				--sample_name	Set[String]	[]	Include genotypes from this sample. Can be specified multiple times
				--select_expressions	ArrayList[String]	[]	One or more criteria to use when selecting the data
				--select_random_fraction	double	0.0	Selects a fraction (a number between 0 and 1) of the total
						variants at random from the variant track
				--select_random_number	int	0	Selects a number of variants at random from the variant track
				--selectTypeToInclude	List[Type]	[]	Select only a certain type of variants from the input file. 
						Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION. Can be specified multiple times
				 
		"""
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		
		refFastaF = refFastaFList[0]
		extraArgumentList = [javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T SelectVariants", "-R", refFastaF, \
						"--variant", inputF, "-o ", outputF]
		if interval:
			extraArgumentList.append("-L %s"%(interval))
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraDependentInputLs.append(inputF)
		extraDependentInputLs = extraDependentInputLs + refFastaFList
		if sampleIDKeepFile:
			extraArgumentList.extend(['--sample_file', sampleIDKeepFile])
			extraDependentInputLs.append(sampleIDKeepFile)
		if snpIDKeepFile:
			extraArgumentList.extend(['--keepIDs', snpIDKeepFile])
			extraDependentInputLs.append(snpIDKeepFile)
		if sampleIDExcludeFile:
			extraArgumentList.extend(['--exclude_sample_file', sampleIDExcludeFile])
			extraDependentInputLs.append(sampleIDExcludeFile)
		
		#register the two idx files so they will be cleaned out
		extraOutputLs = [outputF, File('%s.idx'%(outputF.name)), File('%s.idx'%(inputF.name))]
		
		job = self.addGenericJob(executable=SelectVariantsJava, inputFile=None,\
					outputFile=None, \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					job_max_memory=job_max_memory, sshDBTunnel=None, \
					key2ObjectForJob=None, **keywords)
		return job
	
	def addSelectAlignmentJob(self, executable=None, inputFile=None, \
							outputFile=None, region=None, parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
							extraArguments=None, job_max_memory=2000, needBAMIndexJob=True, **keywords):
		"""
		2012.9.17 copied from vervet/src/AlignmentReadBaseQualityRecalibrationWorkflow.py
		2012.6.27
		"""
		#select reads that are aligned to one region
		extraArgumentList = ['view', '-h', "-b", "-u", "-o", outputFile, inputFile, region]	# -b -u forces uncompressed bam output
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[inputFile]
		else:
			extraDependentInputLs.append(inputFile)
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[outputFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		if needBAMIndexJob:
			# add the index job on the bam file
			bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
						BuildBamIndexFilesJar=self.BuildBamIndexFilesJar, \
						inputBamF=job.output, parentJobLs=[job], \
						transferOutput=transferOutput, job_max_memory=job_max_memory)
		else:
			bamIndexJob = None
		return job, bamIndexJob
	
	def getVCFFileID2path(self, inputDir):
		"""
		2011-12-1
		"""
		sys.stderr.write("Getting all vcf files from %s "%(inputDir))
		vcfFileID2path = {}
		for inputFname in os.listdir(inputDir):
			inputAbsPath = os.path.join(os.path.abspath(inputDir), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False) and not NextGenSeq.isVCFFileEmpty(inputAbsPath):
				fileID = Genome.getChrFromFname(inputFname)
				if fileID in vcfFileID2path:
					sys.stderr.write("fileID %s already has value (%s) in dictionary. but now a 2nd file %s overwrites previous value.\n"%\
									(fileID, vcfFileID2path.get(fileID), inputFname))
				vcfFileID2path[fileID] = inputAbsPath
		sys.stderr.write("  found %s files.\n"%(len(vcfFileID2path)))
		return vcfFileID2path
	
	def addCheckTwoVCFOverlapJob(self, workflow=None, executable=None, vcf1=None, vcf2=None, chromosome=None, chrLength=None, \
					outputF=None, outputFnamePrefix=None, parentJobLs=None, \
					extraDependentInputLs=None, transferOutput=False, extraArguments=None, job_max_memory=1000, \
					perSampleMatchFraction=False, **keywords):
		"""
		2012.8.16
			now perSampleMatchFraction output is in a separate output file.
		2011.12.9
		"""
		if outputF is None and outputFnamePrefix:
			outputF = File('%s.tsv'%(outputFnamePrefix))
		extraOutputLs = []
		extraArgumentList = []
		key2ObjectForJob = {}
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 		
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		if vcf2:
			extraArgumentList.extend(["-j", vcf2])
			extraDependentInputLs.append(vcf2)
		if chromosome:
			extraArgumentList.extend(["-c", chromosome])
		if chrLength:
			extraArgumentList.append("-l %s"%(chrLength))
		if perSampleMatchFraction and outputFnamePrefix:
			extraArgumentList.append("-p")
			suffixAndNameTupleList.append(['_perSample.tsv','perSample'])
		if outputFnamePrefix:
			extraArgumentList.extend(['-O', outputFnamePrefix])
			suffixAndNameTupleList.append(['_overlapSitePos.tsv','overlapSitePos'])
		if extraArguments:
			extraArgumentList.append(extraArguments)

		self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, suffixAndNameTupleList=suffixAndNameTupleList, \
													extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
		job= self.addGenericJob(executable=executable, inputFile=vcf1, outputFile=outputF, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
						key2ObjectForJob=key2ObjectForJob, **keywords)
		#if outputFnamePrefix:
		#	job.overlapSitePosF = overlapSitePosF
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
							snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
							perIndividualDepth=False, perIndividualHeterozygosity=False, \
							perSiteHWE=False, haploLD=False, genoLD=False, minLDr2=0, LDWindowByNoOfSites=None,\
							LDWindowByBP=None, TsTvWindowSize=None, piWindowSize=None, perSitePI=False, \
							SNPDensityWindowSize=None, calculateMissingNess=False, calculateFreq=False, calculateFreq2=False,\
							getSiteDepth=False, getSiteMeanDepth=False, getSiteQuality=False,\
							outputFormat='--recode', parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
							extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.8.1
		2012.5.9
			moved from FilterVCFPipeline.py
			add argument outputFormat to make "--recode" by default and also allow other output formats.
			add argument extraArguments to accept something like "--recode-INFO-all".
			
			this could be just used to output vcf in various formats
			
		2011-11-21
			argument vcftools is replaced with a wrapper, which takes vcftools path as 1st argument
			outputFormat="--recode" instructs vcftools to output a VCF file. Without it, "--recode-INFO-all" will do nothing.
				"--recode-INFO-all" is added to vcftools to output input VCF also in VCF format and recalculate all the INFO fields.
				OR "--recode-INFO string" options only keeps key=string in the INFO field.
				The latter two arguments should be added through extraArguments.
				i.e.
				vcf1KeepGivenSNPByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
								inputVCFF=vcf1, \
								outputFnamePrefix=outputFnamePrefix, \
								parentJobLs=[vcf1_vcftoolsFilterDirJob], \
								snpMisMatchStatFile=keepSNPPosF, \
								minMAC=None, minMAF=None, \
								maxSNPMissingRate=None,\
								extraDependentInputLs=[vcf1.tbi_F], outputFormat='--recode', extraArguments="--recode-INFO-all")
				 
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraDependentInputLs.append(inputVCFF)
		extraArgumentList = [self.vcftoolsPath]
		extraOutputLs = []
		key2ObjectForJob = {}
		if inputVCFF.name[-2:]=='gz':
			extraArgumentList.extend(["--gzvcf", inputVCFF])
		else:
			extraArgumentList.extend(["--vcf", inputVCFF])
		#filter options
		if snpMisMatchStatFile:
			extraArgumentList.extend(["--positions", snpMisMatchStatFile])
			extraDependentInputLs.append(snpMisMatchStatFile)
		if maxSNPMissingRate is not None:
			extraArgumentList.append("--geno %s"%(1-maxSNPMissingRate))
		if minMAF is not None:
			extraArgumentList.append("--maf %s"%(minMAF))
		if minMAC is not None:
			extraArgumentList.append("--mac %s"%(minMAC))
		
		if outputFnamePrefix:
			extraArgumentList.extend(["--out", outputFnamePrefix])
		
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 
		
		#prime output, added into extraOutputLs, ahead of all others.
		if outputFormat:
			extraArgumentList.append(outputFormat)
		
		if outputFormat=='--recode':	#2012.5.9
			outputVCFF = File("%s.recode.vcf"%(outputFnamePrefix))
			extraOutputLs.insert(0, outputVCFF)
				#since outputFile passed to addGenericJob is None, 1st entry of extraOutputLs
				# becomes job.output
		elif outputFormat=='--plink-tped':
			suffixAndNameTupleList.extend([['.tped'], ['.tfam']])
		elif outputFormat=='--plink':
			suffixAndNameTupleList.extend([['.ped'], ['.fam']])
		elif outputFormat=='--IMPUTE':
			suffixAndNameTupleList.extend([['.impute.hap'], ['.impute.hap.legend'], ['.impute.hap.indv']])
		elif outputFormat=='--BEAGLE-GL':
			suffixAndNameTupleList.extend([['.BEAGLE.GL']])
			#output1 = File("%s.BEAGLE.GL"%(outputFnamePrefix))
			#extraOutputLs.extend([output1, output2])
		
		#leave the log file in workflow folder, not staging it out
		vcftoolsLogFile = File('%s.log'%(outputFnamePrefix))
		#2nd-tier output, mostly stats output		
		if perIndividualDepth:
			"""
			output example:
			*.idepth
INDV    N_SITES MEAN_DEPTH
1511_639_1987079_GA_vs_524      2483    21.3363
1512_640_1985088_GA_vs_524      2483    28.8647
1513_641_1986014_GA_vs_524      2483    29.9299
1514_642_1988009_GA_vs_524      2483    27.7483
1515_688_1988086_GA_vs_524      2483    16.1776

			"""
			extraArgumentList.append('--depth')
			suffixAndNameTupleList.extend([['.idepth']])
		if perIndividualHeterozygosity:
			"""
			output example:
INDV    O(HOM)  E(HOM)  N_SITES F
1511_639_1987079_GA_vs_524      1741    1778.7  2475    -0.05418
1512_640_1985088_GA_vs_524      1589    1778.7  2475    -0.27248
1513_641_1986014_GA_vs_524      1588    1778.7  2475    -0.27392
1514_642_1988009_GA_vs_524      1423    1778.7  2475    -0.51089
1515_688_1988086_GA_vs_524      2455    1778.7  2475    0.97128

			"""
			extraArgumentList.append('--het')
			suffixAndNameTupleList.extend([['.het'], ])
		if perSiteHWE:
			"""
			output example:
CHR     POS     OBS(HOM1/HET/HOM2)      E(HOM1/HET/HOM2)        ChiSq   P
Contig966       203     4/8/4   4.00/8.00/4.00  0.000000        1.000000
Contig966       570     15/1/0  15.02/0.97/0.02 0.016649        1.000000
Contig966       1462    5/8/3   5.06/7.88/3.06  0.004031        1.000000
Contig966       3160    15/1/0  15.02/0.97/0.02 0.016649        1.000000
Contig966       3311    15/1/0  15.02/0.97/0.02 0.016649        1.000000
Contig966       3539    15/1/0  15.02/0.97/0.02 0.016649        1.000000


			"""
			extraArgumentList.append('--hardy')
			suffixAndNameTupleList.extend([['.hwe']])
		if haploLD or genoLD:
			if haploLD:
				extraArgumentList.append('--hap-r2')
				output = File('%s.hap.ld'%(outputFnamePrefix))
			elif genoLD:
				"""
				example:
CHR     POS1    POS2    N_INDV  R^2
Contig966       203     1462    16      0.790323
Contig966       203     4101    16      0.790323
Contig966       203     4200    16      0.790323
Contig966       203     4573    16      1
Contig966       203     4984    16      0.882883
Contig966       203     5289    16      0.790323
Contig966       203     7383    16      0.790323
Contig966       203     7403    16      0.882883
Contig966       203     8129    16      0.882883
Contig966       203     8453    16      0.333333
Contig966       203     8508    16      1
Contig966       203     8655    16      0.790323
Contig966       203     8817    16      0.790323
Contig966       203     9862    16      0.790323
Contig966       203     10439   16      1

				"""
				extraArgumentList.append('--geno-r2')
				output = File('%s.geno.ld'%(outputFnamePrefix))
			extraArgumentList.append('--min-r2 %s'%(minLDr2))
			if LDWindowByBP:
				extraArgumentList.append("--ld-window-bp %s"%(LDWindowByBP))
			elif LDWindowByNoOfSites:
				extraArgumentList.append("--ld-window %s"%(LDWindowByNoOfSites))
			extraOutputLs.append(output)
			key2ObjectForJob['ldFile'] = output
		
		if TsTvWindowSize:
			"""
			output example:
			*.TsTv
CHROM   BinStart        SNP_count       Ts/Tv
Contig459       0       96      1.18182
Contig459       20000   118     1.40816
Contig459       40000   97      1.30952
Contig459       60000   91      0.857143
Contig459       80000   91      1.275
Contig459       100000  109     0.786885
Contig459       120000  92      1.55556

			*.TsTv.summary
MODEL   COUNT
AC      1948
AG      3139
AT      696
CG      899
CT      3230
GT      1935
Ts      6369
Tv      5478

			"""
			extraArgumentList.append("--TsTv %s"%(TsTvWindowSize))
			TsTvFile = File('%s.TsTv'%(outputFnamePrefix))
			TsTvSummaryFile = File('%s.TsTv.summary'%(outputFnamePrefix))
			extraOutputLs.extend([TsTvFile, TsTvSummaryFile])
			key2ObjectForJob['TsTvFile'] = TsTvFile
			key2ObjectForJob['TsTvSummaryFile'] = TsTvSummaryFile
		if piWindowSize:
			"""
			output example:
CHROM   BIN_START       N_SNPS  PI
Contig966       0       56      16.9133
Contig966       20000   80      19.1976
Contig966       40000   48      14.9133
Contig966       60000   132     28.8085

			"""
			extraArgumentList.append("--window-pi %s"%(piWindowSize))
			windowPIFile = File('%s.windowed.pi'%(outputFnamePrefix))
			extraOutputLs.append(windowPIFile)
			key2ObjectForJob['windowPIFile'] = windowPIFile
		if perSitePI:
			"""
			example:
CHROM   POS     PI
Contig966       203     0.516129
Contig966       570     0.0625
Contig966       1462    0.508065
Contig966       3160    0.0625
Contig966       3311    0.0625
Contig966       3539    0.0625
Contig966       4101    0.508065
Contig966       4200    0.508065

			"""
			extraArgumentList.append("--site-pi")
			sitePIFile = File('%s.sites.pi'%(outputFnamePrefix))
			extraOutputLs.append(sitePIFile)
			key2ObjectForJob['sitePIFile'] = sitePIFile
		if SNPDensityWindowSize:
			"""
			output example:
CHROM   BIN_START       SNP_COUNT       SNPS/KB
Contig966       0       56      2.8
Contig966       20000   80      4
Contig966       40000   50      2.5
Contig966       60000   132     6.6
Contig966       80000   53      2.65

			"""
			extraArgumentList.append("--SNPdensity %s"%(SNPDensityWindowSize))
			snpDensityFile = File('%s.snpden'%(outputFnamePrefix))
			extraOutputLs.append(snpDensityFile)
			key2ObjectForJob['snpDensityFile'] = snpDensityFile
		if calculateMissingNess:
			"""
			output example:
			*.imiss
INDV    N_DATA  N_GENOTYPES_FILTERED    N_MISS  F_MISS
1511_639_1987079_GA_vs_524      2483    0       0       0
1512_640_1985088_GA_vs_524      2483    0       0       0
1513_641_1986014_GA_vs_524      2483    0       0       0
1514_642_1988009_GA_vs_524      2483    0       0       0
			
			*.lmiss
CHR     POS     N_DATA  N_GENOTYPE_FILTERED     N_MISS  F_MISS
Contig966       203     32      0       0       0
Contig966       570     32      0       0       0
Contig966       1462    32      0       0       0
Contig966       3160    32      0       0       0

			"""
			extraArgumentList.append("--missing")
			imissFile = File('%s.imiss'%(outputFnamePrefix))
			lmissFile = File('%s.lmiss'%(outputFnamePrefix))
			extraOutputLs.extend([imissFile, lmissFile])
			key2ObjectForJob['imissFile'] = imissFile
			key2ObjectForJob['lmissFile'] = lmissFile
		if calculateFreq or calculateFreq2:
			if calculateFreq:
				"""
				output example:
CHROM   POS     N_ALLELES       N_CHR   {ALLELE:FREQ}
Contig966       203     2       32      A:0.5   G:0.5
Contig966       570     2       32      C:0.96875       A:0.03125
Contig966       1462    2       32      T:0.5625        A:0.4375
Contig966       3160    2       32      A:0.96875       C:0.03125
Contig966       3311    2       32      G:0.96875       A:0.03125
Contig966       3539    2       32      C:0.96875       T:0.03125

				"""
				extraArgumentList.append("--freq")
			elif calculateFreq2:
				extraArgumentList.append("--freq2")
				"""
				output example:
				
CHROM   POS     N_ALLELES       N_CHR   {FREQ}
Contig966       203     2       32      0.5     0.5
Contig966       570     2       32      0.96875 0.03125
Contig966       1462    2       32      0.5625  0.4375

				"""
			freqFile = File('%s.frq'%(outputFnamePrefix))
			extraOutputLs.extend([freqFile])
			key2ObjectForJob['freqFile'] = freqFile
		if getSiteDepth:
			"""
			output example:
CHROM   POS     SUM_DEPTH       SUMSQ_DEPTH
Contig966       203     332     7590
Contig966       570     349     8751
Contig966       1462    119     1223
Contig966       3160    273     6331
Contig966       3311    327     7715

			"""
			extraArgumentList.append("--site-depth")
			outputFile = File('%s.ldepth'%(outputFnamePrefix))
			extraOutputLs.append(outputFile)
			key2ObjectForJob['ldepthFile'] = outputFile
		if getSiteMeanDepth:
			"""
			output example:
CHROM   POS     MEAN_DEPTH      VAR_DEPTH
Contig966       203     20.75   46.7333
Contig966       570     21.8125 75.8958
Contig966       1462    7.4375  22.5292
Contig966       3160    17.0625 111.529
Contig966       3311    20.4375 68.7958

			"""
			extraArgumentList.append("--site-mean-depth")
			outputFile = File('%s.ldepth.mean'%(outputFnamePrefix))
			extraOutputLs.append(outputFile)
			key2ObjectForJob['ldpethMeanFile'] = outputFile
		if getSiteQuality:
			"""
			output example:
CHROM   POS     QUAL
Contig966       203     999
Contig966       570     999
Contig966       1462    999
Contig966       3160    50

			"""
			extraArgumentList.append("--site-quality")
			outputFile = File('%s.lqual'%(outputFnamePrefix))
			extraOutputLs.append(outputFile)
			key2ObjectForJob['lqualFile'] = outputFile
			
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
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
		
		job= self.addGenericJob(executable=vcftoolsWrapper, inputFile=None, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
		return job
	
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
			if parentJob:
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
	
	def getReferenceSequence(self, workflow=None, **keywords):
		"""
		2013.1.25 placeholder, usually from database. such as:
			refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
			refFastaFname = os.path.join(self.data_dir, refSequence.path)
			refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
								input_site_handler=self.input_site_handler,\
								checkAffiliateFileExistence=True)
			return refFastaFList
			
		"""
		sys.stderr.write("Getting reference sequences (placeholder) ... nothing \n")
	
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
	
	def getTopNumberOfContigs(self, contigMaxRankBySize=100, contigMinRankBySize=1, tax_id=60711, sequence_type_id=9):
		"""
		2012.8.2
			moved from vervet/src/AlignmentToCallPipeline.py
			call GenomeDB.getTopNumberOfContigs() instead
		2011-11-6
			rename argument maxContigID to contigMaxRankBySize
			add argument contigMinRankBySize
		2011-9-13
			return chr2size instead of a set of ref names
		2011-7-12
			get all the top contigs
		"""
		no_of_contigs_to_fetch = contigMaxRankBySize-contigMinRankBySize+1
		sys.stderr.write("Getting %s contigs with rank (by size) between %s and %s  ..."%\
						(no_of_contigs_to_fetch, contigMinRankBySize, contigMaxRankBySize))
		chr2size = {}
		
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=contigMaxRankBySize, contigMinRankBySize=contigMinRankBySize, \
											tax_id=tax_id, sequence_type_id=sequence_type_id)
		return chr2size
	
	def restrictContigDictionry(self, dc=None, maxContigID=None, minContigID=None):
		"""
		2012.8.2
			if maxContigID is None or zero, no filter. same for minContigID.
		"""
		sys.stderr.write("Restricting a contig dictionary of size %s within minContigID=%s, maxContigID=%s, ... "%\
						(len(dc), maxContigID, minContigID))
		if (maxContigID is not None and maxContigID!=0) and (minContigID is not None and minContigID!=0):
			new_dc = {}
			for contig, data in dc.iteritems():
				try:
					contigID = int(self.getContigIDFromFname(contig))
					included = True
					if (maxContigID is not None and maxContigID!=0) and contigID>maxContigID:
						included = False
					if (minContigID is not None and minContigID!=0) and contigID<minContigID:
						included = False
					if included:
						new_dc[contig] = data
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			dc = new_dc
		sys.stderr.write(" %s contigs left.\n"%(len(dc)))
		return dc
	
	def getChr2IntervalDataLsBySplitBEDFile(self, intervalFname=None, noOfLinesPerUnit=2000, folderName=None, parentJobLs= None):
		"""
		2012.8.9 update it so that the interval encompassing all lines in one block/unit is known.
			good for mpileup to only work on that interval and then "bcftools view" select from sites from the block.
			TODO: offer partitioning by equal-chromosome span, rather than number of sites.
				Some sites could be in far from each other in one block, which could incur long-running mpileup. goal is to skip these deserts.
		2012.8.8 bugfix add -1 to the starting number below cuz otherwise it's included in the next block's start
				blockStopLineNumber = min(startLineNumber+(i+1)*noOfLinesPerUnit-1, stopLineNumber)	
		2012.7.30
			1. intervalFname is in BED format (tab/comma-delimited, chr start stop) and has to be sorted.
				start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.
			2. folderName is the relative path of the folder in the pegasus workflow, that holds intervalFname.
				it'll be created upon file stage-in. no mkdir job for it.
				
			get the number of lines in intervalFname.
			get chr2StartStopDataLsTuple
			for each chr, split its lines into units  that don't exceed noOfLinesPerUnit
				add the split job
				 
		"""
		sys.stderr.write("Splitting %s into blocks of %s lines ... "%(intervalFname, noOfLinesPerUnit))
		#from pymodule import utils
		#noOfLines = utils.getNoOfLinesInOneFileByWC(intervalFname)
		chr2StartStopDataLs = {}
		import csv
		from pymodule import figureOutDelimiter
		inf = open(intervalFname)
		reader = csv.reader(inf, delimiter=figureOutDelimiter(intervalFname))
		lineNumber = 0
		previousChromosome = None
		previousLine = None
		chromosome = None
		for row in reader:
			lineNumber += 1
			chromosome, start, stop = row[:3]
			start = int(start)	#0-based, starting base
			stop = int(stop)	#0-based, stopping base but not inclusive, i.e. [start, stop)
			
			if previousLine is None or chromosome!=previousLine.chromosome:	#first line or different chromosome
				if previousLine is not None and previousLine.chromosome is not None:
					
					prevChrLastStartStopData = chr2StartStopDataLs[previousLine.chromosome][-1]
					if prevChrLastStartStopData.stopLineNumber is None:
						prevChrLastStartStopData.stopLineNumber = previousLine.lineNumber
						prevChrLastStartStopData.stopLineStart = previousLine.start
						prevChrLastStartStopData.stopLineStop = previousLine.stop
					
				if chromosome not in chr2StartStopDataLs:
					StartStopData = PassingData(startLineNumber=lineNumber, startLineStart=start, startLineStop=stop, \
											stopLineNumber=None, stopLineStart=None, stopLineStop=None)
					chr2StartStopDataLs[chromosome] = [StartStopData]
			else:	#same chromosome and not first line
				lastStartStopData = chr2StartStopDataLs[chromosome][-1]
				if lastStartStopData.stopLineNumber is None:	#last block hasn't been closed yet.
					noOfLinesInCurrentBlock = lineNumber - lastStartStopData.startLineNumber +1
					if noOfLinesInCurrentBlock>=noOfLinesPerUnit:	#time to close it
						lastStartStopData.stopLineNumber = lineNumber
						lastStartStopData.stopLineStart = start
						lastStartStopData.stopLineStop = stop
				else:	#generate a new block
					StartStopData = PassingData(startLineNumber=lineNumber, startLineStart=start, startLineStop=stop, \
											stopLineNumber=None, stopLineStart=None, stopLineStop=None)
					chr2StartStopDataLs[chromosome].append(StartStopData)
			previousLine = PassingData(chromosome = chromosome, start=start, stop=stop, lineNumber=lineNumber)
		#final closure
		if previousLine is not None:	#intervalFname is not empty
			lastStartStopData = chr2StartStopDataLs[previousLine.chromosome][-1]
			if lastStartStopData.stopLineNumber is None:	#last block hasn't been closed yet.
				#close it regardless of whether it has enough lines in it or not.
				lastStartStopData.stopLineNumber = previousLine.lineNumber
				lastStartStopData.stopLineStart = previousLine.start
				lastStartStopData.stopLineStop = previousLine.stop
		sys.stderr.write("%s chromosomes out of %s lines.\n"%(len(chr2StartStopDataLs), lineNumber))
		
		intervalFile = self.registerOneInputFile(inputFname=intervalFname, folderName=folderName)
		chr2IntervalDataLs = {}
		counter = 0
		for chr, startStopDataLs in chr2StartStopDataLs.iteritems():
			for startStopData in startStopDataLs:
				blockStartLineNumber = startStopData.startLineNumber
				blockStopLineNumber = startStopData.stopLineNumber
				# 2012.8.9 the large interval that encompasses all BED lines 
				interval = '%s:%s-%s'%(chr, startStopData.startLineStart, startStopData.stopLineStop)
				blockIntervalFile = File(os.path.join(folderName, '%s_line_%s_%s_bed.tsv'%(chr, \
												blockStartLineNumber, blockStopLineNumber)))
				blockIntervalJob = self.addSelectLineBlockFromFileJob(executable=self.SelectLineBlockFromFile, \
						inputFile=intervalFile, outputFile=blockIntervalFile,\
						startLineNumber=blockStartLineNumber, stopLineNumber=blockStopLineNumber, \
						parentJobLs=parentJobLs, extraDependentInputLs=None, \
						transferOutput=False, job_max_memory=500)
				intervalFnameSignature = '%s_%s_%s'%(chr, blockStartLineNumber, blockStopLineNumber)
				if chr not in chr2IntervalDataLs:
					chr2IntervalDataLs[chr] = []
				intervalData = PassingData(file=blockIntervalFile, intervalFnameSignature=intervalFnameSignature, interval=interval,\
										overlapInterval=interval,\
										chr=chr, jobLs=[blockIntervalJob], job=blockIntervalJob)
				chr2IntervalDataLs[chr].append(intervalData)
				counter += 1
		sys.stderr.write("%s intervals and %s SelectLineBlockFromFile jobs.\n"%(counter, counter))
		return chr2IntervalDataLs
	
	def getChr2IntervalDataLsBySplitChrSize(self, chr2size=None, intervalSize=None, intervalOverlapSize=None):
		"""
		2012.7.30
		"""
		sys.stderr.write("Splitting %s references into intervals of %s bp (overlap=%s) ... "%(len(chr2size), intervalSize,\
																		intervalOverlapSize))
		chr2IntervalDataLs = {}
		counter =0
		for chr, refSize in chr2size.iteritems():
			no_of_intervals = max(1, int(math.ceil(refSize/float(intervalSize)))-1)
			for i in range(no_of_intervals):
				originalStartPos = i*intervalSize + 1
				#to render adjacent intervals overlapping because trioCaller uses LD
				startPos = max(1, originalStartPos-intervalOverlapSize)
				if i<no_of_intervals-1:
					originalStopPos = min((i+1)*intervalSize, refSize)
				else:	#last chunk, include bp till the end
					originalStopPos = refSize
				#to render adjacent intervals overlapping because trioCaller uses LD
				stopPos = min(refSize, originalStopPos+intervalOverlapSize)
				
				interval = "%s:%s-%s"%(chr, originalStartPos, originalStopPos)
				intervalFnameSignature = '%s_%s_%s'%(chr, originalStartPos, originalStopPos)
				overlapInterval = "%s:%s-%s"%(chr, startPos, stopPos)
				overlapIntervalFnameSignature = '%s_%s_%s'%(chr, startPos, stopPos)
				if chr not in chr2IntervalDataLs:
					chr2IntervalDataLs[chr] = []
				intervalData = PassingData(overlapInterval=overlapInterval, overlapIntervalFnameSignature=overlapIntervalFnameSignature,\
							interval=interval, intervalFnameSignature=intervalFnameSignature, \
							file=None,\
							chr=chr, start=originalStartPos,\
							stop=originalStopPos, overlapStart=startPos, overlapStop=stopPos, jobLs=[])
				chr2IntervalDataLs[chr].append(intervalData)
				counter += 1
		sys.stderr.write("%s intervals.\n"%(counter))
		return chr2IntervalDataLs
	
	def addPutStuffIntoDBJob(self, workflow=None, executable=None, inputFileLs=[], \
					logFile=None, commit=False, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, sshDBTunnel=0, **keywords):
		"""
		2012.5.8 add sshDBTunnel
		2012.4.3
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		job.addArguments("--drivername", self.drivername, "--hostname", self.hostname, "--dbname", self.dbname, \
						"--db_user", self.db_user, "--db_passwd", self.db_passwd)
		if extraArguments:
			job.addArguments(extraArguments)
		if commit:
			job.addArguments("--commit")
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
	
	def addAddAlignmentFile2DBJob(self, workflow=None, executable=None, inputFile=None, otherInputFileList=None,\
						individual_alignment_id=None, individual_sequence_id=None,\
						ref_sequence_id=None, \
						alignment_method_id=None,\
						parent_individual_alignment_id=None, \
						mask_genotype_method_id=None,\
						individual_sequence_file_raw_id=None,\
						format=None,\
						logFile=False, data_dir=None, \
						parentJobLs=None, \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=True, \
						job_max_memory=2000, sshDBTunnel=False, commit=True, **keywords):
		"""
		2012.9.20
			To specify individual_alignment:
				either individual_alignment_id or (parent_individual_alignment_id + mask_genotype_method_id)
				or others 
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = ['--logFilename', logFile]
		extraOutputLs = []
		key2ObjectForJob = {}
		if otherInputFileList:
			extraDependentInputLs.extend(otherInputFileList)
		if logFile:
			extraOutputLs.append(logFile)
		if individual_alignment_id:
			extraArgumentList.append("--individual_alignment_id %s"%(individual_alignment_id))
		if individual_sequence_id:
			extraArgumentList.append("--individual_sequence_id %s"%(individual_sequence_id))
		if ref_sequence_id:
			extraArgumentList.append("--ref_sequence_id %s"%(ref_sequence_id))
		if alignment_method_id:
			extraArgumentList.append("--alignment_method_id %s"%(alignment_method_id))
		if parent_individual_alignment_id:
			extraArgumentList.append("--parent_individual_alignment_id %s"%(parent_individual_alignment_id))
		if mask_genotype_method_id:
			extraArgumentList.append("--mask_genotype_method_id %s"%(mask_genotype_method_id))
		if individual_sequence_file_raw_id:
			extraArgumentList.append("--individual_sequence_file_raw_id %s"%(individual_sequence_file_raw_id))
		if format:
			extraArgumentList.append("--format %s"%(format))
		if data_dir:
			extraArgumentList.append("--data_dir %s"%(data_dir))
		if commit:
			extraArgumentList.append("--commit")
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		job.logFile = logFile
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		
		#add all input files to the last (after db arguments,) otherwise, it'll mask others (cuz these don't have options).
		if otherInputFileList:
			for inputFile in otherInputFileList:
				if inputFile:
					job.addArguments(inputFile)
		return job
	
	def addAlignmentMergeJob(self, workflow=None, AlignmentJobAndOutputLs=[], outputBamFile=None,samtools=None,\
					java=None, mergeSamFilesJava=None, mergeSamFilesJar=None, \
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					mv=None, parentJobLs=None, namespace='workflow', version='1.0', transferOutput=False,\
					**keywords):
		"""
		2012.9.17 copied from vervet/src/ShortRead2AlignmentPipeline.py
		2012.7.4 bugfix. add job dependency between alignmentJob and merge_sam_job after all have been added to the workflow.
		2012.3.29
			no more threads (only 2 threads at maximum and increase only 20% performance anyway).
			Some nodes' kernels can't handle threads properly and it leads to process hanging forever.
		2011-11-15
			MarkDuplicates will be run after this step. So outputBamFile no longer needs to be transferred out.
		2011-9-4
			add argument transferOutput, default=False, which means leave the output files where they are
		2011-8-28
			merge alignment
			index it
		"""
		if workflow is None:
			workflow = self
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=5000, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		namespace = getattr(workflow, 'namespace', namespace)
		version = getattr(workflow, 'version', version)
		
		if len(AlignmentJobAndOutputLs)>1:
			merge_sam_job = Job(namespace=namespace, name=mergeSamFilesJava.name, version=version)
			merge_sam_job.addArguments(javaMemRequirement, "-jar", mergeSamFilesJar, 'SORT_ORDER=coordinate', \
						'ASSUME_SORTED=true', 'OUTPUT=', outputBamFile, "VALIDATION_STRINGENCY=LENIENT")
			# 'USE_THREADING=true', threading might be causing process hanging forever (sleep).
			merge_sam_job.uses(outputBamFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			yh_pegasus.setJobProperRequirement(merge_sam_job, job_max_memory=job_max_memory)
			workflow.addJob(merge_sam_job)
			for AlignmentJobAndOutput in AlignmentJobAndOutputLs:
				alignmentJob, alignmentOutput = AlignmentJobAndOutput[:2]
				merge_sam_job.addArguments('INPUT=', alignmentOutput)
				merge_sam_job.uses(alignmentOutput, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=alignmentJob, child=merge_sam_job)
		else:	#one input file, no samtools merge. use "mv" to rename it instead. should use "cp", then the input would be cleaned by cleaning job.
			alignmentJob, alignmentOutput = AlignmentJobAndOutputLs[0][:2]
			merge_sam_job = Job(namespace=namespace, name=mv.name, version=version)
			merge_sam_job.addArguments(alignmentOutput, outputBamFile)
			merge_sam_job.uses(alignmentOutput, transfer=True, register=True, link=Link.INPUT)
			merge_sam_job.uses(outputBamFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			workflow.addJob(merge_sam_job)
			workflow.depends(parent=alignmentJob, child=merge_sam_job)	#2012.7.4
			sys.stderr.write(" copy (no merge, only one alignment) from %s to %s.\n"%(alignmentOutput.name, outputBamFile.name))
		#assign output
		merge_sam_job.output = outputBamFile
		#2012.9.21
		if parentJobLs:
			for parentJob in parentJobLs:
				workflow.depends(parent=parentJob, child=merge_sam_job)
						
		# add the index job on the merged bam file
		bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					inputBamF=outputBamFile,\
					parentJobLs=[merge_sam_job], namespace=namespace, version=version,\
					transferOutput=transferOutput, javaMaxMemory=3000)
		return merge_sam_job, bamIndexJob
	
	
	def addMergeVCFReplicateGenotypeColumnsJob(self, workflow=None, executable=None, genomeAnalysisTKJar=None, \
						inputF=None, outputF=None, replicateIndividualTag=None, \
						debugHaplotypeDistanceFile=None, \
						debugMajoritySupportFile=None,\
						refFastaFList=[], parentJobLs=[], extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=2000, analysis_type='MergeVCFReplicateHaplotypes',\
						**keywords):
		"""
		2012.8.15
			add argument analysis_type (could be MergeVCFReplicateHaplotypes, MergeVCFReplicateGenotypeColumns
		2012.7.25
			use self.addGenericJob() and moved from AlignmentToTrioCallPipeline.py
			added "-XX:MaxPermSize=1024m" jvm combat this error:
				java.lang.OutOfMemoryError: Java heap space
		2012.6.1
			change MergeVCFReplicateGenotypeColumns to MergeVCFReplicateHaplotypes
			
		2012.4.2
			java -jar /home/crocea/script/gatk/dist/GenomeAnalysisTK.jar -T MergeVCFReplicateGenotypeColumns 
				-R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
				--variant /tmp/Contig0.vcf -o /tmp/contig0_afterMerge.vcf --onlyKeepBiAllelicSNP --replicateIndividualTag copy
		"""
		#GATK job
		#MaxPermSize= min(35000, max(1024, job_max_memory*9/7))
		#javaMemRequirement = "-Xms%sm -Xmx%sm -XX:PermSize=%sm -XX:MaxPermSize=%sm"%(job_max_memory*95/100, job_max_memory, \
		#																			MaxPermSize*95/100, MaxPermSize)
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		refFastaF = refFastaFList[0]
		extraArgumentList = [javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", analysis_type,\
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