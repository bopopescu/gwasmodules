#!/usr/bin/env python
"""
2011-11-22
	a common class for pegasus workflows that work on NGS (next-gen sequencing) data
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import Executable, File, PFN, Link, Job
from pymodule import Genome, utils
from pymodule import ProcessOptions
from pymodule.Genome import IntervalData
from pymodule.yhio import NextGenSeq
from pymodule.pegasus import yh_pegasus
from pymodule.utils import PassingData
from pymodule.yhio.VCFFile import VCFFile
from pymodule.yhio.MatrixFile import MatrixFile
from pymodule.yhio.AlignmentDepthIntervalFile import AlignmentDepthIntervalFile
from pymodule.yhio.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey
from pymodule.algorithm.RBTree import RBDict
from AbstractWorkflow import AbstractWorkflow



parentClass = AbstractWorkflow
class AbstractNGSWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(parentClass.db_option_dict)
	
	option_default_dict.update({
						('ref_ind_seq_id', 1, int): [None, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("platypus_path", 1, ): ["%s/bin/Platypus/Platypus.py", '', 1, 'haplotype caller from http://www.well.ox.ac.uk/platypus'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/gatk/dist", '', 1, 'GATK folder containing its jar binaries'],\
						("gatk2_path", 1, ): ["%s/script/gatk2/", '', 1, 'GATK version 2 or afterwards, no more source code, just binary jar files.'],\
						('tabixPath', 1, ): ["%s/bin/tabix", '', 1, 'path to the tabix binary', ],\
						('bgzipPath', 1, ): ["%s/bin/bgzip", '', 1, 'path to the bgzip binary', ],\
						('vcftoolsPath', 1, ): ["%s/bin/vcftools", '', 1, 'path to the vcftools binary', ],\
						('vcfSubsetPath', 1, ): ["%s/bin/vcf-subset", '', 1, 'path to the vcf-subset program', ],\
						("ligateVcfPerlPath", 1, ): ["%s/bin/umake/scripts/ligateVcf.pl", '', 1, 'path to ligateVcf.pl'],\
						("vcfsorterPath", 1, ): ["%s/bin/vcfsorter.pl", '', 1, 'path to vcfsorter.pl, http://code.google.com/p/vcfsorter/'],\
						
						
						#to filter chromosomes
						('maxContigID', 0, int): [None, 'x', 1, 'if contig/chromosome(non-sex) ID > this number, it will not be included. If None or 0, no restriction.', ],\
						('minContigID', 0, int): [None, 'V', 1, 'if contig/chromosome(non-sex) ID < this number, it will not be included. If None or 0, no restriction.', ],\
						("contigMaxRankBySize", 1, int): [2500, 'N', 1, 'maximum rank (rank 1=biggest) of a contig/chr to be included in calling'],\
						("contigMinRankBySize", 1, int): [1, 'M', 1, 'minimum rank (rank 1=biggest) of a contig/chr to be included in calling'],\
						('chromosome_type_id', 0, int):[None, '', 1, 'what type of chromosomes to be included, same as table genome.chromosome_type.\n\
	0 or None: all, 1: autosomes, 2: X, 3:Y, 4: mitochondrial '],\
						('ref_genome_tax_id', 0, int):[60711, '', 1, 'used to fetch chromosome info from GenomeDB. column GenomeDB.AnnotAssembly.tax_id'],\
						('ref_genome_sequence_type_id', 0, int):[1, '', 1, 'used to fetch chromosome info from GenomeDB. column GenomeDB.SequenceType.id 1: assembledChromosome, 9: Scaffold'],\
						('ref_genome_version', 0, int):[None, '', 1, 'used to fetch chromosome info from GenomeDB. column GenomeDB.AnnotAssembly.version'],\
						('ref_genome_outdated_index', 0, int):[0, '', 1, 'used to fetch chromosome info from GenomeDB. 0 means not outdated. column GenomeDB.AnnotAssembly.outdated_index'],\
						
						('completedAlignment', 0, int):[None, '', 1, 'a flag requiring whether user chooses alignment that has been completed or not.\n\
	--completedAlignment 0 is same as --skipDoneAlignment. --completedAlignment 1 gets you only the alignments that has been completed. Default (None) has no effect.'],\
						('mask_genotype_method_id', 0, int):[None, '', 1, 'to filter alignments with this field'],\
						('skipDoneAlignment', 0, int):[0, '', 0, 'skip alignment whose db_entry is complete and affiliated file is valid\n\
	(for ShortRead2AlignmentWorkflow or AlignmentReadBaseQualityRecalibrationWorkflow)'],\
						('checkEmptyVCFByReading', 0, int):[0, 'E', 0, 'toggle to check if a vcf file is empty by reading its content'],\
						
						("needFastaIndexJob", 0, int): [0, 'A', 0, 'need to add a reference index job by samtools?'],\
						("needFastaDictJob", 0, int): [0, 'B', 0, 'need to add a reference dict job by picard CreateSequenceDictionary.jar?'],\
						
						#to filter alignment or individual_sequence
						("reduce_reads", 0, int): [None, '', 1, 'To filter which input alignments to fetch from db'],\
						('excludeContaminant', 0, int):[0, '', 0, 'toggle this to exclude alignments or sequences that are from contaminated individuals, \n\
		(IndividualSequence.is_contaminated=1)'],\
						("sequence_filtered", 0, int): [None, 'Q', 1, 'to filter alignments/individual_sequences. None: whatever; 0: unfiltered sequences, 1: filtered sequences: 2: ...'],\
						("site_id_ls", 0, ): ["", 'S', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
						("country_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
						("tax_id_ls", 0, ): ["60711", '', 1, 'comma/dash-separated list of taxonomy IDs. individuals must come from these taxonomies.'],\
						("sequence_type_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of IndividualSequence.sequence_type_id. Empty for no filtering'],\
						("sequencer_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of IndividualSequence.sequencer_id. Empty for no filtering'],\
						("sequence_batch_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of IndividualSequence.sequence_batch_id. Empty for no filtering'],\
						("version_ls", 0, ): ["", '', 1, 'comma/dash-separated list of IndividualSequence.version. Empty for no filtering'],\
						("sequence_max_coverage", 0, float): [None, '', 1, 'max IndividualSequence.coverage. Empty for no filtering'],\
						("sequence_min_coverage", 0, float): [None, '', 1, 'min IndividualSequence.coverage. Empty for no filtering'],\
						
						('alignmentDepthIntervalMethodShortName', 0, ): [None, '', 1, 'fetch intervals from AlignmentDepthIntervalFile table', ],\
						('minAlignmentDepthIntervalLength', 0, int): [1000, '', 1, 'minimum length for a alignment depth interval to be included', ],\
						('alignmentDepthMaxFold', 0, float): [2, '', 1, 'depth of an alignment depth interval have to be within a range of [1/foldChange, foldChange]*medianDepth', ],\
						('alignmentDepthMinFold', 0, float): [0.1, '', 1, 'depth of an alignment depth interval have to be within a range of [1/foldChange, foldChange]*medianDepth', ],\
						
						('intervalOverlapSize', 1, int): [300000, 'U', 1, 'overlap #bps/#loci between adjacent intervals from one contig/chromosome,\
				only used for TrioCaller, not for SAMtools/GATK', ],\
						('intervalSize', 1, int): [5000000, 'Z', 1, '#bps/#loci for adjacent intervals from one contig/chromosome (alignment or VCF)', ],\
						
						('defaultGATKArguments', 1, ): [" --unsafe ALL --validation_strictness SILENT --read_filter BadCigar ", '', 1, 'arguments that will be added to every GATK-related job', ],\
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
	
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		self.pathToInsertHomePathList.extend(['samtools_path', 'platypus_path', 'picard_path', 'gatk_path', 'tabixPath', \
									'bgzipPath', 'gatk2_path', 'ligateVcfPerlPath',\
									'vcftoolsPath', 'vcfSubsetPath', 'vcfsorterPath'])
		#inserted before parentClass.__init__()
		parentClass.__init__(self, **keywords)
		#from pymodule import ProcessOptions
		#self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
		#												class_to_have_attr=self)
#		self.samtools_path = self.insertHomePath(self.samtools_path, self.home_path)
#		self.picard_path =  self.insertHomePath(self.picard_path, self.home_path)
#		self.gatk_path =  self.insertHomePath(self.gatk_path, self.home_path)
#		self.tabixPath =  self.insertHomePath(self.tabixPath, self.home_path)
		
		self.chr_pattern = Genome.chr_pattern
		self.contig_id_pattern = Genome.contig_id_pattern
		
		self.needSplitChrIntervalData = True	#2013.06.21
	
	def extra__init__(self):
		"""
		2013.2.15
		"""
		parentClass.extra__init__(self)
		
		if hasattr(self, 'contigMaxRankBySize') and hasattr(self, 'contigMinRankBySize'):
			#2013.2.6 non-public schema dbs should be connected before the main vervetdb or other db (schema=public) is connected.
			self.chr2size = self.getTopNumberOfContigs(contigMaxRankBySize=self.contigMaxRankBySize, \
							contigMinRankBySize=self.contigMinRankBySize, tax_id=self.ref_genome_tax_id, \
							sequence_type_id=self.ref_genome_sequence_type_id, version=self.ref_genome_version, \
							chromosome_type_id=self.chromosome_type_id, outdated_index=self.ref_genome_outdated_index)
		else:
			self.chr2size = {}
		
		listArgumentName_data_type_ls = [("site_id_ls", int), ('country_id_ls', int), ('tax_id_ls', int),\
								('sequence_type_id_ls', int), ('sequencer_id_ls', int), \
								('sequence_batch_id_ls', int),('version_ls', int)]
		listArgumentName2hasContent = self.processListArguments(listArgumentName_data_type_ls, emptyContent=[])
		
		self.samtoolsExecutableFile = self.registerOneExecutableAsFile(path=self.samtools_path,\
													site_handler=self.input_site_handler)
		self.tabixExecutableFile = self.registerOneExecutableAsFile(path=self.tabixPath)
		self.bgzipExecutableFile = self.registerOneExecutableAsFile(path=self.bgzipPath)
		self.ligateVcfExecutableFile = self.registerOneExecutableAsFile(path=self.ligateVcfPerlPath)
		
		self.vcftoolsExecutableFile = self.registerOneExecutableAsFile(path=self.vcftoolsPath)
		self.vcfSubsetExecutableFile = self.registerOneExecutableAsFile(path=self.vcfSubsetPath)
		self.vcfsorterExecutableFile = self.registerOneExecutableAsFile(path=self.vcfsorterPath)
	
	
	def getAlignments(self, db=None):
		"""
		2013.04.03
			wrapper so that derivatives could call it easily
		"""
		if db is None:
			db = self.db
		alignmentLs = db.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, data_dir=self.data_dir,\
										individual_sequence_file_raw_id_type=self.individual_sequence_file_raw_id_type,\
										country_id_ls=self.country_id_ls, tax_id_ls=self.tax_id_ls,\
										local_realigned=self.local_realigned, outdated_index=self.alignment_outdated_index,\
										mask_genotype_method_id=self.mask_genotype_method_id,\
										completedAlignment=self.completedAlignment, \
										completeAlignmentCheckFunction=self.isThisInputAlignmentComplete,\
										reduce_reads=self.reduce_reads)
		alignmentLs = db.filterAlignments(data_dir=self.data_dir, alignmentLs=alignmentLs, min_coverage=self.sequence_min_coverage,\
						max_coverage=self.sequence_max_coverage, sequence_filtered=self.sequence_filtered, \
						individual_site_id_set=set(self.site_id_ls),\
						mask_genotype_method_id=self.mask_genotype_method_id, parent_individual_alignment_id=None,\
						country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls),\
						excludeContaminant=self.excludeContaminant,local_realigned=self.local_realigned,\
						completedAlignment=self.completedAlignment,\
						completeAlignmentCheckFunction=self.isThisInputAlignmentComplete,\
						reduce_reads=self.reduce_reads)
		return alignmentLs
	
	def registerAlignmentAndItsIndexFile(self, workflow=None, alignmentLs=None, data_dir=None, checkFileExistence=True):
		"""
		2013.04.03 cosmetic
		2012.9.18 copied from AlignmentToCallPipeline.py
		2012.6.12
			add argument checkFileExistence
		2012.1.9
			register the input alignments and return in a data structure usd by several other functions
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Registering %s alignments ..."%(len(alignmentLs)))
		alignmentDataLs = []
		for alignment in alignmentLs:
			inputFname = os.path.join(data_dir, alignment.path)
			inputFile = File(alignment.path)	#relative path, induces symlinking or stage-in
			baiFilepath = '%s.bai'%(inputFname)
			if checkFileExistence and (not os.path.isfile(inputFname) or not os.path.isfile(baiFilepath)):
				if not os.path.isfile(inputFname):
					sys.stderr.write("registerAlignmentAndItsIndexFile() Warning: %s does not exist. skip entire alignment.\n"%(inputFname))
				if not os.path.isfile(baiFilepath):
					sys.stderr.write("registerAlignmentAndItsIndexFile() Warning: %s does not exist. skip entire alignment.\n"%(baiFilepath))
				continue
			inputFile.addPFN(PFN("file://" + inputFname, workflow.input_site_handler))
			workflow.addFile(inputFile)
			baiF = File('%s.bai'%alignment.path)
			baiF.addPFN(PFN("file://" + baiFilepath, workflow.input_site_handler))
			workflow.addFile(baiF)
			alignmentData = PassingData(alignment=alignment, jobLs = [], bamF=inputFile, baiF=baiF, file=inputFile,\
							fileLs = [inputFile, baiF])
			alignmentDataLs.append(alignmentData)
		sys.stderr.write("Done.\n")
		return alignmentDataLs
	
	def registerJars(self, workflow=None, ):
		"""
		2011-11-22
			register jars to be used in the worflow
		"""
		parentClass.registerJars(self)
		#2013.06.23
		self.registerOneJar(name="Beagle4Jar", path=os.path.expanduser('~/bin/Beagle/beagle4.jar'))
		#2013.06.13
		self.registerOneJar(name="BeagleJar", path=os.path.expanduser('~/bin/Beagle/beagle.jar'))
		
		#add the MergeSamFiles.jar file into workflow
		self.registerOneJar(name="MergeSamFilesJar", path=os.path.join(self.picard_path, 'MergeSamFiles.jar'))
		self.registerOneJar(name="BuildBamIndexJar", path=os.path.join(self.picard_path, 'BuildBamIndex.jar'))
		self.registerOneJar(name="SamToFastqJar", path=os.path.join(self.picard_path, 'SamToFastq.jar'))	#2013.04.03
		
		self.registerOneJar(name="CreateSequenceDictionaryJar", path=os.path.join(self.picard_path, 'CreateSequenceDictionary.jar'))
		self.registerOneJar(name="AddOrReplaceReadGroupsAndCleanSQHeaderJar", path=os.path.join(self.picard_path, 'AddOrReplaceReadGroupsAndCleanSQHeader.jar'))
		self.registerOneJar(name="AddOrReplaceReadGroupsJar", path=os.path.join(self.picard_path, 'AddOrReplaceReadGroups.jar'))
		self.registerOneJar(name="MarkDuplicatesJar", path=os.path.join(self.picard_path, 'MarkDuplicates.jar'))
		self.registerOneJar(name="SplitReadFileJar", path=os.path.join(self.picard_path, 'SplitReadFile.jar'))
		
		self.registerOneJar(name="GenomeAnalysisTKJar", path=os.path.join(self.gatk_path, 'GenomeAnalysisTK.jar'))
		
		#have to be a different folder, otherwise name clash with the old gatk jar file
		self.registerOneJar(name="GenomeAnalysisTK2Jar", path=os.path.join(self.gatk2_path, 'GenomeAnalysisTK.jar'),\
						folderName='gatk2Jar')
		self.registerOneJar(name="SortSamJar", path=os.path.join(self.picard_path, 'SortSam.jar'))
		self.registerOneJar(name="SamFormatConverterJar", path=os.path.join(self.picard_path, 'SamFormatConverter.jar'))
		
	def registerCustomJars(self, workflow=None, ):
		"""
		2012.1.9
		"""
		parentClass.registerCustomJars(self, workflow=workflow)

	def registerExecutables(self, workflow=None):
		"""
		2012.8.7 remove noDefaultClustersSizeExecutableList, use executableClusterSizeMultiplierList instead
		2012.1.9 a symlink to registerCommonExecutables()
		"""
		parentClass.registerExecutables(self, workflow=workflow)
		
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		#2014.01.08
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'polymorphism/qc/mapper/FilterLocusBasedOnLocusStatFile.py'), \
									name='FilterLocusBasedOnLocusStatFile', \
									clusterSizeMultipler=0.5)
		#2013.10.07
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.expanduser("~/script/freebayes/bin/freebayes"), \
											name='freebayes', clusterSizeMultipler=0.5)
		#2013.10.2
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='CombineBeagleAndPreBeagleVariantsJava', \
											clusterSizeMultipler=0.6)
		#2013.10.13
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "shell/ligateVcf.sh"), \
									name="ligateVcf", clusterSizeMultipler=1)
		#2013.09.17 updated 
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "polymorphism/qc/CheckTwoVCFOverlapCC"), \
											name='CheckTwoVCFOverlapCC', clusterSizeMultipler=1)
		
		#2013.10.07 does not exist anymore
		#self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "mapper/SelectAndSplitAlignment.py"), \
		#									name='selectAndSplit', clusterSizeMultipler=0)
		
		selectAndSplitFasta = Executable(namespace=namespace, name="SelectAndSplitFastaRecords", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		selectAndSplitFasta.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mapper/SelectAndSplitFastaRecords.py"), site_handler))
		executableClusterSizeMultiplierList.append((selectAndSplitFasta, 0))
		
		BuildBamIndexFilesJava = Executable(namespace=namespace, name="BuildBamIndexFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		BuildBamIndexFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((BuildBamIndexFilesJava, 0.5))
		
		#2012.9.21 same as BuildBamIndexFilesJava, but no clustering
		IndexMergedBamIndexJava =  Executable(namespace=namespace, name="IndexMergedBamIndexJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		IndexMergedBamIndexJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((IndexMergedBamIndexJava, 0))
		
		CreateSequenceDictionaryJava = Executable(namespace=namespace, name="CreateSequenceDictionaryJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		CreateSequenceDictionaryJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((CreateSequenceDictionaryJava, 0))
		
		DOCWalkerJava = Executable(namespace=namespace, name="DOCWalkerJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		#no clusters_size for this because it could run on a whole bam for hours
		#DOCWalkerJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		DOCWalkerJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((DOCWalkerJava, 0.05))
		
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
		executableClusterSizeMultiplierList.append((SelectVariantsJava, 0.5))
		
		#2013.09.04
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='CombineVariantsJava', clusterSizeMultipler=0.3)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='CombineVariantsJavaInReduce', \
															clusterSizeMultipler=0.001)
		
		
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
		
		#2013.06.28 use this function 
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "shell/bgzip_tabix.sh"), \
											name='bgzip_tabix', clusterSizeMultipler=4)
		#bgzip_tabix_in_reduce is used in reduce() functions, on whole-scaffold/chromosome VCFs, less clustering
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "shell/bgzip_tabix.sh"), \
											name='bgzip_tabix_in_reduce', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "shell/vcf_convert.sh"), \
											name='vcf_convert', clusterSizeMultipler=1)
		#vcf_convert_in_reduce is used in reduce() functions, on whole-scaffold/chromosome VCFs,
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "shell/vcf_convert.sh"), \
											name='vcf_convert_in_reduce', clusterSizeMultipler=0.2)
		
		vcf_isec = Executable(namespace=namespace, name="vcf_isec", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_isec.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_isec.sh"), site_handler))
		executableClusterSizeMultiplierList.append((vcf_isec, 1))
		
		vcf_concat = Executable(namespace=namespace, name="vcf_concat", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_concat.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		executableClusterSizeMultiplierList.append((vcf_concat, 1))
		
		#vcfSubsetPath is first argument to vcfSubset
		vcfSubset = Executable(namespace=namespace, name="vcfSubset", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcfSubset.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcfSubset.sh"), site_handler))
		vcfSubset.vcfSubsetPath = self.vcfSubsetPath
		executableClusterSizeMultiplierList.append((vcfSubset, 1))
		
		concatGATK = Executable(namespace=namespace, name="concatGATK", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		concatGATK.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		executableClusterSizeMultiplierList.append((concatGATK, 1))
		
		concatSamtools = Executable(namespace=namespace, name="concatSamtools", version=version, \
									os=operatingSystem, arch=architecture, installed=True)
		concatSamtools.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcf_concat.sh"), site_handler))
		executableClusterSizeMultiplierList.append((concatSamtools, 1))
		
		
		#2011.12.21 moved from FilterVCFPipeline.py
		FilterVCFByDepthJava = Executable(namespace=namespace, name="FilterVCFByDepthJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		FilterVCFByDepthJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((FilterVCFByDepthJava, 1 ))
		
		
		#2012.3.1
		MergeFiles = Executable(namespace=namespace, name="MergeFiles", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		MergeFiles.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "shell/MergeFiles.sh"), site_handler))
		executableClusterSizeMultiplierList.append((MergeFiles, 0 ))
		
		#2013.09.17 updated 
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "mapper/CheckTwoVCFOverlap.py"), \
											name='CheckTwoVCFOverlap', clusterSizeMultipler=1)
		
		#2012.9.6
		AppendInfo2SmartPCAOutput = Executable(namespace=namespace, name="AppendInfo2SmartPCAOutput", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		AppendInfo2SmartPCAOutput.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/AppendInfo2SmartPCAOutput.py"), site_handler))
		executableClusterSizeMultiplierList.append((AppendInfo2SmartPCAOutput, 0))
		
				
		AddAlignmentFile2DB = Executable(namespace=namespace, name="AddAlignmentFile2DB", version=version, os=operatingSystem,\
										arch=architecture, installed=True)
		AddAlignmentFile2DB.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "db/input/AddAlignmentFile2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((AddAlignmentFile2DB, 0))
		
		
		MergeSamFilesJava = Executable(namespace=namespace, name="MergeSamFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		MergeSamFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((MergeSamFilesJava, 0))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		#2013.07.09 in order to run vcfsorter.pl from http://code.google.com/p/vcfsorter/
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'shell/pipeCommandOutput2File.sh'), name='vcfsorterShellPipe', clusterSizeMultipler=1)
		#2013.11.22 generic pipeCommandOutput2File
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'shell/pipeCommandOutput2File.sh'), name='pipeCommandOutput2File', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='GATKJava', clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.samtools_path, name='samtools', clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='genotyperJava', clusterSizeMultipler=0.1)
		
		#clustering is controlled by a separate parameter
		#genotyperJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		
		#2013.07.10
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/modifier/AddMissingInfoDescriptionToVCFHeader.py"), \
									name='AddMissingInfoDescriptionToVCFHeader', clusterSizeMultipler=1)
		
		#2013.06.21
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/splitter/SplitVCFFile.py"), \
									name='SplitVCFFile', clusterSizeMultipler=0.01)
		#2012.7.25
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='MergeVCFReplicateHaplotypesJava', \
											clusterSizeMultipler=0.5)
		
		#2013.06.13
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='BeagleJava', clusterSizeMultipler=0.3)
		#2013.06.12 use this simple function to register vcftoolsWrapper
		#vcftoolsPath is first argument to vcftoolsWrapper
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "shell/vcftoolsWrapper.sh"), \
													name='vcftoolsWrapper', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='SortSamFilesJava', clusterSizeMultipler=1)
		#2013.06.06
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='PrintReadsJava', clusterSizeMultipler=1)
		#2013.04.09 
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='addOrReplaceReadGroupsJava', clusterSizeMultipler=0.5)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='AddOrReplaceReadGroupsJava', clusterSizeMultipler=0.5)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.platypus_path, name='Platypus', clusterSizeMultipler=0)
		
	
	bwaIndexFileSuffixLs = ['amb', 'ann', 'bwt', 'pac', 'sa']
	#, 'nhr', 'nin', 'nsq' are formatdb (blast) output, 2012.10.18 i think
	
	def registerBWAIndexFile(self, refFastaFname=None, input_site_handler=None, folderName=""):
		"""
		2012.10.10
		"""
		if input_site_handler is None:
			input_site_handler = self.input_site_handler
		return yh_pegasus.registerRefFastaFile(workflow=self, refFastaFname=refFastaFname, registerAffiliateFiles=True, \
									input_site_handler=input_site_handler,\
									checkAffiliateFileExistence=True, addPicardDictFile=False, \
									affiliateFilenameSuffixLs=self.bwaIndexFileSuffixLs,\
									folderName=folderName)
	
	def addBAMIndexJob(self, workflow=None, BuildBamIndexFilesJava=None, BuildBamIndexJar=None, \
					inputBamF=None,\
					extraArguments=None, parentJobLs=None, extraDependentInputLs=None, \
					transferOutput=True, job_max_memory=None, javaMaxMemory=2500,\
					walltime=60, **keywords):
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
		if javaMaxMemory is not None and job_max_memory is None:
			job_max_memory = javaMaxMemory
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementData.memRequirement
		javaMemRequirement = memRequirementData.memRequirementInStr
		baiFile = File('%s.bai'%inputBamF.name)		
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', BuildBamIndexJar,\
							"VALIDATION_STRINGENCY=LENIENT", \
							"INPUT=", inputBamF, "OUTPUT=", baiFile]
					#not including 'SORT_ORDER=coordinate'
					#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexJar won't fail.)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputBamF, BuildBamIndexJar])
		
		job= self.addGenericJob(executable=BuildBamIndexFilesJava, inputFile=None,\
							outputFile=None, outputArgumentOption="-o", \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraOutputLs=[baiFile],\
							transferOutput=transferOutput, \
							extraArgumentList=extraArgumentList, \
							job_max_memory=memRequirementData.memRequirement, walltime=walltime, **keywords)
		job.bamFile = inputBamF
		job.baiFile = baiFile
		#job.parentJobLs is where the actual alignment job and its bam/sam output are.
		return job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.1.9
			abstract function
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
	
	
	def addRefFastaFaiIndexJob(self, workflow=None, samtools=None, refFastaF=None, \
							parentJobLs=None, transferOutput=True, job_max_memory=1000,\
							walltime=300, **keywords):
		"""
		2013.3.23 use addGenericJob()
		2011-11-25
			the *.fai file of refFastaF is required for GATK
		"""
		sys.stderr.write("\t Adding SAMtools index fasta job ...")
		refFastaIndexFname = '%s.fai'%(refFastaF.name)	# the .fai file is required for GATK
		refFastaIndexF = File(refFastaIndexFname)
		frontArgumentList = ["faidx", refFastaF]
		fastaIndexJob = self.addGenericJob(executable=samtools, inputFile=None,\
							outputFile=None, outputArgumentOption=None, \
							parentJobLs=parentJobLs, extraDependentInputLs=[refFastaF], \
							extraOutputLs=[refFastaIndexF],\
							frontArgumentList=frontArgumentList,\
							transferOutput=transferOutput, \
							extraArgumentList=None, \
							job_max_memory=job_max_memory, walltime=walltime, **keywords)
		
		
		fastaIndexJob.refFastaIndexF = refFastaIndexF
		sys.stderr.write("\n")
		return fastaIndexJob
	
	def addRefFastaDictJob(self, workflow=None, CreateSequenceDictionaryJava=None, CreateSequenceDictionaryJar=None, \
						refFastaF=None, parentJobLs=None, transferOutput=True, job_max_memory=1000,\
						walltime=120, **keywords):
		"""
		2013.3.23 use addGenericJob()
		2012.9.14 bugfix. add argument CreateSequenceDictionaryJar
		2011-11-25
			# the .dict file is required for GATK
		"""
		sys.stderr.write("\t Adding picard CreateSequenceDictionaryJar job ...")
		refFastaDictFname = '%s.dict'%(os.path.splitext(refFastaF.name)[0])
		refFastaDictF = File(refFastaDictFname)
		if CreateSequenceDictionaryJar is None:
			sys.stderr.write("Error: CreateSequenceDictionaryJar is nothing (%s) inside addRefFastaDictJob(). Quit.\n"%\
							(CreateSequenceDictionaryJar))
			raise
		fastaDictJob = self.addGenericJavaJob(executable=CreateSequenceDictionaryJava, jarFile=CreateSequenceDictionaryJar, \
					inputFile=refFastaF, inputArgumentOption="REFERENCE=", \
					inputFileList=None, argumentForEachFileInInputFileList=None,\
					outputFile=refFastaDictF, outputArgumentOption="OUTPUT=",\
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					frontArgumentList=[], extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=None, no_of_cpus=None, walltime=walltime, sshDBTunnel=None, **keywords)
		
		fastaDictJob.refFastaDictF = refFastaDictF
		sys.stderr.write("\n")
		return fastaDictJob
	
	def addRefFastaJobDependency(self, workflow=None, job=None, refFastaF=None, fastaDictJob=None, refFastaDictF=None, \
								fastaIndexJob = None, refFastaIndexF = None, **keywords):
		"""
		Examples:
			self.addRefFastaJobDependency(job=job, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
					
		2013.3.23 refFastaDictF and refFastaIndexF are optional if fastaDictJob & fastaIndexJob are for sure not None.
			also call self.addJobUse()
		2011-9-14
		"""
		if workflow is None:
			workflow = self
		if fastaIndexJob:	#2011-7-22 if job doesn't exist, don't add it. means this job isn't necessary to run.
			self.depends(parent=fastaIndexJob, child=job)
			if refFastaIndexF is None:
				refFastaIndexF = fastaIndexJob.output
			self.addJobUse(job=job, file=refFastaIndexF, transfer=True, register=True, link=Link.INPUT)
		elif refFastaIndexF:
			self.addJobUse(job=job, file=refFastaIndexF, transfer=True, register=True, link=Link.INPUT)
			
		if fastaDictJob:
			self.depends(parent=fastaDictJob, child=job)
			if refFastaDictF is None:
				refFastaDictF = fastaDictJob.output
			self.addJobUse(job=job, file=refFastaDictF, transfer=True, register=True, link=Link.INPUT)
		elif refFastaDictF:
			self.addJobUse(job=job, file=refFastaDictF, transfer=True, register=True, link=Link.INPUT)
			
		if fastaIndexJob or fastaDictJob:
			self.addJobUse(job=job, file=refFastaF, transfer=True, register=True, link=Link.INPUT)
	
	def addVCFFormatConvertJob(self, workflow=None, vcf_convert=None, inputF=None, outputF=None, \
							parentJob=None, parentJobLs=None, \
							extraDependentInputLs=None, transferOutput=False, job_max_memory=1000,\
							walltime=120, **keywords):
		"""
		2013.06.14
			use addGenericJob()
		2011-11-4
		"""
		if workflow is None:
			workflow = self
		if vcf_convert is None:
			vcf_convert = self.vcf_convert
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = []
		extraOutputLs = []
		key2ObjectForJob = {}
		
		job = self.addGenericJob(executable=vcf_convert, inputFile=inputF, inputArgumentOption="", \
					outputFile=outputF, outputArgumentOption="", \
					parentJob=parentJob, \
					parentJobLs=parentJobLs, \
					extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, \
					sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, \
					job_max_memory=job_max_memory, walltime=walltime,\
					**keywords)
		return job
	
	def addBeagle3Job(self, workflow=None, executable=None, BeagleJar=None, \
					phasedBeagleInputFile=None,\
					likelihoodBeagleInputFile=None, triosBeagleInputFile=None, pairsBeagleInputFile=None,\
					unphasedBeagleInputFile=None,\
					markersBeagleInputFile=None,
					outputFnamePrefix=None, noOfIterations=None, noOfSamplingHaplotypesPerSample=None, \
					parentJobLs=None, transferOutput=True, job_max_memory=2000,\
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=None, no_of_cpus=None, walltime=120, **keywords):
		"""
		i.e.
			
		2013.06.13 a generic function to add Beagle version 3 jobs
		
		crocea@vervetNFS:~/bin/Beagle/example/imputation$ 
		prefix=method_36_Contig791_replicated_phasedByGATK_noMendelError_unphased_minProb0.65;
		java -Xmx15g -XX:PermSize=1000m -XX:MaxPermSize=6000m
			-jar ~/bin/Beagle/beagle.jar like=$prefix\_familySize1.bgl
			pairs=$prefix\_familySize2.bgl trios=$prefix\_familySize3.bgl
			markers=$prefix.markers out=$prefix\_beagled missing=?
			unphased=...bgl
			
		Note:
			"missing=? " can not be used if all input is likelihood data.
		
		Beagle help http://faculty.washington.edu/browning/beagle/beagle.html
		
		niterations=<number of iterations> where <number of iterations> is a positive even 
			integer giving the number of iterations of the phasing algorithm. If an odd integer is 
			specified, the next even integer is used. The niterations argument is optional. The default 
			value is niterations=10. The default value typically gives good accuracy
		
		nsamples=<number of samples> where <number of samples> is positive integer giving 
			the number of haplotype pairs to sample for each individual during each iteration of the 
			phasing algorithm. The nsamples argument is optional. The default value is nsamples=4.
			If you are phasing an extremely large sample (say > 4000 individuals), you may want to 
			use a smaller nsamples parameter (e.g. 1 or 2) to reduce computation time. If you are 
			phasing a small sample (say < 200 individuals), you may want to use a larger nsamples
			parameter (say 10 or 20) to increase accuracy.
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.java
		if BeagleJar is None:
			BeagleJar = self.BeagleJar
		if frontArgumentList is None:
			frontArgumentList = []
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if extraOutputLs is None:
			extraOutputLs = []
		
		#place holder for key output files, corresponding to these input like=, unphased=, trios=, pairs=
		key2ObjectForJob = {"likelihoodPhasedOutputFile":None,\
						"singletonPhasedOutputFile":None,\
						"pairsPhasedOutputFile":None,\
						"triosPhasedOutputFile":None,\
						}
		
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		
		frontArgumentList.extend([javaMemRequirement, '-jar', BeagleJar, "out=%s"%(outputFnamePrefix)])
		extraDependentInputLs.append(BeagleJar)
		
		if noOfIterations is not None:
			frontArgumentList.append("niterations=%s"%(noOfIterations))
		if noOfSamplingHaplotypesPerSample is not None:
			frontArgumentList.append("nsamples=%s"%(noOfSamplingHaplotypesPerSample))
		
		logFile=File('%s.log'%(outputFnamePrefix))
		extraOutputLs.append(logFile)
		key2ObjectForJob['logFile'] = logFile
		
		phasedBeagleFileSuffix = '.bgl.phased.gz'
		onlyLikelihoodInput = True	#determines whether missing=? should be added
		if likelihoodBeagleInputFile:
			frontArgumentList.extend(["like=%s"%(likelihoodBeagleInputFile.name)])
			extraDependentInputLs.append(likelihoodBeagleInputFile)
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(likelihoodBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['likelihoodPhasedOutputFile'] = phasedFile
			#additional output for the likelihood input
			for suffix in ['.bgl.dose.gz', '.bgl.gprobs.gz', '.bgl.r2']:
				_outputFile=File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, suffix))
				extraOutputLs.append(_outputFile)
		
		if triosBeagleInputFile:
			frontArgumentList.extend(["trios=%s"%(triosBeagleInputFile.name)])
			extraDependentInputLs.append(triosBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(triosBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['triosPhasedOutputFile'] = phasedFile
		
		if pairsBeagleInputFile:
			frontArgumentList.extend(["pairs=%s"%(pairsBeagleInputFile.name)])
			extraDependentInputLs.append(pairsBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(pairsBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['pairsPhasedOutputFile'] = phasedFile
		
		if unphasedBeagleInputFile:
			frontArgumentList.extend(["unphased=%s"%(unphasedBeagleInputFile.name)])
			extraDependentInputLs.append(unphasedBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(unphasedBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['singletonPhasedOutputFile'] = phasedFile
			#additional likelihood-related output
			for suffix in ['.bgl.dose.gz', '.bgl.gprobs.gz', '.bgl.r2']:
				_outputFile=File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, suffix))
				extraOutputLs.append(_outputFile)
		
		if phasedBeagleInputFile:	#reference panel
			frontArgumentList.extend(["phased=%s"%(phasedBeagleInputFile.name)])
			extraDependentInputLs.append(phasedBeagleInputFile)
			
		
		if markersBeagleInputFile:
			frontArgumentList.extend(["markers=%s"%(markersBeagleInputFile.name)])
			extraDependentInputLs.append(markersBeagleInputFile)
		if not onlyLikelihoodInput:
			frontArgumentList.append("missing=?")
		
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=None, \
					inputArgumentOption=None,  inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					frontArgumentList=frontArgumentList, extraArguments=extraArguments, \
					extraArgumentList=extraArgumentList, job_max_memory=job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		
		return job
	
	def addBeagle4Job(self, workflow=None, executable=None, BeagleJar=None, \
					inputFile=None, refPanelFile=None, pedFile=None,\
					outputFnamePrefix=None, \
					burninIterations=None, phaseIterations=None, \
					noOfSamplingHaplotypesPerSample=None, singlescale=None, duoscale=None, trioscale=None,\
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=None, parentJobLs=None, transferOutput=True, job_max_memory=2000,\
					no_of_cpus=None, walltime=120, **keywords):
		"""
		i.e.
			
		2013.06.22 a generic function to add Beagle version 4 job
		
		java commandline example:
			...
		
		Beagle help http://faculty.washington.edu/browning/beagle/beagle.html
		
		gt=[file] specifies a VCF file containing a GT (genotype) format field for each marker.
		gl=[file] specifies a VCF file containing a GL or PL (genotype likelihood) format field for
			each marker. If both GL and PL format fields are present for a sample, the GL format will
			be used. See also the maxlr parameter.
		ref=[file] specifies a reference VCF file containing additional samples and phased
			genotypes for each marker. Use of an appropriate reference panel can increase accuracy.
		ped=[file] specifies a Linkage-format pedigree file for specifying family relationships.
			The pedigree file has one line per individual. The first 4 white-space delimited fields of
			each line are 1) pedigree ID, 2) individual ID, 3) father's ID, and 4) mother's ID. A "0" is
			used in column 3 or 4 if the father or mother is unknown. Beagle uses the data in columns
			2-4 to identify parent-offspring duos and trios. Any or all columns of the pedigree file
			after column 4 may be omitted. See also the duoscale and trioscale parameters.
		out=[prefix] specifies the output filename prefix. The prefix may be an absolute or
			relative filename, but it cannot be a directory name.
		excludesamples=[file] specifies a file containing non-reference samples (one sample per
			line) to be excluded from the analysis and output files.
		excluderefsamples=[file] specifies a file containing reference samples (one sample per
			line) to be excluded from the analysis.
		excludemarkers=[file] specifies a file containing markers (one marker per line) to be
			excluded from the analysis and the output files. An excluded marker identifier can either
			be an identifier from the VCF record's ID field or genomic coordinates in the format:
			CHROM:POS.
		chrom=[chrom:start-end] specifies a chromosome or chromosome interval using a
			chromosome identifier in the VCF file and the starting and ending positions of the
			interval. The entire chromosome, the beginning of the chromosome, and the end of a
			chromosome can be specified by chrom=[chrom], chrom=[chrom:-end], and chrom=[chrom:start-] respectively.
		
		window=[positive integer] (default: window=24000).
		overlap=[positive integer] (default: overlap=3000)
		gprobs=[true/false] (default: gprobs=true).
		usephase=[true/false] (default: usephase=false).
		singlescale=[positive number] (default: singlescale=1.0), change the scale to x, program samples x*x faster
		duoscale=[positive number]
		trioscale=[positive number] 
		burnin-its=[non-negative integer] (default: burnin=5).
		phase-its=[non-negative integer] (default: phase-its=5)
		
		
		Advanced options not recommended for general use:
		
		dump=[file] specifies a file containing sample identifiers (one identifier per line). For
			each marker window, all the sampled haplotypes for these individuals which are sampled
			after the burn-in iterations are printed to an output VCF files (dump.[window #].gz).
		nsamples=[positive integer] specifies the number of haplotype pairs to sample for each
			individual during each iteration of the algorithm (default: nsamples=4).
		buildwindow=[positive integer] specifies the number of markers used to build the
			haplotype frequency model at each locus (default: buildwindow=500).
		
		
		Three output files are created whose names begin with the output file prefix specified on
			the command line argument and whose names end with the suffixes: .log, .vcf.gz, and .ibd.

		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.java
		if BeagleJar is None:
			BeagleJar = self.Beagle4Jar
		if frontArgumentList is None:
			frontArgumentList = []
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if extraOutputLs is None:
			extraOutputLs = []
		
		key2ObjectForJob = {}
		
		extraArgumentList.extend(["out=%s"%(outputFnamePrefix)])
		if inputFile:
			extraArgumentList.append("gl=%s"%(inputFile.name))
			extraDependentInputLs.append(inputFile)
		if refPanelFile:
			extraArgumentList.append("ref=%s"%(refPanelFile.name))
			extraDependentInputLs.append(refPanelFile)
		if pedFile:
			extraArgumentList.append("ped=%s"%(pedFile.name))
			extraDependentInputLs.append(pedFile)
		if burninIterations is not None:
			extraArgumentList.append("burnin-its=%s"%(burninIterations))
		if phaseIterations is not None:
			extraArgumentList.append("phase-its=%s"%(phaseIterations))
		if noOfSamplingHaplotypesPerSample is not None:
			extraArgumentList.append("nsamples=%s"%(noOfSamplingHaplotypesPerSample))
		if singlescale is not None:
			extraArgumentList.append("singlescale=%s"%(singlescale))
		if duoscale is not None:
			extraArgumentList.append("duoscale=%s"%(duoscale))
		if trioscale is not None:
			extraArgumentList.append("trioscale=%s"%(trioscale))
		
		outputVCFFile = File("%s.vcf.gz"%(outputFnamePrefix))
		extraOutputLs.append(outputVCFFile)	#this would be accessible through job.output and job.vcfOutputFile 
		key2ObjectForJob['vcfOutputFile'] = outputVCFFile
		
		logFile=File('%s.log'%(outputFnamePrefix))
		extraOutputLs.append(logFile)
		key2ObjectForJob['logFile'] = logFile	#this would be accessible as job.logFile
		
		job = self.addGenericJavaJob(executable=executable, jarFile=BeagleJar, \
					inputFile=None, inputArgumentOption=None, \
					inputFileList=None, argumentForEachFileInInputFileList=None,\
					outputFile=None, outputArgumentOption=None,\
					frontArgumentList=frontArgumentList, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					extraOutputLs=extraOutputLs, \
					extraDependentInputLs=extraDependentInputLs, \
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					key2ObjectForJob=key2ObjectForJob, \
					no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		return job
	
	def convertVCF2Beagle(self, workflow=None, VCFJobData=None, outputDirJob=None, \
						outputFileBasenamePrefix=None,\
						outputPedigreeJob=None, pedigreeSplitStructure=None,\
						transferOutput=False, \
						job_max_memory=None, walltime=None):
		"""
		2013.06.23 deprecated as Beagle v4 is used and it accepts VCF files
		2013.05.01
			convert VCF file into Beagle input format
		"""
		
		#replicate the individuals involved in more than one trio/duo
		#2012.4.2 replicate individuals who appear in more than 1 families
		#note: remove individuals in pedigree file but not in VCF file
		
		#GATK ProduceBeagleInput
		beagleLikelihoodFile = File(os.path.join(outputDirJob.folder, "%s.bgl"%(outputFileBasenamePrefix)))
		produceBeagleInputJob = self.addGATKJob(executable=self.ProduceBeagleInputJava, \
					GATKAnalysisType="ProduceBeagleInput", \
					inputFile=VCFJobData.file, inputArgumentOption="--variant:VCF", \
					refFastaFList=self.registerReferenceData.refFastaFList, \
					outputFile=beagleLikelihoodFile, \
					parentJobLs=VCFJobData.jobLs + [outputDirJob], \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=None, \
					extraOutputLs=None, extraDependentInputLs=None, \
					job_max_memory=job_max_memory, \
					no_of_cpus=None, walltime=walltime)
		
		#a SplitPedigreeVCFIntoBeagleTriosDuosFiles job
		#need pedigreeFile (from outputPedigreeJob)
		outputFnamePrefix = os.path.join(outputDirJob.folder, outputFileBasenamePrefix)
		
		key2File= {'size1File': None, 'size2File': None, 'size3File':None}	#to store output files
			#size1File: output for singletons (likelihood format)
			#size2File: output for duos (Beagle genotype format)
			#size3File: output for trios (Beagle genotype format)
		markerFile = File("%s.markers"%(outputFnamePrefix))
		extraOutputLs = []	#first is the markers file
		extraOutputLs.append(markerFile)
		key2File['markerFile'] = markerFile
		for familySize, familyLs in pedigreeSplitStructure.familySize2familyLs.iteritems():
			if familyLs:	#non-empty
				outputFile = File('%s_familySize%s.bgl'%(outputFnamePrefix, familySize))
				key2File['size%sFile'%(familySize)] = outputFile
				extraOutputLs.append(outputFile)
		
		splitPedigreeVCFIntoBeagleTriosDuosFilesJob = self.addGenericJob(executable=self.SplitPedigreeVCFIntoBeagleTriosDuosFiles, \
						inputFile=VCFJobData.file, \
						outputFile=None, \
						parentJobLs=VCFJobData.jobLs + [produceBeagleInputJob, outputDirJob, outputPedigreeJob], \
						extraDependentInputLs=[produceBeagleInputJob.output, outputPedigreeJob.output], \
						extraOutputLs=extraOutputLs,  transferOutput=transferOutput, \
						extraArguments=None, \
						extraArgumentList=["--gatkPrintBeagleFname", produceBeagleInputJob.output, \
										"--plinkPedigreeFname", outputPedigreeJob.output, \
										"--minProbForValidCall %s"%(self.minProbForValidCall), \
										"--dummyIndividualNamePrefix dummy", \
										"--outputFnamePrefix", outputFnamePrefix], \
						job_max_memory=job_max_memory, \
						sshDBTunnel=None, key2ObjectForJob=key2File, objectWithDBArguments=None, \
						no_of_cpus=None, walltime=walltime)
		
		#
		return splitPedigreeVCFIntoBeagleTriosDuosFilesJob
	
	def addBGZIP_tabix_Job(self, workflow=None, bgzip_tabix=None, parentJob=None, inputF=None, outputF=None, \
							transferOutput=False, parentJobLs=None, tabixArguments=None,\
							extraDependentInputLs=None, job_max_memory=2000, **keywords):
		"""
		Doc:
			access the tbi file via job.tbi_F
		2013.06.14 added argument job_max_memory
		2012.8.17 if transferOutput is None, do not register output files as OUTPUT with transfer flag
			use addGenericJob()
		2011.12.20
			pass additional tabix arguments to bgzip_tabix shell script
		2011-11-4
		
		"""
		if workflow is None:
			workflow = self
		if bgzip_tabix is None:
			bgzip_tabix = self.bgzip_tabix
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
		extraDependentInputLs.append(self.tabixExecutableFile)
		extraDependentInputLs.append(self.bgzipExecutableFile)
		
		job = self.addGenericJob(executable=bgzip_tabix, inputFile=inputF, inputArgumentOption="", \
					outputFile=outputF, outputArgumentOption="", \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=tabixArguments, extraArgumentList=extraArgumentList, \
					job_max_memory=job_max_memory,  sshDBTunnel=None, \
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
	
	def addGATKCombineVariantsJob(self, workflow=None, executable=None, GenomeAnalysisTKJar=None, \
							refFastaFList=None, inputFileList=None, argumentForEachFileInInputFileList="--variant",\
							outputFile=None, genotypeMergeOptions='UNSORTED', \
					parentJobLs=None, transferOutput=True, job_max_memory=2000,walltime=None,\
					extraArguments=None, extraArgumentList=None, extraDependentInputLs=None, **keywords):
		"""
		examples
			add input file to this job via
				gatkUnionJob = self.addGATKCombineVariantsJob(...)
				self.addInputToStatMergeJob(statMergeJob=gatkUnionJob, parentJobLs=[gatk_job], \
										inputArgumentOption="--variant")
			OR through the inputFileList argument
				gatkUnionJob = self.addGATKCombineVariantsJob(.., inputFileList=[gatk_job.output, ...])
		
			concatVCFFilename = os.path.join(outputDirJob.folder, '%s.vcf'%(passingData.fileBasenamePrefix))
			concatVCFFile = File(concatVCFFilename)
			concatJob = self.addGATKCombineVariantsJob(executable=None, GenomeAnalysisTKJar=None, \
								refFastaFList=None, inputFileList=None, argumentForEachFileInInputFileList="--variant",\
								outputFile=concatVCFFile, genotypeMergeOptions='UNSORTED', \
						parentJobLs=[outputDirJob], transferOutput=False, job_max_memory=job_max_memory,\
						walltime=walltime,\
						extraArguments=None, extraArgumentList=['--assumeIdenticalSamples'], extraDependentInputLs=None)
			
			for intervalJob in intervalJobLs:
				self.addInputToStatMergeJob(statMergeJob=concatJob, inputF=intervalJob.output, inputArgumentOption="--variant",\
								parentJobLs=[intervalJob], extraDependentInputLs=intervalJob.outputLs[1:])
		
		2012.2.26 replacing addVCFConcatJob using GATK's CombineVariants analysis-type
		Value for genotypeMergeOptions:
			UNIQUIFY
				Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
			PRIORITIZE
				Take genotypes in priority order (see the priority argument).
			UNSORTED
				Take the genotypes in any order.
			REQUIRE_UNIQUE
				Require that all samples/genotypes be unique between all inputs.
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.CombineVariantsJava
		if GenomeAnalysisTKJar is None:
			GenomeAnalysisTKJar = self.GenomeAnalysisTK2Jar
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		if genotypeMergeOptions:
			extraArgumentList.append("-genotypeMergeOptions %s"%(genotypeMergeOptions))
		
		
		job = self.addGATKJob(workflow=workflow, executable=executable, GenomeAnalysisTKJar=GenomeAnalysisTKJar, \
							GATKAnalysisType="CombineVariants",\
					inputFile=None, inputArgumentOption=None, refFastaFList=refFastaFList, inputFileList=inputFileList,\
					argumentForEachFileInInputFileList=argumentForEachFileInInputFileList,\
					outputFile=outputFile, \
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					frontArgumentList=None, extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					extraOutputLs=None, \
					extraDependentInputLs=extraDependentInputLs, no_of_cpus=None, walltime=walltime, **keywords)
		return job
	
	def addGATKJob(self, workflow=None, executable=None, GenomeAnalysisTKJar=None, GATKAnalysisType=None,\
					inputFile=None, inputArgumentOption=None, refFastaFList=None, inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					interval=None, outputFile=None, outputArgumentOption="--out", \
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=None, \
					parentJobLs=None, transferOutput=True, \
					no_of_cpus=None, job_max_memory=2000, walltime=120, \
					key2ObjectForJob=None, **keywords):
		"""
		i.e.
			indelUnionOutputF = File(os.path.join(gatkIndelDirJob.folder, '%s.indel.vcf'%chromosome))
			selectIndelJob = self.addGATKJob(executable=self.SelectVariantsJava, GATKAnalysisType="SelectVariants",\
				inputFile=gatkUnionJob.output, inputArgumentOption="--variant", refFastaFList=refFastaFList, inputFileList=None,\
				argumentForEachFileInInputFileList=None,\
				interval=None, outputFile=indelUnionOutputF, \
				parentJobLs=[gatkIndelDirJob, gatkUnionJob], transferOutput=False, job_max_memory=job_max_memory,\
				frontArgumentList=None, extraArguments="-selectType INDEL", \
				extraArgumentList=None, extraOutputLs=None, \
				extraDependentInputLs=None, no_of_cpus=None, walltime=None)
		2013.2.26 add gatkVCFIndexOutput to the extraOutputLs if outputFile is a .vcf file
		2013.2.26 a generic function to add GATK jobs
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.java
		if GenomeAnalysisTKJar is None:
			GenomeAnalysisTKJar = self.GenomeAnalysisTK2Jar
		if frontArgumentList is None:
			frontArgumentList = []
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if extraOutputLs is None:
			extraOutputLs = []
		
		refFastaFile = refFastaFList[0]
		extraDependentInputLs.extend(refFastaFList)
		
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		
		frontArgumentList.extend([javaMemRequirement, '-jar', GenomeAnalysisTKJar, "--reference_sequence",  refFastaFile, \
						"-T", GATKAnalysisType,\
						self.defaultGATKArguments])
		extraDependentInputLs.append(GenomeAnalysisTKJar)
		if no_of_cpus:
			frontArgumentList.append("--num_threads %s"%no_of_cpus)
		if interval is not None:
			if isinstance(interval, File):	#2012.7.30
				extraDependentInputLs.append(interval)
				frontArgumentList.extend(["-L:bed", interval])
			else:
				frontArgumentList.extend(["-L", interval])
		gatkVCFIndexOutput = None
		if outputFile:	#2013.06.09 bugfix
			outputFnameSuffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(outputFile.name)
			if outputFnameSuffix=='.vcf':	#2013.2.26 GATK will generate an index file along with the VCF output
				gatkVCFIndexOutputFname = '%s.idx'%(outputFile.name)
				gatkVCFIndexOutput = File(gatkVCFIndexOutputFname)
				extraOutputLs.append(gatkVCFIndexOutput)
		
		if inputFile:	#2013.05.04 this could be None in some cases (addGATKCombineVariantsJob, inputFile is added in later)
			inputFnameSuffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(inputFile.name)
			if inputFnameSuffix=='.vcf':	#2013.04.16 GATK will generate an index file along with the VCF input
				inputGATKVCFIndexFname = '%s.idx'%(inputFile.name)
				inputGATKVCFIndexFile = File(inputGATKVCFIndexFname)
				extraOutputLs.append(inputGATKVCFIndexFile)
			
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=inputFile, \
					inputArgumentOption=inputArgumentOption,  inputFileList=inputFileList,\
					argumentForEachFileInInputFileList=argumentForEachFileInInputFileList,\
					outputFile=outputFile, outputArgumentOption=outputArgumentOption, \
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					frontArgumentList=frontArgumentList, extraArguments=extraArguments, \
					extraArgumentList=extraArgumentList, job_max_memory=job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		
		job.gatkVCFIndexOutput = gatkVCFIndexOutput
		return job
	
	
	def addGATKRealignerTargetCreatorJob(self, workflow=None, executable=None, GenomeAnalysisTKJar=None, \
							refFastaFList=None, inputFile=None, inputArgumentOption="-I", indelVCFFile=None, \
							outputFile=None, interval=None, \
					parentJobLs=None, transferOutput=True, job_max_memory=2000,walltime=60,\
					extraArguments=None, extraArgumentList=None, extraDependentInputLs=None, **keywords):
		"""
		2013.04.08
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.RealignerTargetCreatorJava
		if GenomeAnalysisTKJar is None:
			GenomeAnalysisTKJar = self.GenomeAnalysisTK2Jar
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		if indelVCFFile:
			extraArgumentList.extend(["-known:vcf", indelVCFFile])
			if indelVCFFile not in extraDependentInputLs:
				extraDependentInputLs.append(indelVCFFile)
		job = self.addGATKJob(executable=executable, GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
						GATKAnalysisType='RealignerTargetCreator',\
						inputFile=inputFile, inputArgumentOption=inputArgumentOption, \
						refFastaFList=refFastaFList, inputFileList=None,\
						argumentForEachFileInInputFileList=None,\
						interval=interval, outputFile=outputFile, \
						parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
						frontArgumentList=None, extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
						extraOutputLs=None, \
						extraDependentInputLs=extraDependentInputLs, no_of_cpus=None, walltime=walltime)
		return job

	def addGATKOutputAlignmentJob(self, workflow=None, executable=None, GenomeAnalysisTKJar=None, GATKAnalysisType=None, \
					refFastaFList=None, inputFile=None, inputArgumentOption="-I", \
					inputFileList=None, argumentForEachFileInInputFileList=None,\
					outputFile=None, interval=None, \
					extraArguments=None, extraArgumentList=None, extraDependentInputLs=None, extraOutputLs=None,\
					parentJobLs=None, transferOutput=True, no_of_cpus=None, job_max_memory=2000,walltime=60,\
					needBAMIndexJob=True, **keywords):
		"""
		2013.06.06 a variant of addGATKJob() that outputs a bam file and thus a bam index file might be needed.
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.GATKJava
		if GenomeAnalysisTKJar is None:
			GenomeAnalysisTKJar = self.GenomeAnalysisTK2Jar
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		job = self.addGATKJob(executable=executable, GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
						GATKAnalysisType=GATKAnalysisType,\
						refFastaFList=refFastaFList, \
						inputFile=inputFile, inputArgumentOption=inputArgumentOption, \
						inputFileList=inputFileList, argumentForEachFileInInputFileList=argumentForEachFileInInputFileList,\
						interval=interval, outputFile=outputFile, \
						parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
						frontArgumentList=None, extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
						extraOutputLs=extraOutputLs, \
						extraDependentInputLs=extraDependentInputLs, no_of_cpus=no_of_cpus, walltime=walltime)
		if needBAMIndexJob:
			# add the index job on the bam file
			bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
						BuildBamIndexJar=self.BuildBamIndexJar, \
						inputBamF=job.output, parentJobLs=[job], \
						transferOutput=transferOutput, job_max_memory=2500, walltime=max(50, walltime/2))
		else:
			bamIndexJob = None
		job.bamIndexJob = bamIndexJob
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
	
	def addVCFBeforeAfterFilterStatJob(self, executable=None, chromosome=None, outputF=None, vcf1=None, vcf2=None,\
									lastVCFJob=None, currentVCFJob=None,\
									statMergeJob=None, statMergeJobLs=None, parentJobLs=None):
		"""
		
		examples:
			outputF = File(os.path.join(self.statDirJob.output, '%s.noOfLociAfterFilterLiftover.tsv'%(intervalFileBasenamePrefix)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
								vcf1=vcfSorterJob.output, currentVCFJob=returnData.filterLiftoverVariantsJob, \
								statMergeJobLs=[self.noOfLociChangeAfterFilterLiftOverMergeJob,\
											self.noOfLociPerContigAfterFilterLiftOverMergeJob], \
								parentJobLs=[vcfSorterJob, returnData.filterLiftoverVariantsJob, self.statDirJob])
								
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									currentVCFJob=currentVCFJob, lastVCFJob=lastVCFJob,\
									statMergeJob=filterByMaxSNPMissingRateMergeJob)
		
		2013.07.11 added argument statMergeJobLs
			statMergeJob could be None.
		2013.06.19 moved from vervet/src/qc/FilterVCFPipeline.py
		2013.06.11 renamed old arguments to lastVCFJob, currentVCFJob
		2012.7.30
		"""
		if vcf1 is None and lastVCFJob:
			vcf1 = lastVCFJob.output
		if vcf2 is None and currentVCFJob:
			vcf2 = currentVCFJob.output
		if parentJobLs is None:
			parentJobLs = []
		if lastVCFJob:
			parentJobLs.append(lastVCFJob)
		if currentVCFJob:
			parentJobLs.append(currentVCFJob)
		if executable is None:
			executable = self.CheckTwoVCFOverlap
		if statMergeJobLs is None:
			statMergeJobLs = []
		vcfFilterStatJob = self.addCheckTwoVCFOverlapJob(executable=executable, \
					vcf1=vcf1, vcf2=vcf2, \
					chromosome=chromosome, chrLength=None, \
					outputFile=outputF, parentJobLs=parentJobLs, \
					extraDependentInputLs=None, transferOutput=False, extraArguments=None, job_max_memory=1000, \
					perSampleMismatchFraction=False)
		if statMergeJob:
			self.addInputToStatMergeJob(statMergeJob=statMergeJob, \
							inputF=vcfFilterStatJob.output , \
							parentJobLs=[vcfFilterStatJob])
		for _statMergeJob in statMergeJobLs:
			self.addInputToStatMergeJob(statMergeJob=_statMergeJob, \
							inputF=vcfFilterStatJob.output , \
							parentJobLs=[vcfFilterStatJob])
		return vcfFilterStatJob
	
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
		if workflow is None:
			workflow = self
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
		for extraInputFile in extraDependentInputLs:
			if extraInputFile:
				job.uses(extraInputFile, transfer=True, register=True, link=Link.INPUT)
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
	
	def addSelectVariantsJob(self, workflow=None, SelectVariantsJava=None, GenomeAnalysisTKJar=None, \
					inputF=None, outputF=None, \
					interval=None,\
					refFastaFList=None, sampleIDKeepFile=None, snpIDKeepFile=None, sampleIDExcludeFile=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=True, \
					extraArguments=None, extraArgumentList=None, job_max_memory=2000, walltime=None,\
					**keywords):
		"""
		if some samples in sampleIDKeepFile or sampleIDExcludeFile are not in inputF, then this option should be passed 
			extraArguments="--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES"
		2013.04.16 updated to use addGATKJob
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
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if SelectVariantsJava is None:
			SelectVariantsJava = self.SelectVariantsJava
		
		#2013.05.02 addGATKJob will take care of these two
		#extraDependentInputLs.append(inputF)
		#extraDependentInputLs.append(GenomeAnalysisTKJar)	#2013.2.14
		#extraDependentInputLs = extraDependentInputLs + refFastaFList
		if sampleIDKeepFile:
			extraArgumentList.extend(['--sample_file', sampleIDKeepFile])
			extraDependentInputLs.append(sampleIDKeepFile)
		if snpIDKeepFile:
			extraArgumentList.extend(['--keepIDs', snpIDKeepFile])
			extraDependentInputLs.append(snpIDKeepFile)
		if sampleIDExcludeFile:
			extraArgumentList.extend(['--exclude_sample_file', sampleIDExcludeFile])
			extraDependentInputLs.append(sampleIDExcludeFile)
		
		job = self.addGATKJob(executable=SelectVariantsJava, GenomeAnalysisTKJar=GenomeAnalysisTKJar, \
				GATKAnalysisType="SelectVariants",\
				inputFile=inputF, inputArgumentOption="--variant:vcf", refFastaFList=refFastaFList, inputFileList=None,\
				argumentForEachFileInInputFileList=None,\
				interval=interval, outputFile=outputF, \
				parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
				frontArgumentList=None, extraArguments=extraArguments, \
				extraArgumentList=extraArgumentList, extraOutputLs=None, \
				extraDependentInputLs=extraDependentInputLs, no_of_cpus=None, walltime=walltime,\
				key2ObjectForJob=None, **keywords)
		return job
	
	def addSortAlignmentJob(self, workflow=None, inputBamFile=None, \
					outputBamFile=None,\
					SortSamFilesJava=None, SortSamJar=None,\
					parentJobLs=None, extraDependentInputLs=None, \
					extraArguments=None, job_max_memory = 2500, transferOutput=False, \
					walltime=180, needBAMIndexJob=True, **keywords):
		"""
		2013.04.09 added argument needBAMIndexJob
		2013.04.05 moved from ShortRead2AlignmentWorkflow
		2012.9.19
		"""
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementData.memRequirement
		javaMemRequirement = memRequirementData.memRequirementInStr
		
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', SortSamJar,\
							"SORT_ORDER=coordinate", "I=", inputBamFile, \
							"O=", outputBamFile, "VALIDATION_STRINGENCY=LENIENT"]
					#not including 'SORT_ORDER=coordinate'
					#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexJar won't fail.)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputBamFile, SortSamJar])
		
		job= self.addGenericJob(executable=SortSamFilesJava, inputFile=None,\
							outputFile=None, outputArgumentOption="-o", \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraOutputLs=[outputBamFile],\
							transferOutput=transferOutput, \
							extraArgumentList=extraArgumentList, \
							job_max_memory=memRequirementData.memRequirement, \
							walltime=walltime, **keywords)
		if needBAMIndexJob:
			# add the index job on the bam file
			bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
						BuildBamIndexJar=self.BuildBamIndexJar, \
						inputBamF=job.output, parentJobLs=[job], \
						transferOutput=transferOutput, job_max_memory=2500, walltime=max(50, walltime/2))
		else:
			bamIndexJob = None
		job.bamIndexJob = bamIndexJob
		return job
	
	def addReadGroupInsertionJob(self, workflow=None, individual_alignment=None, inputBamFile=None, \
								outputBamFile=None,\
								addOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None,\
								parentJobLs=None, extraDependentInputLs=None, \
								extraArguments=None, job_max_memory = 2500, transferOutput=False, walltime=180, max_walltime=1200, \
								needBAMIndexJob=True, **keywords):
		"""
		2013.04.09 moved from ShortRead2AlignmentWorkflow.py
			added argument needBAMIndexJob and the bamIndexJob
		2012.9.19 split out of addAlignmentJob()
		"""
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementData.memRequirement
		javaMemRequirement = memRequirementData.memRequirementInStr
		
		# add RG to this bam
		sequencer = individual_alignment.individual_sequence.sequencer
		read_group = individual_alignment.getReadGroup()	#2012.9.19
		if sequencer.short_name=='454':
			platform_id = 'LS454'
		elif sequencer.short_name=='GA':
			platform_id = 'ILLUMINA'
		else:
			platform_id = 'ILLUMINA'
		# the add-read-group job
		
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', AddOrReplaceReadGroupsJar,\
							"INPUT=", inputBamFile,\
							'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
							'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
							'OUTPUT=', outputBamFile, "VALIDATION_STRINGENCY=LENIENT"]
					#not including 'SORT_ORDER=coordinate'
					#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexJar won't fail.)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputBamFile, AddOrReplaceReadGroupsJar])
		
		job= self.addGenericJob(executable=addOrReplaceReadGroupsJava, inputFile=None,\
							outputFile=None, outputArgumentOption="-o", \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraOutputLs=[outputBamFile],\
							transferOutput=transferOutput, \
							extraArgumentList=extraArgumentList, \
							job_max_memory=memRequirementData.memRequirement, walltime=walltime, max_walltime=max_walltime,\
							**keywords)
		if needBAMIndexJob:
			# add the index job on the bam file
			bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
						BuildBamIndexJar=self.BuildBamIndexJar, \
						inputBamF=job.output, parentJobLs=[job], \
						transferOutput=transferOutput, job_max_memory=job_max_memory)
		else:
			bamIndexJob = None
		job.bamIndexJob = bamIndexJob	#2013.04.09
		return job
	
	def addSelectAlignmentJob(self, executable=None, inputFile=None, \
							outputFile=None, region=None, \
							parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
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
		if executable is None:
			executable = self.samtools
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[outputFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		if needBAMIndexJob:
			# add the index job on the bam file
			bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
						BuildBamIndexJar=self.BuildBamIndexJar, \
						inputBamF=job.output, parentJobLs=[job], \
						transferOutput=transferOutput, job_max_memory=job_max_memory)
		else:
			bamIndexJob = None
		job.bamIndexJob = bamIndexJob	#2013.04.09
		return job, bamIndexJob
	
	def getVCFFileID2object(self, inputDir):
		"""
		2013.05.23 bugfix, include the tbi file
		2011-12-1
		"""
		sys.stderr.write("Getting all vcf files from %s "%(inputDir))
		vcfFileID2object = {}
		for inputFname in os.listdir(inputDir):
			inputAbsPath = os.path.join(os.path.abspath(inputDir), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False) and not NextGenSeq.isVCFFileEmpty(inputAbsPath):
				vcfIndexFname = '%s.tbi'%(inputAbsPath)
				fileID = Genome.getChrFromFname(inputFname)
				if not os.path.isfile(vcfIndexFname):	#does not exist, pass on a None structure
					vcfIndexFname = None
				if fileID in vcfFileID2object:
					sys.stderr.write("fileID %s already has value (%s) in dictionary. but now a 2nd file %s overwrites previous value.\n"%\
									(fileID, vcfFileID2object.get(fileID), inputFname))
				vcfFileID2object[fileID] = PassingData(vcfFilePath=inputAbsPath, vcfIndexFilePath=vcfIndexFname)
		sys.stderr.write("  found %s files.\n"%(len(vcfFileID2object)))
		return vcfFileID2object
	
	def constructGenomeFileRBTreeByFilenameInterval(self, jobDataStructure=None, chr2size=None):
		"""
		2013.09.18
		"""
		sys.stderr.write("Constructing genome file RBtree from %s genome file data ... "%(len(jobDataStructure.jobDataLs)))
		genomeFileRBDict = RBDict()
		for jobData in jobDataStructure.jobDataLs:
			filenameParseData = Genome.parseChrStartStopFromFilename(filename=jobData.file.name, chr2size=chr2size)
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=filenameParseData.chromosome, \
							span_ls=[filenameParseData.start, filenameParseData.stop], \
							min_reciprocal_overlap=0.000000000000001,)	#any overlap is an overlap
			if genomeFileRBDict.has_key(segmentKey):
				sys.stderr.write("\n \t Error: current jobData file %s overlaps node taken by previous jobData file %s.\n"%\
								(jobData.file.name, genomeFileRBDict[segmentKey].file.name))
				raise
			genomeFileRBDict[segmentKey] = jobData
		sys.stderr.write(" %s nodes in the tree. "%(len(genomeFileRBDict)))
		
		return genomeFileRBDict
	
	
	def extractVCFSitesIntoBEDFile(self, inputVCFFolder=None, outputFname=None):
		"""
		2013.05.29
			this function is for genotyping at existing sites.
			BED format is 0-based, tab-delimited, no header, stop position is not inclusive.
				i.e. to describe a single SNP at Contig994 and position=3446. 
					Contig994	3445	3446
				
			
		"""
		sys.stderr.write("Extracting sites from vcf folder %s into file %s ...\n"%(inputVCFFolder, outputFname))
		no_of_vcfFiles = 0
		no_of_loci = 0
		outputF = MatrixFile(inputFname=outputFname, openMode='w', delimiter='\t')
		for inputFname in os.listdir(inputVCFFolder):
			inputAbsPath = os.path.join(os.path.abspath(inputVCFFolder), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False):
				vcfFile= VCFFile(inputFname=inputAbsPath)
				no_of_vcfFiles += 1
				for vcfRecord in vcfFile:
					no_of_loci += 1
					outputF.writerow([vcfRecord.chromosome, vcfRecord.position-1, vcfRecord.position])
				sys.stderr.write("%s%s\t%s "%('\x08'*40, no_of_vcfFiles, no_of_loci))
		outputF.close()
		sys.stderr.write(" %s files and %s loci.\n"%(no_of_vcfFiles, no_of_loci))
	
	
	def addCheckTwoVCFOverlapJob(self, workflow=None, executable=None, vcf1=None, vcf2=None, chromosome=None, chrLength=None, \
					outputFile=None, perSampleConcordanceOutputFile=None, overlapSiteOutputFile=None, parentJobLs=None, \
					extraDependentInputLs=None, transferOutput=False, extraArguments=None, job_max_memory=1000, \
					**keywords):
		"""
		2013.09.17 deal with CheckTwoVCFOverlapCC
		2013.09.10 three output files have explicit arguments now.
		2013.05.20 bugfix
		2012.8.16
			now perSampleMatchFraction output is in a separate output file.
		2011.12.9
		"""
		if workflow is None:
			workflow = self
		extraOutputLs = []
		extraArgumentList = []
		key2ObjectForJob = {}
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 		
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		if vcf2:
			if executable==workflow.CheckTwoVCFOverlapCC:
				extraArgumentList.extend(["-i", vcf2])
			else:	#the python version CheckTwoVCFOverlap.py
				extraArgumentList.extend(["-j", vcf2])
			extraDependentInputLs.append(vcf2)
		if chromosome:
			extraArgumentList.extend(["--chromosome", chromosome])
		if chrLength:
			extraArgumentList.append("--chrLength %s"%(chrLength))
		if perSampleConcordanceOutputFile:
			extraArgumentList.extend(["--perSampleConcordanceOutputFname", perSampleConcordanceOutputFile])
			#suffixAndNameTupleList.append(['_perSample.tsv','perSample'])
			key2ObjectForJob['perSampleFile'] = perSampleConcordanceOutputFile
			extraOutputLs.append(perSampleConcordanceOutputFile)
		if overlapSiteOutputFile:
			extraArgumentList.extend(["--overlappingSitesOutputFname", overlapSiteOutputFile])
			#suffixAndNameTupleList.append(['_overlapSitePos.tsv','overlapSitePos'])
			key2ObjectForJob['overlapSitePosFile'] = overlapSiteOutputFile
			extraOutputLs.append(overlapSiteOutputFile)
		if extraArguments:
			extraArgumentList.append(extraArguments)

		#self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, suffixAndNameTupleList=suffixAndNameTupleList, \
		#											extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
		job= self.addGenericJob(executable=executable, inputFile=vcf1, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
						key2ObjectForJob=key2ObjectForJob, **keywords)
		#if outputFnamePrefix:
		#	job.overlapSitePosF = overlapSitePosF
		return job
	
	def addFilterVCFByDepthJob(self, workflow=None, FilterVCFByDepthJava=None, GenomeAnalysisTKJar=None, \
							refFastaFList=None, inputVCFF=None, outputVCFF=None, outputSiteStatF=None, \
							parentJobLs=[], alnStatForFilterF=None, \
							job_max_memory=1000, extraDependentInputLs=[], onlyKeepBiAllelicSNP=False, \
							namespace=None, version=None, transferOutput=False, **keywords):
		"""
		2013.04.09 added the VCFIndexFile in output (transfer=False)
		2011.12.20
			moved from FilterVCFPipeline.py
			add argument transferOutput, outputSiteStatF
		"""
		# Add a mkdir job for any directory.
		filterByDepthJob = Job(namespace=getattr(self, 'namespace', namespace), name=FilterVCFByDepthJava.name, \
							version=getattr(self, 'version', version))
		refFastaF = refFastaFList[0]
		filterByDepthJob.addArguments("-Xmx%sm"%(job_max_memory), "-jar", GenomeAnalysisTKJar, "-R", refFastaF, "-T FilterVCFByDepth", \
							"--variant", inputVCFF, "-depthFname", alnStatForFilterF)
		self.addJobUse(filterByDepthJob, file=GenomeAnalysisTKJar, transfer=True, register=True, link=Link.INPUT)
		if outputVCFF:
			filterByDepthJob.addArguments("-o", outputVCFF)
			filterByDepthJob.uses(outputVCFF, transfer=transferOutput, register=True, link=Link.OUTPUT)
			filterByDepthJob.output = outputVCFF
			
			VCFIndexFile = File('%s.idx'%(outputVCFF.name))
			self.addJobUse(job=filterByDepthJob, file=VCFIndexFile, transfer=False, register=True, link=Link.OUTPUT)
		
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
		for inputF in extraDependentInputLs:
			filterByDepthJob.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		self.addJob(filterByDepthJob)
		for parentJob in parentJobLs:
			self.depends(parent=parentJob, child=filterByDepthJob)
		yh_pegasus.setJobProperRequirement(filterByDepthJob, job_max_memory=job_max_memory)
		self.no_of_jobs += 1
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
		
		extraArgumentList = [self.vcftoolsExecutableFile]
		extraDependentInputLs.append(self.vcftoolsExecutableFile)
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
							parentJobLs=None, job_max_memory=100, extraDependentInputLs=None, \
							transferOutput=False, **keywords):
		"""
		Examples:
			#tabix retrieve job
			outputF = File(os.path.join(trioInconsistencyDir, '%s.trioInconsistency.tsv'%chromosome))
			tabixRetrieveJob = self.addTabixRetrieveJob(executable=workflow.tabixRetrieve, \
							tabixPath=self.tabixPath, \
							inputF=trioInconsistencyByPosistionF, outputF=outputF, \
							regionOfInterest=chromosome, includeHeader=True,\
							parentJobLs=[trioInconsistencyDirJob], job_max_memory=100, \
							extraDependentInputLs=[trioInconsistencyByPosistion_tbi_F], \
							transferOutput=False)
			
			#tabix index job
			#index .vcf.gz, output of beagle, without index, GATK can't work on gzipped vcf
			tabixIndexFile = File('%s.tbi'%(beagleJob.output.name))
			tabixJob = self.addGenericJob(executable=self.tabix, \
							inputFile=beagleJob.output, inputArgumentOption="",\
							outputFile=None, outputArgumentOption="-o", extraDependentInputLs=None, \
							extraOutputLs=[beagleJob.output, tabixIndexFile], transferOutput=False, \
							frontArgumentList=["-p vcf"], extraArguments=None, extraArgumentList=None, \
							parentJobLs=[beagleJob, outputDirJob], no_of_cpus=None, \
						job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
								baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
								minJobPropertyValue=2000, maxJobPropertyValue=4000).value,\
						walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
								baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
								minJobPropertyValue=60, maxJobPropertyValue=180).value)

		
		2011.12.20
			The executable should be tabixRetrieve (a tabix shell wrapper).
				vervet/src/shell/tabixRetrieve.sh
			
			http://samtools.sourceforge.net/tabix.shtml
			run something like below to extract data from regionOfInterest out of bgzipped&tabix-indexed file.
				tabix sorted.gff.gz chr1:10,000,000-20,000,000
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = workflow.tabixRetrieve
		extraArgumentList=[]
		extraArgumentList.append(regionOfInterest)
		if includeHeader:
			extraArgumentList.append("-h")
		
		job = self.addGenericJob(executable=executable, frontArgumentList=[tabixPath], \
					inputFile=inputF, inputArgumentOption=None,\
					outputFile=outputF, outputArgumentOption=None, \
					extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=[], transferOutput=transferOutput, \
					extraArgumentList=extraArgumentList, parentJobLs=parentJobLs,\
					no_of_cpus=None, job_max_memory = job_max_memory,\
					walltime= 60)
		
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
	
	def getTopNumberOfContigs(self, contigMaxRankBySize=100, contigMinRankBySize=1, tax_id=60711, sequence_type_id=9,\
							version=None, chromosome_type_id=0, outdated_index=0):
		"""
		2013.3.14 added argument version
		2013.2.15 added argument chromosome_type_id
		chromosome_type_id: what type of chromosomes to be included, same as table genome.chromosome_type.
			0: all, 1: autosomes, 2: X, 3:Y, 4: mitochondrial
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
		
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=contigMaxRankBySize, contigMinRankBySize=contigMinRankBySize, \
							tax_id=tax_id, sequence_type_id=sequence_type_id, \
							version=version, chromosome_type_id=chromosome_type_id, outdated_index=outdated_index)
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
					sys.stderr.write("Error in handling contig %s.\n"%(contig))
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			dc = new_dc
		sys.stderr.write(" %s contigs left.\n"%(len(dc)))
		return dc
	
	def getChr2IntervalDataLsBySplitBEDFile(self, intervalFname=None, noOfLinesPerUnit=2000, folderName=None, parentJobLs= None):
		"""
		2013.05.29 added span and chromosomeSize to final returned data
		2012.08.09 update it so that the interval encompassing all lines in one block/unit is known.
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
		sys.stderr.write("Splitting %s into blocks, each block with %s lines ... "%(intervalFname, noOfLinesPerUnit))
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
		chr2MaxStopPosition = {}	#2013.05.29
		for row in reader:
			lineNumber += 1
			chromosome, start, stop = row[:3]
			start = int(start)	#0-based, starting base
			stop = int(stop)	#0-based, stopping base but not inclusive, i.e. [start, stop)
			#2013.05.29
			if chromosome not in chr2MaxStopPosition:
				chr2MaxStopPosition[chromosome] = stop
			elif stop>chr2MaxStopPosition[chromosome]:
				chr2MaxStopPosition[chromosome] = stop
			
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
				intervalFileBasenameSignature = '%s_%s_%s'%(chr, blockStartLineNumber, blockStopLineNumber)
				if chr not in chr2IntervalDataLs:
					chr2IntervalDataLs[chr] = []
				span = blockStopLineNumber - blockStartLineNumber + 1	#how many loci, rather than the region size
				intervalData = PassingData(file=blockIntervalFile, intervalFileBasenameSignature=intervalFileBasenameSignature, interval=interval,\
										overlapInterval=interval, span=span, chromosomeSize=chr2MaxStopPosition.get(chr), \
										chr=chr, jobLs=[blockIntervalJob], job=blockIntervalJob)
				chr2IntervalDataLs[chr].append(intervalData)
				counter += 1
		sys.stderr.write("%s intervals and %s SelectLineBlockFromFile jobs.\n"%(counter, counter))
		return chr2IntervalDataLs
	
	def getChr2IntervalDataLsBySplitChrSize(self, chr2size=None, intervalSize=None, intervalOverlapSize=None):
		"""
		2013.04.09 added interval-span, chromosomeSize in returning data, 
		2012.7.30
		"""
		sys.stderr.write("Splitting %s references into intervals of %s bp (overlap=%s) ... "%(len(chr2size), intervalSize,\
																		intervalOverlapSize))
		chr2IntervalDataLs = {}
		counter =0
		for chromosome, chromosomeSize in chr2size.iteritems():
			no_of_intervals = max(1, int(math.ceil(chromosomeSize/float(intervalSize)))-1)
			for i in range(no_of_intervals):
				originalStartPos = i*intervalSize + 1
				#to render adjacent intervals overlapping because trioCaller uses LD
				overlapStart = max(1, originalStartPos-intervalOverlapSize)
				if i<no_of_intervals-1:
					originalStopPos = min((i+1)*intervalSize, chromosomeSize)
				else:	#last chunk, include bp till the end
					originalStopPos = chromosomeSize
				#to render adjacent intervals overlapping because trioCaller uses LD
				overlapStop = min(chromosomeSize, originalStopPos+intervalOverlapSize)
				
				if chromosome not in chr2IntervalDataLs:
					chr2IntervalDataLs[chromosome] = []
				intervalData = IntervalData(chromosome=chromosome, chromosomeSize=chromosomeSize,\
								start=originalStartPos, stop=originalStopPos,\
								overlapStart=overlapStart, overlapStop=overlapStop)
				chr2IntervalDataLs[chromosome].append(intervalData)
				counter += 1
		sys.stderr.write("%s intervals.\n"%(counter))
		return chr2IntervalDataLs
	
	def getChr2IntervalDataLsFromDBAlignmentDepthInterval(self, db=None, intervalSize=None, intervalOverlapSize=None,\
									alignmentDepthIntervalMethodShortName=None, alignmentDepthMinFold=0.1, alignmentDepthMaxFold=2, \
									minAlignmentDepthIntervalLength=1000, maxContigID=None, minContigID=None):
		"""
		2013.10.27 bugfix
		2013.08.29
		"""
		method = db.getAlignmentDepthIntervalMethod(short_name =alignmentDepthIntervalMethodShortName)
		minDepth = method.sum_median_depth*float(alignmentDepthMinFold)
		maxDepth = method.sum_median_depth*alignmentDepthMaxFold
		sys.stderr.write("Getting alignment depth intervals (depth %sX median - %sX median depth)=(%s-%s) from db, maxContigID=%s, minContigID=%s  ... \n "%\
					(alignmentDepthMinFold, alignmentDepthMaxFold, minDepth, maxDepth, maxContigID, minContigID))
		counter =0
		
		noOfRawIntervals = 0
		chr2alignmentDepthIntervalData = {}
		for alignment_depth_interval_file in method.alignment_depth_interval_file_ls:
			chromosome = alignment_depth_interval_file.chromosome
			included = True
			if (maxContigID is not None and maxContigID!=0) or (minContigID is not None and minContigID!=0):	#2012.10.27 bugfix, "and" changed to "or"
				try:
					contigID = int(self.getContigIDFromFname(chromosome))
					if (maxContigID is not None and maxContigID!=0) and contigID>maxContigID:
						included = False
					if (minContigID is not None and minContigID!=0) and contigID<minContigID:
						included = False
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
					included = False
			if not included:	#skip this chromosome
				continue
			chromosomeSize = alignment_depth_interval_file.chromosome_size
			inputFname = os.path.join(self.data_dir, alignment_depth_interval_file.path)
			sys.stderr.write("\t %s ..."%(inputFname))
			if os.path.isfile(inputFname):
				reader = AlignmentDepthIntervalFile(inputFname)
				intervalLs = reader.getAllIntervalWithinDepthRange(minDepth=minDepth, maxDepth=maxDepth)
				reader.close()
				noOfRawIntervals += len(intervalLs)
				if chromosome not in chr2alignmentDepthIntervalData:
					chr2alignmentDepthIntervalData[chromosome] = PassingData(chromosomeSize=chromosomeSize, intervalLs=intervalLs)
				else:
					sys.stderr.write("\t Warning: chromosome %s sees new alignment depth interval data from %s.\n"%(chromosome, inputFname))
					chr2alignmentDepthIntervalData[chromosome].intervalLs.extend(intervalLs)
			else:
				sys.stderr.write("\t Warning: interval file %s of chromosome %s does not exist. Ignore.\n"%(inputFname, chromosome))
			sys.stderr.write(" noOfRawIntervals=%s.\n"%(noOfRawIntervals))
		sys.stderr.write(" %s alignment depth intervals covering %s chromosomes.\n"%(noOfRawIntervals, len(chr2alignmentDepthIntervalData)))
		
		sys.stderr.write("Splitting alignment depth intervals (depth %sX median - %sX median depth), minAlignmentDepthIntervalLength=%s) into size=%s intervals, overlap=%s,  ... "%\
					(alignmentDepthMinFold, alignmentDepthMaxFold, minAlignmentDepthIntervalLength, intervalSize, intervalOverlapSize))
		chr2IntervalDataLs = {}
		for chromosome, alignmentDepthIntervalData in chr2alignmentDepthIntervalData.iteritems():
			intervalLs = alignmentDepthIntervalData.intervalLs
			chromosomeSize = alignmentDepthIntervalData.chromosomeSize
			for interval  in intervalLs:
				if interval.length<minAlignmentDepthIntervalLength:	#skip short intervals
					continue
				if intervalSize is None:
					no_of_intervals = 1
				else:
					no_of_intervals = max(1, int(math.ceil(interval.length/float(intervalSize)))-1)
				for i in range(no_of_intervals):
					originalStartPos = interval.start + i*intervalSize 
					#to render adjacent intervals overlapping because trioCaller uses LD
					overlapStart = max(interval.start, originalStartPos-intervalOverlapSize)
					if i<no_of_intervals-1:
						originalStopPos = min(interval.start + (i+1)*intervalSize-1, interval.stop)
					else:	#last chunk, include bp till the end
						originalStopPos = interval.stop
					#to render adjacent intervals overlapping because trioCaller uses LD
					overlapStop = min(interval.stop, originalStopPos+intervalOverlapSize)
					#2013.09.09 overlap is not considered because considering overlap would go into those highly repetitive inter-interval regions.
					
					if chromosome not in chr2IntervalDataLs:
						chr2IntervalDataLs[chromosome] = []
					intervalData = IntervalData(chromosome=chromosome, chromosomeSize=chromosomeSize,\
								start=originalStartPos, stop=originalStopPos)
					chr2IntervalDataLs[chromosome].append(intervalData)
					counter += 1
		sys.stderr.write("%s intervals.\n"%(counter))
		return chr2IntervalDataLs
	
	def addPutStuffIntoDBJob(self, workflow=None, executable=None, \
					inputFile=None, inputArgumentOption="-i", \
					inputFileList=None, \
					outputFile=None, outputArgumentOption="-o",\
					logFile=None, commit=False, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=True, \
					extraArguments=None, extraArgumentList=None,\
					job_max_memory=10, sshDBTunnel=0, **keywords):
		"""
		2013.08.06 added extraArgumentList
		2013.3.24 use addGenericFile2DBJob()
		2012.5.8 add sshDBTunnel
		2012.4.3
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		if extraArgumentList is None:
			extraArgumentList = []
		key2ObjectForJob = {}
		
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if inputFileList:
			extraDependentInputLs.extend(inputFileList)
		job = self.addGenericFile2DBJob(workflow=workflow, executable=executable, \
					inputFile=inputFile, inputArgumentOption=inputArgumentOption, \
					inputFileList=inputFileList, \
					outputFile=outputFile, outputArgumentOption=outputArgumentOption, \
					data_dir=None, logFile=logFile, commit=commit,\
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=None, \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					sshDBTunnel=sshDBTunnel, \
					key2ObjectForJob=key2ObjectForJob, objectWithDBArguments=self, **keywords)
		return job
	
	def addSamtoolsFlagstatJob(self, workflow=None, executable=None, samtoolsExecutableFile=None, \
					inputFile=None, outputFile=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, walltime=120, **keywords):
		"""
		2013.03.25 use pipeCommandOutput2File to get output piped into outputF
		2013.3.24 use addGenericJob()
		2012.4.3
			samtools (sam_stat.c) has been modified so that it could take one more optional argument to store the
				stats that are usually directed to stdout. 
			inputF is bam file. outputF is to store the output.
			
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraDependentInputLs.append(inputFile)
		job = self.addGenericPipeCommandOutput2FileJob(workflow=workflow, executable=executable, \
					executableFile=samtoolsExecutableFile, \
					outputFile=outputFile, \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=['flagstat', inputFile], sshDBTunnel=None,\
					job_max_memory=job_max_memory, walltime=walltime)
		return job
	
	def addAddAlignmentFile2DBJob(self, workflow=None, executable=None, inputFile=None, otherInputFileList=None,\
						individual_alignment_id=None, individual_sequence_id=None,\
						ref_sequence_id=None, \
						alignment_method_id=None,\
						parent_individual_alignment_id=None, \
						mask_genotype_method_id=None,\
						individual_sequence_file_raw_id=None,\
						format=None, local_realigned=0,\
						logFile=False, data_dir=None, \
						parentJobLs=None, \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=True, \
						job_max_memory=2000, walltime=180, sshDBTunnel=False, commit=True, **keywords):
		"""
		2013.04.05 added argument local_realigned
		2012.9.20
			To specify individual_alignment:
				either individual_alignment_id or (parent_individual_alignment_id + mask_genotype_method_id)
				or others 
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = []
		extraOutputLs = []
		key2ObjectForJob = {}
		if otherInputFileList:
			extraDependentInputLs.extend(otherInputFileList)
		if logFile:
			extraArgumentList.extend(['--logFilename', logFile])
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
		if local_realigned is not None:
			extraArgumentList.append("--local_realigned %s"%(local_realigned))
		
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, walltime=walltime, **keywords)
		job.logFile = logFile
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		
		#add all input files to the last (after db arguments,) otherwise, it'll mask others (cuz these don't have options).
		if otherInputFileList:
			for inputFile in otherInputFileList:
				if inputFile:
					job.addArguments(inputFile)
		return job
	
	def addAlignmentMergeJob(self, workflow=None, AlignmentJobAndOutputLs=None, outputBamFile=None,\
					MergeSamFilesJava=None, \
					BuildBamIndexFilesJava=None, \
					mv=None, parentJobLs=None, namespace=None, version=None, transferOutput=False,\
					job_max_memory=7000, walltime=680, **keywords):
		"""
		not certain, but it looks like MergeSamFilesJar does not require the .bai (bam index) file.
		
		2013.03.31 AlignmentJobAndOutputLs has changed it cell structure
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
		if MergeSamFilesJava is None:
			MergeSamFilesJava = workflow.MergeSamFilesJava
		if BuildBamIndexFilesJava is None:
			BuildBamIndexFilesJava = workflow.BuildBamIndexFilesJava
		if mv is None:
			mv = workflow.mv
		
		namespace = getattr(workflow, 'namespace', namespace)
		version = getattr(workflow, 'version', version)
		
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		
		
		if len(AlignmentJobAndOutputLs)>1:
			# 'USE_THREADING=true', threading might be causing process hanging forever (sleep).
			alignmentJobLs = []
			alignmentOutputLs = []
			for AlignmentJobAndOutput in AlignmentJobAndOutputLs:
				if AlignmentJobAndOutput.jobLs:
					alignmentJobLs.extend(AlignmentJobAndOutput.jobLs)
				alignmentOutput = AlignmentJobAndOutput.file
				if alignmentOutput:
					alignmentOutputLs.append(alignmentOutput)
						#2013.04.01
			merge_sam_job = self.addGenericJavaJob(executable=MergeSamFilesJava, jarFile=workflow.MergeSamFilesJar, \
					inputFile=None, inputArgumentOption=None, \
					inputFileList=alignmentOutputLs, argumentForEachFileInInputFileList="INPUT=",\
					outputFile=outputBamFile, outputArgumentOption="OUTPUT=",\
					parentJobLs=alignmentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					frontArgumentList=None, extraArguments=None, extraArgumentList=['SORT_ORDER=coordinate', \
						'ASSUME_SORTED=true', "VALIDATION_STRINGENCY=LENIENT"], extraOutputLs=None, \
					extraDependentInputLs=None, no_of_cpus=None, walltime=walltime, sshDBTunnel=None, **keywords)
		elif len(AlignmentJobAndOutputLs)==1:	#one input file, no samtools merge. use "mv" to rename it instead. should use "cp", then the input would be cleaned by cleaning job.
			AlignmentJobAndOutput = AlignmentJobAndOutputLs[0]
			alignmentJobLs = AlignmentJobAndOutput.jobLs
			alignmentOutput = AlignmentJobAndOutput.file
			#2012.7.4
			sys.stderr.write(" copy (instead of merging small alignment files) due to only one alignment file, from %s to %s.\n"%\
							(alignmentOutput.name, outputBamFile.name))
			#2013.04.01
			merge_sam_job = self.addGenericJob(executable=mv, inputFile=alignmentOutput, inputArgumentOption=None, \
					outputFile=outputBamFile, outputArgumentOption=None, \
					inputFileList=None, argumentForEachFileInInputFileList=None, \
					parentJob=None, parentJobLs=alignmentJobLs, extraDependentInputLs=None, \
					extraOutputLs=None, transferOutput=transferOutput, \
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, objectWithDBArguments=None, no_of_cpus=None, walltime=40)
		else:
			sys.stderr.write("Error: no input for MergeSamFilesJar to output %s.\n"%(outputBamFile))
			raise
		#assign output
		merge_sam_job.output = outputBamFile
		#2012.9.21
		if parentJobLs:
			for parentJob in parentJobLs:
				workflow.depends(parent=parentJob, child=merge_sam_job)
						
		# add the index job on the merged bam file
		bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, \
					BuildBamIndexJar=workflow.BuildBamIndexJar, \
					inputBamF=outputBamFile,\
					parentJobLs=[merge_sam_job], namespace=namespace, version=version,\
					transferOutput=transferOutput, job_max_memory=3000, walltime=max(180, int(walltime/3)))
		merge_sam_job.bamIndexJob = bamIndexJob
		return merge_sam_job, bamIndexJob
	
	
	def addMergeVCFReplicateGenotypeColumnsJob(self, workflow=None, executable=None, GenomeAnalysisTKJar=None, \
						inputF=None, outputF=None, replicateIndividualTag=None, \
						debugHaplotypeDistanceFile=None, \
						debugMajoritySupportFile=None,\
						refFastaFList=None, parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=2000, analysis_type='MergeVCFReplicateHaplotypes',\
						walltime=None, **keywords):
		"""
		2013.06.21 use addGATKJob() instead
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
		extraArgumentList = ['--onlyKeepBiAllelicSNP', "--replicateIndividualTag %s"%(replicateIndividualTag)]
		if debugHaplotypeDistanceFile:
			extraArgumentList.extend(["--debugHaplotypeDistanceFname", debugHaplotypeDistanceFile])
		if debugMajoritySupportFile:
			extraArgumentList.extend(["--debugMajoritySupportFname", debugMajoritySupportFile])
		
		job = self.addGATKJob(executable=executable, GenomeAnalysisTKJar=GenomeAnalysisTKJar, \
						GATKAnalysisType=analysis_type,\
						inputFile=inputF, inputArgumentOption="--variant:VCF", \
						refFastaFList=refFastaFList, inputFileList=None,\
						argumentForEachFileInInputFileList=None,\
						interval=None, outputFile=outputF, \
						parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
						frontArgumentList=None, extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
						extraOutputLs=None, \
						extraDependentInputLs=extraDependentInputLs, no_of_cpus=None, walltime=walltime)
		
		return job
	
	def isThisAlignmentComplete(self, individual_alignment=None, data_dir=None, returnFalseIfInexitentFile=False,**keywords):
		"""
		this is to check whether the new and to-be-generated alignment is completed already or not.
			in contrast to isThisInputAlignmentComplete().
		2013.05.04 added argument returnFalseIfInexitentFile
			
		2013.04.09 for subsequent children to override
		"""
		return self.db.isThisAlignmentComplete(individual_alignment=individual_alignment, data_dir=data_dir,\
											returnFalseIfInexitentFile=returnFalseIfInexitentFile, **keywords)

	def isThisInputAlignmentComplete(self, individual_alignment=None, data_dir=None, returnFalseIfInexitentFile=True, \
									**keywords):
		"""
		2013.05.04 this is used to check whether an input (to be worked on by downstream programs) is completed or not.
			watch returnFalseIfInexitentFile is True (because you need the file for an input)
		
		"""
		return self.db.isThisAlignmentComplete(individual_alignment=individual_alignment, data_dir=data_dir, \
											returnFalseIfInexitentFile=returnFalseIfInexitentFile, **keywords)
	
if __name__ == '__main__':
	main_class = AbstractNGSWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()