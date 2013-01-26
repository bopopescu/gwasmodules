#!/usr/bin/env python
"""
2013.1.25 an abstract class for pegasus workflows that work on alignment files (already aligned).
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, utils
from AbstractNGSWorkflow import AbstractNGSWorkflow

class AbstractAlignmentWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	#option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('intervalOverlapSize', 1, int): [300000, 'U', 1, 'overlap #bps/#loci between adjacent intervals from one contig/chromosome,\
				only used for TrioCaller, not for SAMtools/GATK', ],\
						('intervalSize', 1, int): [5000000, 'Z', 1, '#bps/#loci for adjacent intervals from one contig/chromosome (alignment or VCF)', ],\
						})
	commonAlignmentWorkflowOptionDict = {
						('ind_seq_id_ls', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						("sequence_filtered", 0, int): [None, 'Q', 1, 'To filter alignments. None: whatever; 0: unfiltered sequences, 1: filtered sequences: 2: ...'],\
						("alignment_method_id", 0, int): [None, 'G', 1, 'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
						("site_id_ls", 0, ): ["", 'S', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
						("country_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
						("tax_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
						('defaultSampleAlignmentDepth', 1, int): [10, '', 1, "when database doesn't have median_depth info for one alignment, use this number instead.", ],\
						('individual_sequence_file_raw_id_type', 1, int): [1, '', 1, "1: only all-library-fused libraries,\n\
		2: only library-specific alignments,\n\
		3: both all-library-fused and library-specific alignments", ],\
						}
	option_default_dict.update(commonAlignmentWorkflowOptionDict)
	partitionWorkflowOptionDict= {
						("needFastaIndexJob", 0, int): [0, 'A', 0, 'need to add a reference index job by samtools?'],\
						("needFastaDictJob", 0, int): [0, 'B', 0, 'need to add a reference dict job by picard CreateSequenceDictionary.jar?'],\
						("selectedRegionFname", 0, ): ["", 'R', 1, 'the file is in bed format, tab-delimited, chr start stop.\
		used to restrict SAMtools/GATK to only make calls at this region. \
		start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.\
		This overrides the contig/chromosome selection approach defined by --contigMaxRankBySize and --contigMinRankBySize. \
		This file would be split into maxNoOfRegionsPerJob lines.'],\
						('maxNoOfRegionsPerJob', 1, int): [5000, 'K', 1, 'Given selectedRegionFname, this dictates the maximum number of regions each job would handle,\
		The actual number could be lower because the regions are first grouped into chromosomes. If one chromosome has <maxNoOfRegionsPerJob, then that job handles less.', ],\
						}
	
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		listArgumentName_data_type_ls = [('ind_seq_id_ls', int), ("ind_aln_id_ls", int), \
								("site_id_ls", int), ('country_id_ls', int), ('tax_id_ls', int)]
		listArgumentName2hasContent = self.processListArguments(listArgumentName_data_type_ls, emptyContent=[])
	
	
	
	def preReduce(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachChromosome(self, workflow=None, alignmentData=None, chromosome=None,\
				VCFFile=None, passingData=None, reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def map(self, workflow=None, alignmentData=None, intervalData=None,\
		VCFFile=None, passingData=None, mapEachChromosomeData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachInterval(self, **keywords):
		"""
		2012.9.22 link to map()
		"""
		return self.map(**keywords)
	
	
	def linkMapToReduce(self, workflow=None, mapEachIntervalData=None, preReduceReturnData=None, passingData=None, transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData

	def mapEachAlignment(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.22
			similar to reduceBeforeEachAlignmentData() but for mapping programs that run on one alignment each.
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def reduceAfterEachChromosome(self, workflow=None, chromosome=None, passingData=None, transferOutput=True, \
								mapEachIntervalDataLs=None, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
		return returnData
	
	def reduceBeforeEachAlignment(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9 setup some reduce jobs before loop over all intervals of one alignment begins.
			these reduce jobs will collect stuff from each map() job.
			the link will be established in linkMapToReduce().
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def reduceAfterEachAlignment(self, workflow=None, passingData=None, mapEachChromosomeDataLs=None,\
								reduceAfterEachChromosomeDataLs=None,\
								transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachChromosomeDataLs = mapEachChromosomeDataLs
		returnData.reduceAfterEachChromosomeDataLs = reduceAfterEachChromosomeDataLs
		return returnData

	def reduce(self, workflow=None, passingData=None, reduceAfterEachAlignmentDataLs=None,
			transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.reduceAfterEachAlignmentDataLs = reduceAfterEachAlignmentDataLs
		return returnData
	
	def mapReduceWholeAlignment(self, workflow=None, alignmentData=None, passingData=None, \
						chrIDSet=None, chr2IntervalDataLs=None, chr2VCFFile=None, \
						outputDirPrefix=None, transferOutput=False, **keywords):
		"""
		2013.1.25
		"""
		pass
	
	def prePreprocess(self, chr2IntervalDataLs=None, **keywords):
		"""
		2013.1.25
			chr2VCFFile is None.
		"""
		chrIDSet = set(chr2IntervalDataLs.keys())
		return PassingData(chrIDSet=chrIDSet, chr2VCFFile=None)
	
	def addAllJobs(self, workflow=None, alignmentDataLs=None, chr2IntervalDataLs=None, samtools=None, \
				genomeAnalysisTKJar=None, \
				mergeSamFilesJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mv=None, \
				refFastaFList=None, \
				needFastaIndexJob=False, needFastaDictJob=False, \
				dataDir=None, no_of_gatk_threads = 1, \
				outputDirPrefix="", transferOutput=True, **keywords):
		"""
		2012.7.26
		"""
		prePreprocessData = self.prePreprocess(chr2IntervalDataLs=chr2IntervalDataLs, **keywords)
		chrIDSet = prePreprocessData.chrIDSet
		chr2VCFFile = prePreprocessData.chr2VCFFile
		
		sys.stderr.write("Adding alignment (+VCF)-related jobs for %s chromosomes/contigs ..."%(len(chrIDSet)))
		refFastaF = refFastaFList[0]
		
		topOutputDir = "%sMap"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		if needFastaDictJob:	# the .dict file is required for GATK
			fastaDictJob = self.addRefFastaDictJob(workflow, createSequenceDictionaryJava=createSequenceDictionaryJava, \
												refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
		else:
			fastaDictJob = None
			refFastaDictF = None
		
		if needFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
		else:
			fastaIndexJob = None
			refFastaIndexF = None
		
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		#2012.9.22 AlignmentJobAndOutputLs is a relic.
		#	but it's similar to mapEachIntervalDataLs but designed for addAlignmentMergeJob(),
		#	so AlignmentJobAndOutputLs gets re-set for every alignment.
		# 	mapEachAlignmentDataLs is never reset.
		#	mapEachChromosomeDataLs is reset right after a new alignment is chosen.
		#	mapEachIntervalDataLs is reset right after each chromosome is chosen.
		#	all reduce dataLs never gets reset.
		passingData = PassingData(AlignmentJobAndOutputLs=[], bamFnamePrefix=None, topOutputDirJob=topOutputDirJob,\
					outputDirPrefix=outputDirPrefix, refFastaFList=refFastaFList, \
					
					mapEachAlignmentData = None,\
					mapEachChromosomeData=None, \
					mapEachIntervalData=None,\
					reduceBeforeEachAlignmentData = None, \
					reduceAfterEachAlignmentData=None,\
					reduceAfterEachChromosomeData=None,\
					
					mapEachAlignmentDataLs = [],\
					mapEachChromosomeDataLs=[], \
					mapEachIntervalDataLs=[],\
					reduceBeforeEachAlignmentDataLs = [], \
					reduceAfterEachAlignmentDataLs=[],\
					reduceAfterEachChromosomeDataLs=[],\
					
					gzipReduceAfterEachChromosomeFolderJob=None,\
					gzipReduceBeforeEachAlignmentFolderJob = None,\
					gzipReduceAfterEachAlignmentFolderJob = None,\
					gzipPreReduceFolderJob = None,\
					gzipReduceFolderJob=None,\
					)
		preReduceReturnData = self.preReduce(workflow=workflow, passingData=passingData, transferOutput=False,\
											**keywords)
		passingData.preReduceReturnData = preReduceReturnData
		
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			parentJobLs = alignmentData.jobLs
			bamF = alignmentData.bamF
			baiF = alignmentData.baiF
			
			bamFnamePrefix = alignment.getReadGroup()
			
			passingData.AlignmentJobAndOutputLs = []
			passingData.bamFnamePrefix = bamFnamePrefix
			passingData.individual_alignment = alignment
			
			mapEachAlignmentData = self.mapEachAlignment(workflow=workflow, passingData=passingData, transferOutput=False, \
												preReduceReturnData=preReduceReturnData, **keywords)
			passingData.mapEachAlignmentDataLs.append(mapEachAlignmentData)
			passingData.mapEachAlignmentData = mapEachAlignmentData
			
			reduceBeforeEachAlignmentData = self.reduceBeforeEachAlignment(workflow=workflow, passingData=passingData, \
													preReduceReturnData=preReduceReturnData, transferOutput=False, \
													**keywords)
			passingData.reduceBeforeEachAlignmentData = reduceBeforeEachAlignmentData
			passingData.reduceBeforeEachAlignmentDataLs.append(reduceBeforeEachAlignmentData)
			
			
			self.mapReduceWholeAlignment(workflow=workflow, alignmentData=alignmentData, passingData=passingData, \
							chrIDSet=chrIDSet, chr2IntervalDataLs=chr2IntervalDataLs, chr2VCFFile=chr2VCFFile, \
							outputDirPrefix=outputDirPrefix, transferOutput=transferOutput)
			
			reduceAfterEachAlignmentData = self.reduceAfterEachAlignment(workflow=workflow, \
												mapEachAlignmentData=mapEachAlignmentData,\
												mapEachChromosomeDataLs=passingData.mapEachChromosomeDataLs,\
												reduceAfterEachChromosomeDataLs=passingData.reduceAfterEachChromosomeDataLs,\
												passingData=passingData, \
												transferOutput=False, dataDir=dataDir, **keywords)
			passingData.reduceAfterEachAlignmentData = reduceAfterEachAlignmentData
			passingData.reduceAfterEachAlignmentDataLs.append(reduceAfterEachAlignmentData)
			
			gzipReduceBeforeEachAlignmentData = self.addGzipSubWorkflow(workflow=workflow, \
						inputData=reduceBeforeEachAlignmentData, transferOutput=transferOutput,\
						outputDirPrefix="%sreduceBeforeEachAlignment"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceBeforeEachAlignmentFolderJob, report=False)
			passingData.gzipReduceBeforeEachAlignmentFolderJob = gzipReduceBeforeEachAlignmentData.topOutputDirJob
			
			gzipReduceAfterEachAlignmentData = self.addGzipSubWorkflow(workflow=workflow, \
						inputData=reduceAfterEachAlignmentData, transferOutput=transferOutput,\
						outputDirPrefix="%sreduceAfterEachAlignment"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceAfterEachAlignmentFolderJob, \
						report=False)
			passingData.gzipReduceAfterEachAlignmentFolderJob = gzipReduceAfterEachAlignmentData.topOutputDirJob
		reduceReturnData = self.reduce(workflow=workflow, passingData=passingData, \
							mapEachAlignmentData=mapEachAlignmentData, \
							reduceAfterEachAlignmentDataLs=passingData.reduceAfterEachAlignmentDataLs,\
							**keywords)
		passingData.reduceReturnData = reduceReturnData
		
		
		#2012.9.18 gzip the final output
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=preReduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%spreReduce"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipPreReduceFolderJob, \
						report=False)
		passingData.gzipPreReduceFolderJob = newReturnData.topOutputDirJob
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=reduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%sreduce"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceFolderJob, \
						report=False)
		passingData.gzipReduceFolderJob = newReturnData.topOutputDirJob
		
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return returnData
	
	def registerAlignmentAndItsIndexFile(self, workflow=None, alignmentLs=None, dataDir=None, checkFileExistence=True):
		"""
		2012.9.18 copied from AlignmentToCallPipeline.py
		2012.6.12
			add argument checkFileExistence
		2012.1.9
			register the input alignments and return in a data structure usd by several other functions
		"""
		sys.stderr.write("Registering %s alignments ..."%(len(alignmentLs)))
		returnData = []
		for alignment in alignmentLs:
			inputFname = os.path.join(dataDir, alignment.path)
			input = File(alignment.path)	#relative path, induces symlinking or stage-in
			baiFilepath = '%s.bai'%(inputFname)
			if checkFileExistence and (not os.path.isfile(inputFname) or not os.path.isfile(baiFilepath)):
				continue
			input.addPFN(PFN("file://" + inputFname, workflow.input_site_handler))
			workflow.addFile(input)
			baiF = File('%s.bai'%alignment.path)
			baiF.addPFN(PFN("file://" + baiFilepath, workflow.input_site_handler))
			workflow.addFile(baiF)
			alignmentData = PassingData(alignment=alignment, jobLs = [], bamF=input, baiF=baiF)
			returnData.append(alignmentData)
		sys.stderr.write("Done.\n")
		return returnData
	
	
	def registerCustomExecutables(self, workflow=None):
		
		"""
		"""
		AbstractVervetWorkflow.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		executableClusterSizeMultiplierList = []
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	
	def run(self):
		"""
		2013.1.25
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		if not self.dataDir:
			self.dataDir = self.db.data_dir
		
		if not self.localDataDir:
			self.localDataDir = self.db.data_dir
		
		chr2size = self.getTopNumberOfContigs(contigMaxRankBySize=self.contigMaxRankBySize, contigMinRankBySize=self.contigMinRankBySize)
		#chr2size = set(['Contig149'])	#temporary when testing Contig149
		#chr2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		chrLs = chr2size.keys()
		chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=chr2size, \
													intervalSize=self.intervalSize, \
													intervalOverlapSize=self.intervalOverlapSize)
		
		alignmentLs = self.db.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, dataDir=self.localDataDir,\
										individual_sequence_file_raw_id_type=self.individual_sequence_file_raw_id_type,\
										country_id_ls=self.country_id_ls, tax_id_ls=self.tax_id_ls)
		alignmentLs = self.db.filterAlignments(alignmentLs, sequence_filtered=self.sequence_filtered, \
									individual_site_id_set=set(self.site_id_ls),\
									mask_genotype_method_id=None, parent_individual_alignment_id=None,\
									country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls))
		
		workflow = self.initiateWorkflow()
		
		refFastaFList = self.getReferenceSequence()
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables()
		alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs=alignmentLs, dataDir=self.dataDir)
		
		self.addAllJobs(workflow=workflow, alignmentDataLs=alignmentDataLs, \
				chr2IntervalDataLs=chr2IntervalDataLs, samtools=workflow.samtools, \
				genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
				mergeSamFilesJar=workflow.mergeSamFilesJar, \
				createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar,\
				mv=workflow.mv, \
				refFastaFList=refFastaFList,\
				needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
				dataDir=self.dataDir, no_of_gatk_threads = 1, transferOutput=True,\
				outputDirPrefix=self.pegasusFolderName)
		
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)

if __name__ == '__main__':
	main_class = AbstractAlignmentWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()