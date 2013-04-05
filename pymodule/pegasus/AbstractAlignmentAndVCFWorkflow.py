#!/usr/bin/env python
"""
2013.1.25 an abstract class for pegasus workflows that work on alignment & VCF files.
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, utils
from AbstractAlignmentWorkflow import AbstractAlignmentWorkflow
from AbstractVCFWorkflow import AbstractVCFWorkflow

class AbstractAlignmentAndVCFWorkflow(AbstractAlignmentWorkflow, AbstractVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractAlignmentWorkflow.option_default_dict)
	option_default_dict.update({
						('inputDir', 0, ): ['', 'L', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('minDepth', 0, float): [0, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						
						('intervalOverlapSize', 1, int): [300000, 'U', 1, 'overlap #bps/#loci between adjacent intervals from one contig/chromosome,\
				only used for TrioCaller, not for SAMtools/GATK', ],\
						('intervalSize', 1, int): [5000000, 'Z', 1, '#bps/#loci for adjacent intervals from one contig/chromosome (alignment or VCF)', ],\
						("ligateVcfPerlPath", 1, ): ["%s/bin/umake/scripts/ligateVcf.pl", '', 1, 'path to ligateVcf.pl'],\
						})
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractAlignmentWorkflow.__init__(self, **keywords)
	
	registerAllInputFiles = AbstractVCFWorkflow.registerAllInputFiles
	
	def mapReduceOneAlignment(self, workflow=None, alignmentData=None, passingData=None, \
						chrIDSet=None, chr2IntervalDataLs=None, chr2VCFFile=None, \
						outputDirPrefix=None, transferOutput=False, **keywords):
		"""
		2013.1.25
		"""
		mapEachChromosomeDataLs = passingData.mapEachChromosomeDataLs
		mapEachChromosomeDataLs = []
		reduceBeforeEachAlignmentData = passingData.reduceBeforeEachAlignmentData
		mapEachAlignmentData = passingData.mapEachAlignmentData
		preReduceReturnData = passingData.preReduceReturnData
		
		for chromosome in chrIDSet:
			intervalDataLs = chr2IntervalDataLs.get(chromosome)
			VCFFile = chr2VCFFile.get(chromosome)
			if VCFFile is None:
				if self.report:
					sys.stderr.write("WARNING: no VCFFile for chromosome %s. no base-quality recalibration. only local realignment.\n"%\
									(chromosome))
				#continue
			mapEachChromosomeData = self.mapEachChromosome(workflow=workflow, alignmentData=alignmentData, chromosome=chromosome, \
								VCFFile=VCFFile, passingData=passingData, reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
								mapEachAlignmentData=mapEachAlignmentData,\
								transferOutput=False, **keywords)
			passingData.mapEachChromosomeData = mapEachChromosomeData
			mapEachChromosomeDataLs.append(mapEachChromosomeData)
			
			mapEachIntervalDataLs = passingData.mapEachIntervalDataLs
			mapEachIntervalDataLs = []
			
			for intervalData in intervalDataLs:
				if intervalData.file:
					mpileupInterval = intervalData.interval
					bcftoolsInterval = intervalData.file
				else:
					mpileupInterval = intervalData.interval
					bcftoolsInterval = intervalData.interval
				intervalFnameSignature = intervalData.intervalFnameSignature
				overlapInterval = intervalData.overlapInterval
				overlapFilenameSignature = intervalData.overlapIntervalFnameSignature
				
				mapEachIntervalData = self.mapEachInterval(workflow=workflow, alignmentData=alignmentData, intervalData=intervalData,\
									VCFFile=VCFFile, passingData=passingData, reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
									mapEachAlignmentData=mapEachAlignmentData,\
									mapEachChromosomeData=mapEachChromosomeData, transferOutput=False, **keywords)
				passingData.mapEachIntervalData = mapEachIntervalData
				mapEachIntervalDataLs.append(mapEachIntervalData)
				
				linkMapToReduceData = self.linkMapToReduce(workflow=workflow, mapEachIntervalData=mapEachIntervalData, \
									preReduceReturnData=preReduceReturnData, \
									reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
									mapEachAlignmentData=mapEachAlignmentData,\
									passingData=passingData, \
									**keywords)
			
			reduceAfterEachChromosomeData = self.reduceAfterEachChromosome(workflow=workflow, chromosome=chromosome, \
								passingData=passingData, \
								mapEachIntervalDataLs=passingData.mapEachIntervalDataLs,\
								transferOutput=False, data_dir=self.data_dir, \
								**keywords)
			passingData.reduceAfterEachChromosomeData = reduceAfterEachChromosomeData
			passingData.reduceAfterEachChromosomeDataLs.append(reduceAfterEachChromosomeData)
			
			gzipReduceAfterEachChromosomeData = self.addGzipSubWorkflow(workflow=workflow, \
					inputData=reduceAfterEachChromosomeData, transferOutput=transferOutput,\
					outputDirPrefix="%sreduceAfterEachChromosome"%(outputDirPrefix), \
					topOutputDirJob=passingData.gzipReduceAfterEachChromosomeFolderJob, report=False)
			passingData.gzipReduceAfterEachChromosomeFolderJob = gzipReduceAfterEachChromosomeData.topOutputDirJob
	
	def setup(self, inputVCFData=None, chr2IntervalDataLs=None, **keywords):
		"""
		2013.04.01 derive chr2VCFFile only when inputVCFData is available
		2013.1.25
		"""
		#2012.8.26 so that each recalibration will pick up the right vcf
		chr2VCFFile = {}
		if inputVCFData:
			for jobData in inputVCFData.jobDataLs:
				inputF = jobData.file
				chromosome = self.getChrFromFname(os.path.basename(inputF.name))
				chr2VCFFile[chromosome] = inputF
		chrIDSet = set(chr2VCFFile.keys())&set(chr2IntervalDataLs.keys())
		return PassingData(chrIDSet=chrIDSet, chr2VCFFile=chr2VCFFile)
	
	def registerCustomExecutables(self, workflow=None):
		
		"""
		"""
		AbstractAlignmentWorkflow.registerCustomExecutables(self, workflow=workflow)
		AbstractVCFWorkflow.registerCustomExecutables(self, workflow=workflow)
		
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
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if not self.data_dir:
			self.data_dir = self.db.data_dir
		
		if not self.local_data_dir:
			self.local_data_dir = self.db.data_dir
		
		chrLs = self.chr2size.keys()
		chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=self.chr2size, \
													intervalSize=self.intervalSize, \
													intervalOverlapSize=self.intervalOverlapSize)
		
		alignmentLs = self.db.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, data_dir=self.local_data_dir,\
										individual_sequence_file_raw_id_type=self.individual_sequence_file_raw_id_type,\
										country_id_ls=self.country_id_ls, tax_id_ls=self.tax_id_ls)
		alignmentLs = self.db.filterAlignments(alignmentLs, min_coverage=self.sequence_min_coverage,\
						max_coverage=self.sequence_max_coverage, sequence_filtered=self.sequence_filtered, \
						individual_site_id_set=set(self.site_id_ls),\
						mask_genotype_method_id=None, parent_individual_alignment_id=None,\
						country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls))
		
		workflow = self.initiateWorkflow()
		
		registerReferenceData = self.getReferenceSequence(workflow)
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables()
		alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, data_dir=self.data_dir)
		
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName=self.pegasusFolderName)
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("No VCF files in this folder , %s.\n"%self.inputDir)
			sys.exit(0)
		#adding inputVCFData=... is the key difference from the parent class
		self.addAllJobs(workflow=workflow, inputVCFData=inputData, alignmentDataLs=alignmentDataLs, \
					chr2IntervalDataLs=chr2IntervalDataLs, samtools=workflow.samtools, \
				GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
				MergeSamFilesJar=workflow.MergeSamFilesJar, \
				CreateSequenceDictionaryJava=workflow.CreateSequenceDictionaryJava, CreateSequenceDictionaryJar=workflow.CreateSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar,\
				mv=workflow.mv, \
				registerReferenceData=registerReferenceData,\
				needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
				data_dir=self.data_dir, no_of_gatk_threads = 1, transferOutput=True,\
				outputDirPrefix=self.pegasusFolderName)
		
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)

if __name__ == '__main__':
	main_class = AbstractAlignmentAndVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()