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
		pdata = self.setup_run()
		workflow = pdata.workflow
				
		inputData = self.registerAllInputFiles(inputDir=self.inputDir, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName=self.pegasusFolderName)
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("Error: No VCF files in the input VCF folder %s.\n"%self.inputDir)
			raise
		#adding inputVCFData=... is the key difference from the parent class
		self.addAllJobs(workflow=workflow, inputVCFData=inputData, alignmentDataLs=pdata.alignmentDataLs, \
					chr2IntervalDataLs=pdata.chr2IntervalDataLs, samtools=workflow.samtools, \
				GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
				MergeSamFilesJar=workflow.MergeSamFilesJar, \
				CreateSequenceDictionaryJava=workflow.CreateSequenceDictionaryJava, CreateSequenceDictionaryJar=workflow.CreateSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar,\
				mv=workflow.mv, \
				registerReferenceData=pdata.registerReferenceData,\
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