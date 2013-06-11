#!/usr/bin/env python
"""
2013.1.25 an abstract class for pegasus workflows that work on alignment files (already aligned).
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, utils
from AbstractNGSWorkflow import AbstractNGSWorkflow

class AbstractAlignmentWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractNGSWorkflow.option_default_dict)
	#option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('intervalOverlapSize', 1, int): [300000, 'U', 1, 'overlap #bps/#loci between adjacent intervals from one contig/chromosome,\
				only used for TrioCaller, not for SAMtools/GATK', ],\
						('intervalSize', 1, int): [5000000, 'Z', 1, '#bps/#loci for adjacent intervals from one contig/chromosome (alignment or VCF)', ],\
						})
	commonAlignmentWorkflowOptionDict = {
						('ind_seq_id_ls', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['', '', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						('alignment_outdated_index', 0, int): [0, '', 1, 'filter based on value of IndividualAlignment.outdated_index.', ],\
						("alignment_method_id", 0, int): [None, 'G', 1, 'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
						("local_realigned", 0, int): [None, '', 1, 'To filter which input alignments to fetch from db (i.e. AlignmentReadBaseQualityRecalibrationWorkflow.py)\
	OR to instruct whether local_realigned should be applied (i.e. ShortRead2AlignmentWorkflow.py)'],\
						('defaultSampleAlignmentDepth', 1, int): [10, '', 1, "when database doesn't have median_depth info for one alignment, use this number instead.", ],\
						('individual_sequence_file_raw_id_type', 1, int): [1, '', 1, "1: only all-library-fused libraries,\n\
		2: only library-specific alignments,\n\
		3: both all-library-fused and library-specific alignments", ],\
						}
	option_default_dict.update(commonAlignmentWorkflowOptionDict)
	partitionWorkflowOptionDict= {
						("selectedRegionFname", 0, ): ["", 'R', 1, 'the file is in bed format, tab-delimited, chr start stop.\
		used to restrict SAMtools/GATK to only make calls at this region. \
		start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.\
		This overrides the contig/chromosome selection approach defined by --contigMaxRankBySize and --contigMinRankBySize. \
		This file would be split into maxNoOfRegionsPerJob lines.'],\
						('maxNoOfRegionsPerJob', 1, int): [5000, 'K', 1, 'Given selectedRegionFname, this dictates the maximum number of regions each job would handle,\
		The actual number could be lower because the regions are first grouped into chromosomes. If one chromosome has <maxNoOfRegionsPerJob, then that job handles less.', ],\
						}
	option_default_dict.update(partitionWorkflowOptionDict)
	
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)	#extra__init__() will be executed inside __init__()
		
	def extra__init__(self):
		"""
		2013.2.14
		"""
		AbstractNGSWorkflow.extra__init__(self)
		
		listArgumentName_data_type_ls = [('ind_seq_id_ls', int), ("ind_aln_id_ls", int)]
		listArgumentName2hasContent = self.processListArguments(listArgumentName_data_type_ls, emptyContent=[])
	
	def addAlignmentAsInputToJobLs(self, workflow=None, alignmentDataLs=None, jobLs=[], jobInputOption=""):
		"""
		2012.1.9
			used in addGenotypeCallJobs() to add alignment files as input to calling jobs
		"""
		if workflow is None:
			workflow = self
			
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			parentJobLs = alignmentData.jobLs
			bamF = alignmentData.bamF
			baiF = alignmentData.baiF
			for job in jobLs:
				if jobInputOption:
					job.addArguments(jobInputOption)
				job.addArguments(bamF)
				#it's either symlink or stage-in
				job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
				job.uses(baiF, transfer=True, register=True, link=Link.INPUT)
				for parentJob in parentJobLs:
					if parentJob:
						self.addJobDependency(workflow=workflow, parentJob=parentJob, childJob=job)
		
	def addAlignmentAsInputToPlatypusJobLs(self, workflow=None, alignmentDataLs=None, jobLs=[], jobInputOption="--bamFiles"):
		"""
		2013.05.21 bugfix: pegasus/condor would truncate long single-argument.
		2013.05.16 different from addAlignmentAsInputToJobLs, that platypus is like this:
				--bamFiles=1.bam,2.bam,3.bam
			used in addGenotypeCallJobs() to add alignment files as input to calling jobs
		"""
		if workflow is None:
			workflow = self
		for job in jobLs:
			if jobInputOption:
				job.addArguments(jobInputOption)
			fileArgumentLs = []
			alignmentFileFolder = None
			for alignmentData in alignmentDataLs:
				alignment = alignmentData.alignment
				parentJobLs = alignmentData.jobLs
				bamF = alignmentData.bamF
				baiF = alignmentData.baiF
				
				fileArgumentLs.append(bamF.name)
				if alignmentFileFolder is None:
					alignmentFileFolder = os.path.split(bamF.name)[0]
				#it's either symlink or stage-in
				job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
				job.uses(baiF, transfer=True, register=True, link=Link.INPUT)
				for parentJob in parentJobLs:
					if parentJob:
						self.addJobDependency(workflow=workflow, parentJob=parentJob, childJob=job)
			#if alignmentFileFolder:	#2013.05.21 pegasus/condor would truncate long single-argument.
			#	job.addArguments('%s/*.bam'%(alignmentFileFolder))
			#else:
			job.addArguments(','.join(fileArgumentLs))
	
	def addAddRG2BamJobsAsNeeded(self, workflow=None, alignmentDataLs=None, site_handler=None, input_site_handler=None, \
							addOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None, \
							BuildBamIndexFilesJava=None, BuildBamIndexJar=None, \
							mv=None, namespace="workflow", version="1.0", \
							data_dir=None, tmpDir="/tmp"):
		"""
		2012.4.5
			fix some bugs here
		2011-9-15
			add a read group only when the alignment doesn't have it according to db record
			DBVervet.pokeBamReadGroupPresence() from misc.py helps to fill in db records if it's unclear.
		2011-9-14
			The read-group adding jobs will have a "move" part that overwrites the original bam&bai if site_handler and input_site_handler is same.
			For those alignment files that don't need to. It doesn't matter. pegasus will transfer/symlink them.
		"""
		sys.stderr.write("Adding add-read-group2BAM jobs for %s alignments if read group is not detected ..."%(len(alignmentDataLs)))
		if workflow is None:
			workflow = self
		job_max_memory = 3500	#in MB
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		indexJobMaxMem=2500
		
		addRG2BamDir = None
		addRG2BamDirJob = None
		
		no_of_rg_jobs = 0
		returnData = []
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			parentJobLs = alignmentData.jobLs
			bamF = alignmentData.bamF
			baiF = alignmentData.baiF
			if alignment.read_group_added!=1:
				if addRG2BamDir is None:
					addRG2BamDir = "addRG2Bam"
					addRG2BamDirJob = self.addMkDirJob(outputDir=addRG2BamDir)
				
				# add RG to this bam
				sequencer = alignment.individual_sequence.sequencer
				#read_group = '%s_%s_%s_%s_vs_%s'%(alignment.id, alignment.ind_seq_id, alignment.individual_sequence.individual.code, \
				#						sequencer, alignment.ref_ind_seq_id)
				read_group = alignment.getReadGroup()	##2011-11-02
				if sequencer=='454':
					platform_id = 'LS454'
				elif sequencer=='GA':
					platform_id = 'ILLUMINA'
				else:
					platform_id = 'ILLUMINA'
				
				# the add-read-group job
				#addRGJob = Job(namespace=namespace, name=addRGExecutable.name, version=version)
				addRGJob = Job(namespace=namespace, name=addOrReplaceReadGroupsJava.name, version=version)
				outputRGSAM = File(os.path.join(addRG2BamDir, os.path.basename(alignment.path)))
				
				addRGJob.addArguments(javaMemRequirement, '-jar', AddOrReplaceReadGroupsJar, \
									"INPUT=", bamF,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM, 'SORT_ORDER=coordinate', "VALIDATION_STRINGENCY=LENIENT")
									#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexJar won't fail.)
				self.addJobUse(addRGJob, file=AddOrReplaceReadGroupsJar, transfer=True, register=True, link=Link.INPUT)
				if tmpDir:
					addRGJob.addArguments("TMP_DIR=%s"%tmpDir)
				addRGJob.uses(bamF, transfer=True, register=True, link=Link.INPUT)
				addRGJob.uses(baiF, transfer=True, register=True, link=Link.INPUT)
				addRGJob.uses(outputRGSAM, transfer=True, register=True, link=Link.OUTPUT)
				yh_pegasus.setJobProperRequirement(addRGJob, job_max_memory=job_max_memory)
				for parentJob in parentJobLs:
					if parentJob:
						workflow.depends(parent=parentJob, child=addRGJob)
				workflow.addJob(addRGJob)
				
				
				index_sam_job = self.addBAMIndexJob(workflow, BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar, \
					inputBamF=outputRGSAM, parentJobLs=[addRGJob], transferOutput=True, javaMaxMemory=2000)
				newAlignmentData = PassingData(alignment=alignment)
				newAlignmentData.jobLs = [index_sam_job, addRGJob]
				newAlignmentData.bamF = index_sam_job.bamFile
				newAlignmentData.baiF = index_sam_job.baiFile
				"""
				# add the index job to the bamF (needs to be re-indexed)
				index_sam_job = Job(namespace=namespace, name=BuildBamIndexFilesJava.name, version=version)
				
				if input_site_handler==site_handler:	#on the same site. overwrite the original file without RG
					mvJob = Job(namespace=namespace, name=mv.name, version=version)
					mvJob.addArguments(outputRGSAM, inputFname)	#watch, it's inputFname, not input. input is in relative path.
					#samToBamJob.uses(outputRG, transfer=False, register=True, link=Link.OUTPUT)	#don't register it here
					workflow.addJob(mvJob)
					workflow.depends(parent=addRGJob, child=mvJob)
					bai_output = File('%s.bai'%inputFname)	#in absolute path, don't register it to the job
				else:
					##on different site, input for index should be outputRGSAM and register it as well
					mvJob = addRGJob
					bamF = outputRGSAM	
					addRGJob.uses(outputRGSAM, transfer=True, register=True, link=Link.OUTPUT)
					bai_output = File('%s.bai'%outputRGSAMFname)
					index_sam_job.uses(bai_output, transfer=True, register=False, link=Link.OUTPUT)
				
				index_sam_job
				
				yh_pegasus.setJobProperRequirement(index_sam_job, job_max_memory=indexJobMaxMem)
				
				workflow.addJob(index_sam_job)
				workflow.depends(parent=mvJob, child=index_sam_job)
				alignmentId2RGJobDataLs[alignment.id]= [index_sam_job, inputFile, bai_output]
				"""
				no_of_rg_jobs += 1
			else:
				newAlignmentData = alignmentData
			returnData.append(newAlignmentData)
		sys.stderr.write(" %s alignments need read-group addition. Done\n"%(no_of_rg_jobs))
		return returnData
	
	def preReduce(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
			setup additional mkdir folder jobs, before mapEachAlignment, mapEachChromosome, mapReduceOneAlignment 
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachChromosome(self, workflow=None, alignmentData=None, chromosome=None,\
				VCFJobData=None, passingData=None, reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def map(self, workflow=None, alignmentData=None, intervalData=None,\
		VCFJobData=None, passingData=None, mapEachChromosomeData=None, transferOutput=True, **keywords):
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

	def mapEachAlignment(self, workflow=None, alignmentData=None,  passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.22
			similar to reduceBeforeEachAlignmentData() but for mapping programs that run on one alignment each.
			
			passingData.AlignmentJobAndOutputLs = []
			passingData.bamFnamePrefix = bamFnamePrefix
			passingData.individual_alignment = alignment
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		refFastaF = passingData.refFastaFList[0]
		
		alignment = alignmentData.alignment
		parentJobLs = alignmentData.jobLs
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		
		bamFnamePrefix = alignment.getReadGroup()
		
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
	
	def mapReduceOneAlignment(self, workflow=None, alignmentData=None, passingData=None, \
						chrIDSet=None, chrSizeIDList=None, chr2IntervalDataLs=None, chr2VCFJobData=None, \
						outputDirPrefix=None, transferOutput=False, skipChromosomeIfVCFMissing=False, **keywords):
		"""
		2013.04.11 moved from AbstractAlignmentAndVCFWorkflow.py
		2013.04.08, added skipChromosomeIfVCFMissing
		2013.1.25
		"""
		returnData = PassingData()
		mapEachChromosomeDataLs = passingData.mapEachChromosomeDataLs
		mapEachChromosomeDataLs = []
		reduceBeforeEachAlignmentData = passingData.reduceBeforeEachAlignmentData
		mapEachAlignmentData = passingData.mapEachAlignmentData
		preReduceReturnData = passingData.preReduceReturnData
		
		for chromosomeSize, chromosome in chrSizeIDList:
			intervalDataLs = chr2IntervalDataLs.get(chromosome)
			if chr2VCFJobData:
				VCFJobData = chr2VCFJobData.get(chromosome)
			else:
				VCFJobData = None
			if VCFJobData is None:
				if self.report:
					sys.stderr.write("WARNING: no VCFJobData for chromosome %s.\n"%(chromosome))
				if skipChromosomeIfVCFMissing:
					continue
				VCFJobData = PassingData(job=None, jobLs=[],\
									vcfFile=None, tbi_F=None, file=None, fileLs=[])
				VCFFile = None
			else:
				VCFFile = VCFJobData.file
				if VCFFile is None:
					if self.report:
						sys.stderr.write("WARNING: no VCFFile for chromosome %s.\n"%(chromosome))
					if skipChromosomeIfVCFMissing:
						continue
			passingData.chromosome = chromosome	#2013.04.08
			mapEachChromosomeData = self.mapEachChromosome(workflow=workflow, alignmentData=alignmentData, chromosome=chromosome, \
								VCFJobData=VCFJobData, passingData=passingData, reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
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
									chromosome=chromosome,\
									VCFJobData=VCFJobData, passingData=passingData, reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
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
		return returnData
	
	def setup(self, chr2IntervalDataLs=None, **keywords):
		"""
		2013.04.09 added chrSizeIDList in return
		2013.1.25
			chr2VCFJobData is None.
		"""
		chrIDSet = set(chr2IntervalDataLs.keys())
		chrSizeIDList = [(chromosomeSize, chromosome) for chromosome, chromosomeSize in self.chr2size.iteritems()]
		chrSizeIDList.sort()
		chrSizeIDList.reverse()	#from big to small
		return PassingData(chrIDSet=chrIDSet, chr2VCFJobData=None, chrSizeIDList=chrSizeIDList)
	
	def addAllJobs(self, workflow=None, alignmentDataLs=None, chr2IntervalDataLs=None, samtools=None, \
				GenomeAnalysisTKJar=None, \
				MergeSamFilesJar=None, \
				CreateSequenceDictionaryJava=None, CreateSequenceDictionaryJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexJar=None,\
				mv=None, skipDoneAlignment=False,\
				registerReferenceData=None, \
				needFastaIndexJob=False, needFastaDictJob=False, \
				data_dir=None, no_of_gatk_threads = 1, \
				outputDirPrefix="", transferOutput=True, **keywords):
		"""
		2012.7.26
		"""
		prePreprocessData = self.setup(chr2IntervalDataLs=chr2IntervalDataLs, **keywords)
		chrIDSet = prePreprocessData.chrIDSet
		chrSizeIDList = prePreprocessData.chrSizeIDList
		chr2VCFJobData = prePreprocessData.chr2VCFJobData
		
		sys.stderr.write("Adding jobs that work on %s alignments (& possibly VCFs) for %s chromosomes/contigs ..."%\
						(len(alignmentDataLs), len(chrIDSet)))
		refFastaFList = registerReferenceData.refFastaFList
		refFastaF = refFastaFList[0]

		topOutputDir = "%sMap"%(outputDirPrefix)
		topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
		
		plotOutputDir = "%sPlot"%(outputDirPrefix)
		plotOutputDirJob = self.addMkDirJob(outputDir=plotOutputDir)
		
		reduceOutputDir = "%sReduce"%(outputDirPrefix)
		reduceOutputDirJob = self.addMkDirJob(outputDir=reduceOutputDir)
		
		if needFastaDictJob or registerReferenceData.needPicardFastaDictJob:
			fastaDictJob = self.addRefFastaDictJob(CreateSequenceDictionaryJava=CreateSequenceDictionaryJava, \
										CreateSequenceDictionaryJar=CreateSequenceDictionaryJar, refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
		else:
			fastaDictJob = None
			refFastaDictF = registerReferenceData.refPicardFastaDictF
		
		if needFastaIndexJob or registerReferenceData.needSAMtoolsFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
		else:
			fastaIndexJob = None
			refFastaIndexF = registerReferenceData.refSAMtoolsFastaIndexF
			
		
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		#2012.9.22 AlignmentJobAndOutputLs is a relic.
		#	but it's similar to mapEachIntervalDataLs but designed for addAlignmentMergeJob(),
		#	so AlignmentJobAndOutputLs gets re-set for every alignment.
		# 	mapEachAlignmentDataLs is never reset.
		#	mapEachChromosomeDataLs is reset right after a new alignment is chosen.
		#	mapEachIntervalDataLs is reset right after each chromosome is chosen.
		#	all reduce dataLs never gets reset.
		passingData = PassingData(AlignmentJobAndOutputLs=[], \
					alignmentDataLs = alignmentDataLs,\
					bamFnamePrefix=None, \
					
					outputDirPrefix=outputDirPrefix, \
					topOutputDirJob=topOutputDirJob,\
					plotOutputDirJob=plotOutputDirJob,\
					reduceOutputDirJob = reduceOutputDirJob,\
					
					refFastaFList=refFastaFList, \
					registerReferenceData= registerReferenceData,\
					refFastaF=refFastaFList[0],\
					
					fastaDictJob = fastaDictJob,\
					refFastaDictF = refFastaDictF,\
					fastaIndexJob = fastaIndexJob,\
					refFastaIndexF = refFastaIndexF,\
					
					chromosome=None,\
					chrIDSet=chrIDSet,\
					chrSizeIDList = chrSizeIDList,\
					chr2IntervalDataLs=chr2IntervalDataLs,\
					
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
		no_of_alignments_worked_on= 0
		for alignmentData in passingData.alignmentDataLs:
			alignment = alignmentData.alignment
			parentJobLs = alignmentData.jobLs
			bamF = alignmentData.bamF
			baiF = alignmentData.baiF
			
			bamFnamePrefix = alignment.getReadGroup()
			
			passingData.AlignmentJobAndOutputLs = []
			passingData.bamFnamePrefix = bamFnamePrefix
			passingData.individual_alignment = alignment
			passingData.alignmentData = alignmentData
			
			if skipDoneAlignment and self.isThisAlignmentComplete(individual_alignment=alignment, data_dir=data_dir):
				continue
			no_of_alignments_worked_on += 1
			mapEachAlignmentData = self.mapEachAlignment(workflow=workflow, alignmentData=alignmentData, passingData=passingData, \
								transferOutput=False, \
								preReduceReturnData=preReduceReturnData, **keywords)
			passingData.mapEachAlignmentDataLs.append(mapEachAlignmentData)
			passingData.mapEachAlignmentData = mapEachAlignmentData
			
			reduceBeforeEachAlignmentData = self.reduceBeforeEachAlignment(workflow=workflow, passingData=passingData, \
													preReduceReturnData=preReduceReturnData, transferOutput=False, \
													**keywords)
			passingData.reduceBeforeEachAlignmentData = reduceBeforeEachAlignmentData
			passingData.reduceBeforeEachAlignmentDataLs.append(reduceBeforeEachAlignmentData)
			
			
			mapReduceOneAlignmentReturnData = self.mapReduceOneAlignment(workflow=workflow, alignmentData=alignmentData, \
							passingData=passingData, \
							chrIDSet=chrIDSet, chrSizeIDList=chrSizeIDList, \
							chr2IntervalDataLs=chr2IntervalDataLs, chr2VCFJobData=chr2VCFJobData, \
							outputDirPrefix=outputDirPrefix, transferOutput=transferOutput)
			
			reduceAfterEachAlignmentData = self.reduceAfterEachAlignment(workflow=workflow, \
												mapEachAlignmentData=mapEachAlignmentData,\
												mapEachChromosomeDataLs=passingData.mapEachChromosomeDataLs,\
												reduceAfterEachChromosomeDataLs=passingData.reduceAfterEachChromosomeDataLs,\
												passingData=passingData, \
												transferOutput=False, data_dir=data_dir, **keywords)
			passingData.reduceAfterEachAlignmentData = reduceAfterEachAlignmentData
			passingData.reduceAfterEachAlignmentDataLs.append(reduceAfterEachAlignmentData)
			
			gzipReduceBeforeEachAlignmentData = self.addGzipSubWorkflow(workflow=workflow, \
						inputData=reduceBeforeEachAlignmentData, transferOutput=transferOutput,\
						outputDirPrefix="%sReduceBeforeEachAlignment"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceBeforeEachAlignmentFolderJob, report=False)
			passingData.gzipReduceBeforeEachAlignmentFolderJob = gzipReduceBeforeEachAlignmentData.topOutputDirJob
			
			gzipReduceAfterEachAlignmentData = self.addGzipSubWorkflow(workflow=workflow, \
						inputData=reduceAfterEachAlignmentData, transferOutput=transferOutput,\
						outputDirPrefix="%sReduceAfterEachAlignment"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceAfterEachAlignmentFolderJob, \
						report=False)
			passingData.gzipReduceAfterEachAlignmentFolderJob = gzipReduceAfterEachAlignmentData.topOutputDirJob
		reduceReturnData = self.reduce(workflow=workflow, passingData=passingData, \
							mapEachAlignmentData=passingData.mapEachAlignmentData, \
							reduceAfterEachAlignmentDataLs=passingData.reduceAfterEachAlignmentDataLs,\
							**keywords)
		passingData.reduceReturnData = reduceReturnData
		
		
		#2012.9.18 gzip the final output
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=preReduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%sGzipPreReduce"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipPreReduceFolderJob, \
						report=False)
		passingData.gzipPreReduceFolderJob = newReturnData.topOutputDirJob
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=reduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%sGzipReduce"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceFolderJob, \
						report=False)
		passingData.gzipReduceFolderJob = newReturnData.topOutputDirJob
		
		sys.stderr.write("%s alignments to be worked on. %s jobs.\n"%(no_of_alignments_worked_on, self.no_of_jobs))
		return returnData
	
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
		returnData = []
		for alignment in alignmentLs:
			inputFname = os.path.join(data_dir, alignment.path)
			input = File(alignment.path)	#relative path, induces symlinking or stage-in
			baiFilepath = '%s.bai'%(inputFname)
			if checkFileExistence and (not os.path.isfile(inputFname) or not os.path.isfile(baiFilepath)):
				if not os.path.isfile(inputFname):
					sys.stderr.write("registerAlignmentAndItsIndexFile() Warning: %s does not exist. skip entire alignment.\n"%(inputFname))
				if not os.path.isfile(baiFilepath):
					sys.stderr.write("registerAlignmentAndItsIndexFile() Warning: %s does not exist. skip entire alignment.\n"%(baiFilepath))
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
		AbstractNGSWorkflow.registerCustomExecutables(self, workflow=workflow)
		
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
		#self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name="exampleJava", clusterSizeMultipler=0.3)
	
	def setup_run(self):
		"""
		2013.04.07 wrap all standard pre-run() related functions into this function.
			setting up for run(), called by run()
		"""
		pdata = AbstractNGSWorkflow.setup_run(self)
		workflow = pdata.workflow
		
		chrLs = self.chr2size.keys()
		chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=self.chr2size, \
													intervalSize=self.intervalSize, \
													intervalOverlapSize=self.intervalOverlapSize)
		
		alignmentLs = self.getAlignments()
		
		registerReferenceData = self.getReferenceSequence()
		
		alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow=workflow, alignmentLs=alignmentLs, data_dir=self.data_dir)
		pdata.alignmentLs = alignmentLs
		pdata.alignmentDataLs = alignmentDataLs
		pdata.chr2IntervalDataLs = chr2IntervalDataLs
		pdata.registerReferenceData = registerReferenceData
		return pdata
	
	def run(self):
		"""
		2013.1.25
		"""
		
		pdata = self.setup_run()
		workflow = pdata.workflow
		
		self.addAllJobs(workflow=workflow, alignmentDataLs=pdata.alignmentDataLs, \
				chr2IntervalDataLs=pdata.chr2IntervalDataLs, samtools=workflow.samtools, \
				GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
				MergeSamFilesJar=workflow.MergeSamFilesJar, \
				CreateSequenceDictionaryJava=workflow.CreateSequenceDictionaryJava, \
				CreateSequenceDictionaryJar=workflow.CreateSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar,\
				mv=workflow.mv, skipDoneAlignment=self.skipDoneAlignment,\
				registerReferenceData=pdata.registerReferenceData,\
				needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
				data_dir=self.data_dir, no_of_gatk_threads = 1, transferOutput=True,\
				outputDirPrefix=self.pegasusFolderName)
		
		self.end_run()

if __name__ == '__main__':
	main_class = AbstractAlignmentWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()