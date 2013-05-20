#!/usr/bin/env python
"""
Examples:
	#2013.04.12 test on one outdated alignment (--local_realigned 0 --alignment_outdated_index 1)
	%s --ind_aln_id_ls 552 --ref_ind_seq_id 524  --ref_genome_outdated_index 1 -o dags/ReduceReads/ReduceReadsAln552.xml
		-l condorpool -j condorpool -z uclaOffice -u yh --intervalSize 20000000
		--intervalOverlapSize 0 --contigMaxRankBySize 1000  --clusters_size 5
		-J ~/bin/jdk/bin/java --commit --skipDoneAlignment --local_realigned 0 --alignment_outdated_index 1
	
	# 2013.04.12 use sequence coverage to filter alignments
	%s 	--sequence_min_coverage 0 --sequence_max_coverage 2  --ind_seq_id_ls 632-3230  --ind_aln_id_ls 32
		--ref_ind_seq_id 3280 -o dags/ReduceReads/ReduceReads_ISQ632_3230_coverage0_2.xml
		-l hcondor -j hcondor -z localhost -u yh --intervalSize 20000000 --intervalOverlapSize 0
		-e /u/home/eeskin/polyacti --contigMaxRankBySize 250
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--clusters_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java
		--commit --skipDoneAlignment
		# --ref_genome_version 2 #(optional, as by default, it gets the outdated_index=0 reference chromosomes from GenomeDB)
		# --ref_genome_outdated_index 1 #to get old reference. incompatible here as alignment is based on 3280, new ref.
	
Description:
	#2013.04.11 workflow that carries out ReduceReads for all alignments in map-reduce fashion.
		preferably on local_realigned=1 alignments. It is set to 1 by default.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, \
	figureOutDelimiter, getColName2IndexFromHeader, utils
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule import VCFFile
from vervet.src import AbstractVervetAlignmentWorkflow

parentClass = AbstractVervetAlignmentWorkflow
class AlignmentReduceReadsWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(parentClass.partitionWorkflowOptionDict.copy())
	option_default_dict.update({
							})
	option_default_dict[('intervalSize', 1, int)][0] = 20000000
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('local_realigned', 0, int)][0] = 1
	option_default_dict[('completedAlignment', 0, int)][0]=1	#2013.05.03
	
	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
		self.chr2IndelVCFJobData = None	#2013.04.04 mark this variable. setup in setup()
		self.candidateCountCovariatesJob = None	#2013.04.09 this BQSR count-variates job encompasses one of the top big intervals.
			# replacing equivalent jobs for small intervals (not accurate if intervals are too small)
		#AlignmentToCallPipeline.__init__(self, **keywords)
		#self.inputDir = os.path.abspath(self.inputDir)
	
	def mapEachInterval(self, workflow=None, alignmentData=None, intervalData=None, chromosome=None, \
							VCFJobData=None, passingData=None, reduceBeforeEachAlignmentData=None,\
							mapEachChromosomeData=None, transferOutput=False, \
							**keywords):
		"""
		2013.03.31 use VCFJobData to decide whether to add BQSR jobs, called in ShortRead2AlignmentWorkflow.py
		2012.9.17
		"""
		if workflow is None:
			workflow = self
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		
		alignment = alignmentData.alignment
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		bamFnamePrefix = passingData.bamFnamePrefix
		
		#SNPVCFFile = VCFJobData.file
		#if SNPVCFFile is None or VCFJobData is None:	#2013.04.09	BQSR requires a VCF input regardless of the chromosome
		#	VCFJobData = self.randomSNPVCFJobDataForBQSR
		
		#SNPVCFFile = VCFJobData.file
		#SNPVCFJobLs = VCFJobData.jobLs
		
		if intervalData.file:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.file
		else:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.interval
		intervalFnameSignature = intervalData.intervalFnameSignature
		overlapInterval = intervalData.overlapInterval
		overlapFilenameSignature = intervalData.overlapIntervalFnameSignature
		span = intervalData.span
		
		if chromosome is None:
			chromosome = getattr(passingData, 'chromosome', None)
		
		
		median_depth = getattr(alignment, 'median_depth', 4)
		readSpace = median_depth * span
		#base is 4X coverage in 20Mb region => 120 minutes
		reduceReadsJobWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=readSpace, \
							baseInputVolume=4*20000000, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value
		#base is 4X, => 5000M
		reduceReadsJobMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=median_depth, \
							baseInputVolume=4, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=8000).value
							
		reduceReadsBamFile = File(os.path.join(topOutputDirJob.output, '%s_%s.reduceReads.bam'%\
											(bamFnamePrefix, overlapFilenameSignature)))
		#Default downsampling setting is 40 in GATK 2.4.9
		extraArgumentList= ["--downsample_to_coverage 60", "--downsampling_type BY_SAMPLE"]
		reduceReadsJob = self.addGATKJob(executable=self.ReduceReadsJava, GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
					GATKAnalysisType='ReduceReads',\
					inputFile=bamF, inputArgumentOption="-I", refFastaFList=passingData.refFastaFList, inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					interval=overlapInterval, outputFile=reduceReadsBamFile, \
					parentJobLs=alignmentData.jobLs, transferOutput=False, \
					job_max_memory=reduceReadsJobMaxMemory,\
					frontArgumentList=None, extraArguments=None, \
					extraArgumentList=extraArgumentList, \
					extraOutputLs=[], \
					extraDependentInputLs=[baiF], no_of_cpus=None, \
					walltime=reduceReadsJobWalltime)
		indexBamJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
										BuildBamIndexJar=self.BuildBamIndexJar, \
					inputBamF=reduceReadsJob.output,\
					parentJobLs=[reduceReadsJob], \
					transferOutput=False, job_max_memory=3000, \
					walltime=max(120, int(reduceReadsJobWalltime/3)))
		passingData.AlignmentJobAndOutputLs.append(PassingData(parentJobLs=[reduceReadsJob, indexBamJob], \
															file=reduceReadsJob.output))
		return returnData
	
	def reduceAfterEachAlignment(self, workflow=None, passingData=None, transferOutput=False, data_dir=None, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		if workflow is None:
			workflow = self
		AlignmentJobAndOutputLs = getattr(passingData, 'AlignmentJobAndOutputLs', [])
		bamFnamePrefix = passingData.bamFnamePrefix
		topOutputDirJob = passingData.topOutputDirJob
		individual_alignment = passingData.individual_alignment
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		if len(AlignmentJobAndOutputLs)>0:	#2012.3.29	merge alignment output only when there is something to merge!
			#2013.04.09 create a new child alignment local_realigned =1, etc.
			new_individual_alignment = self.db.copyParentIndividualAlignment(parent_individual_alignment_id=individual_alignment.id,\
										data_dir=self.data_dir, local_realigned=individual_alignment.local_realigned,\
										reduce_reads=1)
			
			#2013.04.09 replace read_group with the new one to each alignment job
			NewAlignmentJobAndOutputLs = []
			for AlignmentJobAndOutput in AlignmentJobAndOutputLs:
				#2012.9.19 add a AddReadGroup job
				alignmentJob, indexAlignmentJob = AlignmentJobAndOutput.parentJobLs
				fname_prefix = os.path.splitext(alignmentJob.output.name)[0]
				outputRGBAM = File("%s.isq_RG.bam"%(fname_prefix))
				addRGJob = self.addReadGroupInsertionJob(workflow=workflow, individual_alignment=new_individual_alignment, \
									inputBamFile=alignmentJob.output, \
									outputBamFile=outputRGBAM,\
									addOrReplaceReadGroupsJava=self.addOrReplaceReadGroupsJava, \
									AddOrReplaceReadGroupsJar=self.AddOrReplaceReadGroupsJar,\
									parentJobLs=[alignmentJob, indexAlignmentJob], extraDependentInputLs=None, \
									extraArguments=None, job_max_memory = 2500, transferOutput=False)
				
				NewAlignmentJobAndOutputLs.append(PassingData(parentJobLs=[addRGJob], file=addRGJob.output))
			#
			mergedBamFile = File(os.path.join(reduceOutputDirJob.output, '%s.merged.bam'%(bamFnamePrefix)))
			alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=NewAlignmentJobAndOutputLs, \
					outputBamFile=mergedBamFile, \
					samtools=workflow.samtools, java=workflow.java, \
					MergeSamFilesJava=workflow.MergeSamFilesJava, MergeSamFilesJar=workflow.MergeSamFilesJar, \
					BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexJar=workflow.BuildBamIndexJar, \
					mv=workflow.mv, parentJobLs=[reduceOutputDirJob], \
					transferOutput=False)
			#2012.9.19 add/copy the alignment file to db-affliated storage
			#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
			logFile = File(os.path.join(reduceOutputDirJob.output, '%s_2db.log'%(bamFnamePrefix)))
			alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
								inputFile=alignmentMergeJob.output, otherInputFileList=[],\
								individual_alignment_id=new_individual_alignment.id, \
								logFile=logFile, data_dir=data_dir, \
								parentJobLs=[alignmentMergeJob, bamIndexJob], \
								extraDependentInputLs=[bamIndexJob.output], \
								extraArguments=None, transferOutput=transferOutput, \
								job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel, commit=True)
			self.no_of_jobs += 1
			returnData.jobDataLs.append(PassingData(jobLs=[alignment2DBJob], file=alignment2DBJob.logFile, \
											fileLs=[alignment2DBJob.logFile]))
		return returnData
	
	
	def isThisAlignmentComplete(self, individual_alignment=None, data_dir=None):
		"""
		2013.05.04
			this is to check whether the new and to-be-generated alignment is completed already or not.
			in contrast to isThisInputAlignmentComplete()
			
			different from usual: reduce_reads =1
		"""
		new_individual_alignment = self.db.copyParentIndividualAlignment(parent_individual_alignment_id=individual_alignment.id,\
										mask_genotype_method_id=individual_alignment.mask_genotype_method_id,\
										data_dir=self.data_dir, local_realigned=individual_alignment.local_realigned,\
										reduce_reads=1)
		return self.db.isThisAlignmentComplete(individual_alignment=new_individual_alignment, data_dir=data_dir)

	def isThisInputAlignmentComplete(self, individual_alignment=None, data_dir=None, returnFalseIfInexitentFile=True, \
									**keywords):
		"""
		2013.05.04 
			 this checks whether
			#. an alignment file exists
			#. file_size not null in db,
			#. median_depth is not null from db
			
		this is used to check whether an input (to be worked on by downstream programs) is completed or not.
			watch returnFalseIfInexitentFile is True (because you need the file for input)
			
		"""
		stockAnswer = self.db.isThisAlignmentComplete(individual_alignment=individual_alignment, data_dir=data_dir,\
													returnFalseIfInexitentFile=returnFalseIfInexitentFile, **keywords)
		if stockAnswer is True:
			if individual_alignment.median_depth is None:
				#change stockAnswer if median_depth is null
				stockAnswer = False
		return stockAnswer	
			
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		self.setOrChangeExecutableClusterSize(executable=workflow.samtools, clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='ReduceReadsJava', \
															clusterSizeMultipler=1)


if __name__ == '__main__':
	main_class = AlignmentReduceReadsWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
