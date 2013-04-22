#!/usr/bin/env python
"""
Examples:
	# 2012.9.21 run base quality recalibration on VRC alignments (-S 447), individual_sequence_id from 639-642 (--ind_seq_id_ls ...)
	# filtered sequences (-Q 1), alignment method 2 (-G 2)
	# --contigMaxRankBySize 1000 (top 1000 contigs)
	#  --intervalSize 10000000 (10 million bp for each interval) --intervalOverlapSize 30000 (30kb overlap between intervals),
	%s --inputDir ~/NetworkData/vervet/db/genotype_file/method_17/ --ind_seq_id_ls 639-642
		-S 447 -u yh -z localhost --sequence_filtered 1 --alignment_method_id 2
		-a 524 -o dags/BaseQualityRecalibration/BaseQualityRecalibration_VRC447_vsMethod17.xml
		-l hcondor -j hcondor -z localhost -u yh --contigMaxRankBySize 1000 
		--intervalSize 10000000 --intervalOverlapSize 30000
		--indelVCFFolder ...
		-e /u/home/eeskin/polyacti
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 17
		 --commit --skipDoneAlignment

	# 2012.9.18
	%s  --inputDir ~/NetworkData/vervet/db/genotype_file/method_41 --ind_seq_id_ls 633,634,635,636,637,638 
		--ref_ind_seq_id 524
		-o dags/BaseQualityRecalibration/BaseQualityRecalibration_ISQ633_638_vsMethod41.xml -l hcondor
		-j hcondor -z localhost -u yh --intervalSize 10000000 --intervalOverlapSize 30000
		-e /u/home/eeskin/polyacti
		--indelVCFFolder ...
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--clusters_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 41
		 --commit --skipDoneAlignment
	
	# 2013.3.19 use sequence coverage to filter alignments
	%s  --inputDir ~/NetworkData/vervet/db/genotype_file/method_41
		--sequence_min_coverage 0 --sequence_max_coverage 2  --ind_seq_id_ls 632-3230
		--ref_ind_seq_id 3280 -o dags/BaseQualityRecalibration/BaseQualityRecalibration_ISQ632_3230_coverage0_2_vsMethod41.xml
		-l hcondor -j hcondor -z localhost -u yh --intervalSize 10000000 --intervalOverlapSize 30000
		-e /u/home/eeskin/polyacti --contigMaxRankBySize 250
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--clusters_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 41
		--indelVCFFolder ~/NetworkData/vervet/db/genotype_file/method_88 --commit --skipDoneAlignment
		# --ref_genome_version 2 #(optional, as by default, it gets the outdated_index=0 reference chromosomes from GenomeDB)
		# --ref_genome_outdated_index 1 #to get old reference. incompatible here as alignment is based on 3280, new ref.
		# --needFastaDictJob --needFastaIndexJob
	
Description:
	#2013.04.09  recalibration needs ~500K reads to get accurate estimate. So adjust the --intervalSize according to coverage.  
	Input: a VCF folder, list of alignment IDs (or use site_id, ref_ind_seq_id to filter )
	It could run local_realignment with or without indel VCF and BQSR with or without SNPVCF. Without a SNPVCF from particular chromosome, it uses a random chromosome. If the interval-span is too short , it use BQSR result from long-intervals.
	
	add alignment file into db (input: parent alignment ID, alignment file path)
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, \
	figureOutDelimiter, getColName2IndexFromHeader, utils
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule import VCFFile
#from AlignmentToCallPipeline import AlignmentToCallPipeline
#from AbstractVervetWorkflow import AbstractVervetWorkflow
#from vervet.src.pegasus.AbstractVervetAlignmentAndVCFWorkflow import AbstractVervetAlignmentAndVCFWorkflow
from vervet.src import AbstractVervetAlignmentAndVCFWorkflow

parentClass = AbstractVervetAlignmentAndVCFWorkflow
class AlignmentReadBaseQualityRecalibrationWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(parentClass.partitionWorkflowOptionDict.copy())
	option_default_dict.update({
				('mask_genotype_method_id', 0, int):[None, '', 1, 'which genotype method is used to mask out polymorphic sites for recalibration'],\
				('indelVCFFolder', 0, ): [None, '', 1, 'folder that contains in-del vcf or vcf.gz files that will be used for GATK2 indel-realigner,\
	required if argument local_realigned is non-zero', ],\
							})
	option_default_dict[('intervalSize', 1, int)][0] = 10000000
	option_default_dict[('local_realigned', 0, int)][0] = 0
	"""
	option_default_dict.update({
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						('inputDir', 0, ): [None, 'i', 1, 'folder containing vcf or vcf.gz files that will be used by CountCovariates to exclude variant sites', ],\
						("sequence_filtered", 0, int): [None, 'Q', 1, 'To filter alignments. None: whatever; 0: unfiltered sequences, 1: filtered sequences'],\
						("alignment_method_id", 0, int): [None, 'G', 1, 'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
						('intervalOverlapSize', 1, int): [3000, '', 1, 'extension of an interval on each side. overlap size is actually 2X of this number.\
							Its value should be bigger than maximum read length, to cover reads that are partially aligned within one interval.', ],\
						('intervalSize', 1, int): [5000000, '', 1, 'size for adjacent intervals from one contig/chromosome', ],\
						("site_id_ls", 0, ): ["", 'S', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
						('run_type', 1, int): [1, 'y', 1, '', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	"""
	
	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
		self.chr2IndelVCFJobData = None	#2013.04.04 mark this variable. setup in setup()
		self.candidateCountCovariatesJob = None	#2013.04.09 this BQSR count-variates job encompasses one of the top big intervals.
			# replacing equivalent jobs for small intervals (not accurate if intervals are too small)
		#AlignmentToCallPipeline.__init__(self, **keywords)
		#self.inputDir = os.path.abspath(self.inputDir)
	
	def mapEachAlignment(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.22
			similar to reduceBeforeEachAlignmentData() but for mapping programs that run on one alignment each.
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachChromosome(self, workflow=None, alignmentData=None, chromosome=None,\
				VCFJobData=None, passingData=None, reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		if workflow is None:
			workflow = self
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		
		alignment = alignmentData.alignment
		parentJobLs = alignmentData.jobLs
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		bamFnamePrefix = passingData.bamFnamePrefix
		
		
		"""
		#2012.9.21 perhaps a downsampling job
		outputFname = os.path.join(topOutputDirJob.output, '%s_%s.bam'%(bamFnamePrefix, overlapFilenameSignature))
		outputFile = File(outputFname)
		selectAlignmentJob, bamIndexJob1 = self.addSelectAlignmentJob(executable=workflow.samtools, inputFile=bamF, \
				outputFile=outputFile, region=overlapInterval, parentJobLs=[topOutputDirJob] + parentJobLs, \
				extraDependentInputLs=[baiF], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, needBAMIndexJob=True)
		"""
		
		"""
		#2012.9.21 count covariates job is moved to map()
		recalFile = File(os.path.join(topOutputDirJob.output, '%s_%s.recal_data.csv'%(bamFnamePrefix, chromosome)))
		countCovariatesJob = self.addGATKBaseRecalibratorJob(GenomeAnalysisTKJar=workflow.GenomeAnalysisTK2Jar, inputFile=bamF, \
								VCFFile=VCFFile, interval=chromosome, outputFile=recalFile, \
								refFastaFList=passingData.refFastaFList, parentJobLs=[topOutputDirJob]+parentJobLs, 
								extraDependentInputLs=[baiF, VCFFile.tbi_F], \
								transferOutput=False, \
								extraArguments=None, job_max_memory=4000)
		
		self.no_of_jobs += 1
		returnData.countCovariatesJob = countCovariatesJob
		returnData.jobDataLs.append(PassingData(jobLs=[countCovariatesJob], file=countCovariatesJob.recalFile, \
											fileLs=[countCovariatesJob.recalFile]))
		"""
		
		return returnData
	
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
		
		SNPVCFFile = VCFJobData.file
		if SNPVCFFile is None or VCFJobData is None:	#2013.04.09	BQSR requires a VCF input regardless of the chromosome
			VCFJobData = self.randomSNPVCFJobDataForBQSR
		
		SNPVCFFile = VCFJobData.file
		SNPVCFJobLs = VCFJobData.jobLs
		
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
		
		"""
		outputFname = os.path.join(topOutputDirJob.output, '%s_%s.bam'%(bamFnamePrefix, overlapFilenameSignature))
		outputFile = File(outputFname)
		selectAlignmentJob, bamIndexJob1 = self.addSelectAlignmentJob(executable=workflow.samtools, inputFile=bamF, \
				outputFile=outputFile, region=overlapInterval, parentJobLs=[topOutputDirJob] + parentJobLs, \
				extraDependentInputLs=[baiF], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, needBAMIndexJob=True)
		"""
		
		"""
		#2013.3.18 local realignment
		java -Xmx2g -jar GenomeAnalysisTK.jar -I input.bam -R ref.fasta  -T RealignerTargetCreator
			-o forIndelRealigner.intervals [--known /path/to/indels.vcf]
			
		java -Xmx4g -jar GenomeAnalysisTK.jar -I input.bam -R ref.fasta -T IndelRealigner 
			-targetIntervals forIndelRealigner.intervals \
			-o realignedBam.bam \
			[-known /path/to/indels.vcf] \
			[-compress 0]    (this argument recommended to speed up the process *if* this is only a temporary file; otherwise, use the default value)

		"""
		#get the indel VCF file
		indelVCFJobData = self.chr2IndelVCFJobData.get(chromosome)
		if not indelVCFJobData:
			if self.report:
				sys.stderr.write("Warning from mapEachInterval(): no indel VCF for local realignment for chromosome %s. resort to ab-initial local-realignment.\n \
	(not ideal because intervals with too many reads will be skipped.).\n"%(chromosome))
			indelVCFFile = None
			indelVCFFileLs = []
			indelVCFFileJobLs = []
		else:
			indelVCFFile = indelVCFJobData.file
			indelVCFFileLs = indelVCFJobData.fileLs
			indelVCFFileJobLs = indelVCFJobData.jobLs
		
		median_depth = getattr(alignment, 'median_depth', 4)
		readSpace = median_depth * span
		#base is 4X coverage in 10Mb region => 120 minutes
		indelRealignmentJobWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=readSpace, \
							baseInputVolume=4*10000000, baseJobPropertyValue=120, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value
		#base is 4X, => 5000M
		indelRealignmentJobMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=median_depth, \
							baseInputVolume=4, baseJobPropertyValue=5000, \
							minJobPropertyValue=4000, maxJobPropertyValue=10000).value
		realignerTargetIntervalFile = File(os.path.join(topOutputDirJob.output, '%s_%s.forIndelRealigner.intervals'%\
													(bamFnamePrefix, overlapFilenameSignature)))
		realignerTargetIntervalJob = self.addGATKRealignerTargetCreatorJob(executable=self.RealignerTargetCreatorJava, \
															GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
						refFastaFList=passingData.refFastaFList, inputFile=bamF, inputArgumentOption="-I", indelVCFFile=indelVCFFile, \
						outputFile=realignerTargetIntervalFile, interval=overlapInterval, \
				parentJobLs=[topOutputDirJob]+alignmentData.jobLs+indelVCFFileJobLs, transferOutput=False, \
				job_max_memory=max(4000, indelRealignmentJobMaxMemory/2),walltime=indelRealignmentJobWalltime/2,\
				extraArguments=None, extraArgumentList=None, extraDependentInputLs=[baiF]+indelVCFFileLs)
		
		realignedBamFile = File(os.path.join(topOutputDirJob.output, '%s_%s.indelRealigned.bam'%\
											(bamFnamePrefix, overlapFilenameSignature)))
		#2013.04.09 GATK generates this file. it is not .bam.bai but just .bai. 
		realignedBaiFile = File('%s.bai'%(os.path.splitext(realignedBamFile.name)[0]))
		extraArgumentList=['-targetIntervals',realignerTargetIntervalJob.output,\
						'--read_filter NotPrimaryAlignment', \
						'--maxReadsForConsensuses 250', '--maxReadsForRealignment 90000', '--maxReadsInMemory 300000']
		if indelVCFFile:
			extraArgumentList.extend(["-known:vcf", indelVCFFile])	#"--consensusDeterminationModel KNOWNS_ONLY" is not added since vervet indels are not clear
		
		indelRealignmentJob = self.addGATKJob(executable=self.IndelRealignerJava, GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
					GATKAnalysisType='IndelRealigner',\
					inputFile=bamF, inputArgumentOption="-I", refFastaFList=passingData.refFastaFList, inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					interval=overlapInterval, outputFile=realignedBamFile, \
					parentJobLs=[realignerTargetIntervalJob]+alignmentData.jobLs+indelVCFFileJobLs, transferOutput=False, \
					job_max_memory=indelRealignmentJobMaxMemory,\
					frontArgumentList=None, extraArguments=None, \
					extraArgumentList=extraArgumentList, \
					extraOutputLs=[realignedBaiFile], \
					extraDependentInputLs=[realignerTargetIntervalJob.output, baiF]+indelVCFFileLs, no_of_cpus=None, \
					walltime=indelRealignmentJobWalltime)
		"""
		# 2013.04.07 Sun
		--maxReadsForConsensuses / -greedy ( int with default value 120 )
		max reads used for finding the alternate consensuses (necessary to improve performance in deep coverage). For expert users only! If you need to find the optimal solution regardless of running time, use a higher number.
		
		--maxReadsForRealignment / -maxReads ( int with default value 20000 )
		max reads allowed at an interval for realignment. For expert users only! If this value is exceeded at a given interval, realignment is not attempted and the reads are passed to the output file(s) as-is. If you need to allow more reads (e.g. with very deep coverage) regardless of memory, use a higher number.
		
		--maxReadsInMemory / -maxInMemory ( int with default value 150000 )
		max reads allowed to be kept in memory at a time by the SAMFileWriter. For expert users only! To minimize memory consumption you can lower this number (but then the tool may skip realignment on regions with too much coverage; and if the number is too low, it may generate errors during realignment). Just make sure to give Java enough memory! 4Gb should be enough with the default value.
		
		"""
		
		"""
		# 2013.04.05 Fri
		# '--read_filter NotPrimaryAlignment' is due to this error
		Error caching SAM record, ..., which is usually caused by malformed SAM/BAMfiles in which multiple identical copies of a read are present.
		"""
		
		sortBamF = File(os.path.join(topOutputDirJob.output, '%s_%s.indelRealigned.sorted.bam'%\
											(bamFnamePrefix, overlapFilenameSignature)))
		sortAlignmentJob = self.addSortAlignmentJob(workflow=workflow, inputBamFile=indelRealignmentJob.output, \
					outputBamFile=sortBamF,\
					SortSamFilesJava=self.SortSamFilesJava, SortSamJar=self.SortSamJar,\
					parentJobLs=[indelRealignmentJob], extraDependentInputLs=indelRealignmentJob.outputLs[1:], \
					extraArguments=None, job_max_memory =max(3000, indelRealignmentJobMaxMemory/3), \
					walltime=max(120, indelRealignmentJobWalltime/3), \
					transferOutput=False)
		"""
		# 2013.03.31 add the index job on bam file
		indexRealignedBamJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
										BuildBamIndexJar=self.BuildBamIndexJar, \
					inputBamF=sortAlignmentJob.output,\
					parentJobLs=[sortAlignmentJob], \
					transferOutput=transferOutput, job_max_memory=3000, \
					walltime=max(120, int(indelRealignmentJobWalltime/3)))
		"""		
		countCovariatesParentJobLs = [sortAlignmentJob, sortAlignmentJob.bamIndexJob]
		countCovariatesJobInput = sortAlignmentJob.output
		countCovariatesJobExtraDependentInputLs = [sortAlignmentJob.bamIndexJob.output]
		#else:
		#	countCovariatesParentJobLs = parentJobLs
		#	countCovariatesJobInput = bamF
			
		
		if span<1000000 and self.candidateCountCovariatesJob:	#2013.04.09 use the candidate job (interval big enough)
			countCovariatesJob = self.candidateCountCovariatesJob
		else:
			recalFile = File(os.path.join(topOutputDirJob.output, '%s_%s.recal_data.grp'%(bamFnamePrefix, overlapFilenameSignature)))
			countCovariatesJob = self.addGATKBaseRecalibratorJob(GenomeAnalysisTKJar=workflow.GenomeAnalysisTK2Jar, \
								inputFile=countCovariatesJobInput, \
								VCFFile=SNPVCFFile, interval=overlapInterval, outputFile=recalFile, \
								refFastaFList=passingData.refFastaFList, parentJobLs=[topOutputDirJob] + countCovariatesParentJobLs + SNPVCFJobLs, 
								extraDependentInputLs=[SNPVCFFile.tbi_F] + countCovariatesJobExtraDependentInputLs, \
								transferOutput=False, \
								extraArguments=None, job_max_memory=max(4000, indelRealignmentJobMaxMemory/2), \
								walltime=indelRealignmentJobWalltime/2)
		if span>self.intervalSize and self.candidateCountCovariatesJob is None:	#big chromosomes are first encountered so this should happen in 1st call()
			self.candidateCountCovariatesJob = countCovariatesJob
		
		"""
		countCovariatesJob = mapEachChromosomeData.countCovariatesJob
		"""
		
		recalBamFile = File(os.path.join(topOutputDirJob.output, '%s_%s.recal_data.bam'%(bamFnamePrefix, overlapFilenameSignature)))
		#2013.04.09 GATK generates this file. it is not .bam.bai but just .bai. 
		recalBaiFile = File('%s.bai'%(os.path.splitext(recalBamFile.name)[0]))
		selectAlignmentParentJob, selectAlignmentParentBamIndexJob = self.addGATKPrintRecalibratedReadsJob(GenomeAnalysisTKJar=workflow.GenomeAnalysisTK2Jar, \
							inputFile=bamF, \
							recalFile=countCovariatesJob.recalFile, interval=overlapInterval, outputFile=recalBamFile, \
							refFastaFList=passingData.refFastaFList, parentJobLs=[countCovariatesJob,], \
							extraDependentInputLs=[baiF], transferOutput=False, \
							extraArguments=None, job_max_memory=max(4000, indelRealignmentJobMaxMemory/3), \
							needBAMIndexJob=True, walltime=indelRealignmentJobWalltime/2)
		selectAlignmentParentJobLs = [selectAlignmentParentJob, selectAlignmentParentBamIndexJob]
		selectAlignmentParentJobOutput = selectAlignmentParentJob.output
		selectAlignmentJobExtraDependentInputLs = [selectAlignmentParentBamIndexJob.outputLs[1:]]
		"""
		else:
			sys.stderr.write("Warning from mapEachInterval(): no SNP VCF for BQSR for chromosome %s. skip BQSR.\n"%(chromosome))
			if countCovariatesParentJobLs:
				selectAlignmentParentJobLs = countCovariatesParentJobLs
				selectAlignmentParentJobOutput = countCovariatesJobInput
				selectAlignmentJobExtraDependentInputLs = countCovariatesJobExtraDependentInputLs
			else:
				selectAlignmentParentJobLs = parentJobLs
				selectAlignmentParentJobOutput = bamF
				selectAlignmentJobExtraDependentInputLs = [baiF]
		"""
		nonOverlapBamFile = File(os.path.join(topOutputDirJob.output, '%s_%s.bam'%(bamFnamePrefix, intervalFnameSignature)))
		selectAlignmentJob, indexSelectedAlignmentJob = self.addSelectAlignmentJob(executable=workflow.samtools, \
				inputFile=selectAlignmentParentJobOutput, \
				outputFile=nonOverlapBamFile, region=mpileupInterval, parentJobLs=selectAlignmentParentJobLs, \
				extraDependentInputLs=selectAlignmentJobExtraDependentInputLs, transferOutput=False, \
				extraArguments=None, job_max_memory=2000, needBAMIndexJob=False, \
				walltime=indelRealignmentJobWalltime/2)	#the next mergeBamJob doesn't need bai files.
		
		passingData.AlignmentJobAndOutputLs.append(PassingData(parentJobLs=[selectAlignmentJob, indexSelectedAlignmentJob], \
															file=selectAlignmentJob.output))
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
										mask_genotype_method_id=self.mask_genotype_method_id,\
										data_dir=self.data_dir, local_realigned=1)
			
			baseCoverage = 4	#baseline
			actualCoverage = getattr(individual_alignment.individual_sequence, 'coverage', baseCoverage)
			minMergeAlignmentWalltime = 240	#in minutes, 4 hours, when coverage is defaultCoverage
			maxMergeAlignmentWalltime = 2880	#in minutes, 2 days
			minMergeAlignmentMaxMemory = 5000	#in MB, when coverage is defaultCoverage
			maxMergeAlignmentMaxMemory = 12000	#in MB
			
			mergeAlignmentWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=actualCoverage, \
									baseInputVolume=baseCoverage, baseJobPropertyValue=minMergeAlignmentWalltime*2, \
									minJobPropertyValue=minMergeAlignmentWalltime, maxJobPropertyValue=maxMergeAlignmentWalltime).value
			mergeAlignmentMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=actualCoverage, \
									baseInputVolume=baseCoverage, baseJobPropertyValue=minMergeAlignmentMaxMemory, \
									minJobPropertyValue=minMergeAlignmentMaxMemory, maxJobPropertyValue=maxMergeAlignmentMaxMemory).value
			
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
									extraArguments=None, job_max_memory = 2500, transferOutput=False,\
									walltime=max(180, mergeAlignmentWalltime/20))
				
				NewAlignmentJobAndOutputLs.append(PassingData(parentJobLs=[addRGJob], file=addRGJob.output))
			#
			
			
			mergedBamFile = File(os.path.join(reduceOutputDirJob.output, '%s_recal.bam'%(bamFnamePrefix)))
			alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=NewAlignmentJobAndOutputLs, \
					outputBamFile=mergedBamFile, \
					samtools=workflow.samtools, java=workflow.java, \
					MergeSamFilesJava=workflow.MergeSamFilesJava, MergeSamFilesJar=workflow.MergeSamFilesJar, \
					BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexJar=workflow.BuildBamIndexJar, \
					mv=workflow.mv, parentJobLs=[reduceOutputDirJob], walltime=mergeAlignmentWalltime,\
					job_max_memory=mergeAlignmentMaxMemory, transferOutput=False)
			#2012.9.19 add/copy the alignment file to db-affliated storage
			#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
			logFile = File(os.path.join(reduceOutputDirJob.output, '%s_2db.log'%(bamFnamePrefix)))
			alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
								inputFile=alignmentMergeJob.output, otherInputFileList=[],\
								individual_alignment_id=new_individual_alignment.id, \
								mask_genotype_method_id=self.mask_genotype_method_id,\
								logFile=logFile, data_dir=data_dir, \
								parentJobLs=[alignmentMergeJob, bamIndexJob], \
								extraDependentInputLs=[bamIndexJob.output], \
								extraArguments=None, transferOutput=transferOutput, \
								job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel, commit=True,\
								walltime=max(180, mergeAlignmentWalltime/2))
			self.no_of_jobs += 1
			returnData.jobDataLs.append(PassingData(jobLs=[alignment2DBJob], file=alignment2DBJob.logFile, \
											fileLs=[alignment2DBJob.logFile]))
		return returnData

	
	
	def addGATKBaseRecalibratorJob(self, workflow=None, executable=None, GenomeAnalysisTKJar=None, inputFile=None, \
								VCFFile=None, interval=None, outputFile=None, \
					refFastaFList=[], parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1, **keywords):
		"""
		2013.2.14 upgraded to GATK2's "-T BaseRecalibrator"
		http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr
			java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -I my_reads.bam -R resources/Homo_sapiens_assembly18.fasta
				-knownSites bundle/hg18/dbsnp_132.hg18.vcf \
				-knownSites another/optional/setOfSitesToMask.vcf \
				-o recal_data.grp
		2011-12-4
			inputFile is a bam file.
			outputFile is recalFile (csv file).
		
			java -jar ~/script/vervet/bin/GenomeAnalysisTK-1.0.4705/GenomeAn alysisTK.jar -l INFO
				-R ../../../../NCBI/hs_genome.fasta -I 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.bam
				-T CountCovariates   -cov ReadGroupCovariate    -cov QualityScoreCovariate  
				-cov CycleCovariate    -cov DinucCovariate
				-recalFile  454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal_data.csv
				-B:mask,VCF 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK.vcf
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = workflow.BaseRecalibratorJava
		if GenomeAnalysisTKJar is None:
			GenomeAnalysisTKJar = workflow.GenomeAnalysisTK2Jar
		#GATK job
		#MaxPermSize= min(35000, max(1024, job_max_memory*7/9))
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory)
		#javaMemRequirement = "-Xms%sm -Xmx%sm -XX:PermSize=%sm -XX:MaxPermSize=%sm"%(job_max_memory*50/100, job_max_memory, \
		#																			MaxPermSize*50/100, MaxPermSize)
		refFastaFile = refFastaFList[0]
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', GenomeAnalysisTKJar, "-T BaseRecalibrator",\
						"-I", inputFile, "-R", refFastaFile,\
						"-L", interval, self.defaultGATKArguments, \
						"-knownSites:vcf", VCFFile,\
						"--out", outputFile]
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputFile, VCFFile, GenomeAnalysisTKJar] + refFastaFList)
		
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[outputFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, \
						job_max_memory=memRequirementData.memRequirement, **keywords)
		
		job.recalFile = outputFile
		return job
	
	def addGATKPrintRecalibratedReadsJob(self, workflow=None, executable=None, GenomeAnalysisTKJar=None, inputFile=None, \
								recalFile=None, interval=None, outputFile=None, extraOutputLs=None, \
					refFastaFList=[], parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1, needBAMIndexJob=False,\
					walltime=None,**keywords):
		"""
		2013.04.09 added extraOutputLs
		2013.2.14 upgraded to GATK2's "-T PrintReads"
			java -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fasta -I input.bam -BQSR recalibration_report.grp 
					-o output.bam
		2012.7.27
			java -jar ~/script/vervet/bin/GenomeAnalysisTK-1.0.4705/GenomeAnalysisTK.jar -l INFO 
				-R ../../../../NCBI/hs_genome.fasta -I 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.bam 
				-T TableRecalibration  --out 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal.bam 
				-recalFile 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal_data.csv
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = workflow.PrintRecalibratedReadsJava
		if GenomeAnalysisTKJar is None:
			GenomeAnalysisTKJar = workflow.GenomeAnalysisTK2Jar
		if extraOutputLs is None:
			extraOutputLs = []
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory)
		#MaxPermSize= min(35000, max(1024, job_max_memory*7/9))
		javaMemRequirement = memRequirementData.memRequirementInStr
		job_max_memory = memRequirementData.memRequirement
		
		#"-Xms%sm -Xmx%sm -XX:PermSize=%sm -XX:MaxPermSize=%sm"%(job_max_memory*50/100, job_max_memory, \
		#																			MaxPermSize*50/100, MaxPermSize)
		refFastaFile = refFastaFList[0]
		extraArgumentList = [javaMemRequirement, '-jar', GenomeAnalysisTKJar, "-T PrintReads",\
						"-I", inputFile, "-R", refFastaFile,\
						"-L", interval, self.defaultGATKArguments,\
						"-BQSR", recalFile, "--out", outputFile]
		if extraArguments:
			extraArgumentList.append(extraArguments)
		extraDependentInputLs.extend([inputFile, recalFile, GenomeAnalysisTKJar] + refFastaFList)
		extraOutputLs.append(outputFile)
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		if needBAMIndexJob:
			# add the index job on the bam file
			bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
						BuildBamIndexJar=self.BuildBamIndexJar, \
						inputBamF=job.output, parentJobLs=[job], \
						transferOutput=transferOutput, job_max_memory=job_max_memory, walltime=walltime/2)
		else:
			bamIndexJob = None
		return job, bamIndexJob
	
	def setup(self, inputVCFData=None, chr2IntervalDataLs=None, **keywords):
		"""
		2013.04.07
			register the indel vcfs
		"""
		pdata = parentClass.setup(self, inputVCFData=inputVCFData, chr2IntervalDataLs=chr2IntervalDataLs, **keywords)
		
		##2013.04.09 used in case a corresponding VCFFile for one interval is not available for BQSR
		biggestChromosome = pdata.chrSizeIDList[0][1]
		self.randomSNPVCFJobDataForBQSR = pdata.chr2VCFJobData[biggestChromosome]
		
		self.indelVCFInputData = self.registerAllInputFiles(inputDir=self.indelVCFFolder, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName="%sIndelVCF"%self.pegasusFolderName)
		#2012.8.26 so that each recalibration will pick up the right vcf
		chr2IndelVCFJobData = {}
		if self.indelVCFInputData:
			for jobData in self.indelVCFInputData.jobDataLs:
				inputF = jobData.file
				chromosome = self.getChrFromFname(os.path.basename(inputF.name))
				chr2IndelVCFJobData[chromosome] = jobData
		pdata.chr2IndelVCFJobData = chr2IndelVCFJobData
		self.chr2IndelVCFJobData = chr2IndelVCFJobData
		
		return pdata
	
	def isThisAlignmentComplete(self, individual_alignment=None, data_dir=None):
		"""
		2013.04.09 this is more complicated as it tests the local_realigned version of individual_alignment is complete or not.
			not individual_alignment itself
		"""
		new_individual_alignment = self.db.copyParentIndividualAlignment(parent_individual_alignment_id=individual_alignment.id,\
										mask_genotype_method_id=self.mask_genotype_method_id,\
										data_dir=self.data_dir, local_realigned=1)
		return self.db.isThisAlignmentComplete(individual_alignment=new_individual_alignment, data_dir=data_dir)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		#samtools is only used for select alignment, which is very fast, increase the clustering 
		self.setOrChangeExecutableClusterSize(executable=workflow.samtools, clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='BaseRecalibratorJava', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='PrintRecalibratedReadsJava', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='RealignerTargetCreatorJava', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='IndelRealignerJava', clusterSizeMultipler=1)


if __name__ == '__main__':
	main_class = AlignmentReadBaseQualityRecalibrationWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
