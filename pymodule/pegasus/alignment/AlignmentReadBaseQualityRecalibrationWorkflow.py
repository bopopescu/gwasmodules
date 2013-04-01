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
		-e /u/home/eeskin/polyacti
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 17

	# 2012.9.18
	%s  --inputDir ~/NetworkData/vervet/db/genotype_file/method_41 --ind_seq_id_ls 633,634,635,636,637,638 
		--ref_ind_seq_id 524
		-o dags/BaseQualityRecalibration/BaseQualityRecalibration_ISQ633_638_vsMethod41.xml -l hcondor
		-j hcondor -z localhost -u yh --intervalSize 10000000 --intervalOverlapSize 30000
		-e /u/home/eeskin/polyacti
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--clusters_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 41
	
	# 2013.3.19 use sequence coverage to filter alignments
	%s  --inputDir ~/NetworkData/vervet/db/genotype_file/method_41
		--sequence_min_coverage 0 --sequence_max_coverage 2  --ind_seq_id_ls 632-3230
		--ref_ind_seq_id 3280 -o dags/BaseQualityRecalibration/BaseQualityRecalibration_ISQ632_3230_coverage0_2_vsMethod41.xml
		-l hcondor -j hcondor -z localhost -u yh --intervalSize 10000000 --intervalOverlapSize 30000
		-e /u/home/eeskin/polyacti --contigMaxRankBySize 250
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--clusters_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 41
		# --ref_genome_version 2 #(optional, as by default, it gets the outdated_index=0 reference chromosomes from GenomeDB)
		# --ref_genome_outdated_index 1 #to get old reference. incompatible here as alignment is based on 3280, new ref.
		# --needFastaDictJob --needFastaIndexJob
	
Description:
	#2012.9.21  recalibration needs ~500K reads to get accurate estimate. So adjust the --intervalSize according to coverage.  
	Input: a VCF folder, list of alignment IDs (or use site_id, ref_ind_seq_id to filter )
	
	for each whole-genome alignment
	select a whole-genome alignment into several interval alignment (5Mb of each chr/contig)
			using samtools view -h ... (AlignmentToCallPipeline.addSelectAndSplitBamJobs())
		each sub-alignment is paired with the VCF from the same contig/chromosome.
		run CountCovariates job (skeleton in ShortRead2AlignmentPipeline.py)
		TableRecalibration job
	Merge all sub-alignment into one
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
from vervet.src import VervetDB, AbstractVervetAlignmentAndVCFWorkflow

parentClass = AbstractVervetAlignmentAndVCFWorkflow
class AlignmentReadBaseQualityRecalibrationWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(parentClass.partitionWorkflowOptionDict.copy())
	option_default_dict.update({
				('mask_genotype_method_id', 0, int):[None, '', 1, 'which genotype method is used to mask out polymorphic sites for recalibration'],\
							})
	option_default_dict[('intervalSize', 1, int)][0] = 10000000
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
				VCFFile=None, passingData=None, reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
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
											fileList=[countCovariatesJob.recalFile]))
		"""
		
		return returnData
	
	def mapEachInterval(self, workflow=None, alignmentData=None, intervalData=None,\
							VCFFile=None, passingData=None, reduceBeforeEachAlignmentData=None,\
							mapEachChromosomeData=None, transferOutput=False, \
							**keywords):
		"""
		2013.03.31 use VCFFile to decide whether to add BQSR jobs, called in ShortRead2AlignmentWorkflow.py
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
		
		
		if intervalData.file:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.file
		else:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.interval
		intervalFnameSignature = intervalData.intervalFnameSignature
		overlapInterval = intervalData.overlapInterval
		overlapFilenameSignature = intervalData.overlapIntervalFnameSignature
		
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
		realignerTargetIntervalFile = File(os.path.join(topOutputDirJob.output, '%s_%s.forIndelRealigner.intervals'%\
													(bamFnamePrefix, overlapFilenameSignature)))
		realignerTargetIntervalJob = self.addGATKJob(executable=self.RealignerTargetCreatorJava, GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
					GATKAnalysisType='RealignerTargetCreator',\
					inputFile=bamF, inputArgumentOption="-I", refFastaFList=passingData.refFastaFList, inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					interval=overlapInterval, outputFile=realignerTargetIntervalFile, \
					parentJobLs=[topOutputDirJob]+parentJobLs, transferOutput=False, job_max_memory=4000,\
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=[baiF], no_of_cpus=None, walltime=300)
		
		realignedBamFile = File(os.path.join(topOutputDirJob.output, '%s_%s.indelRealigned.bam'%\
											(bamFnamePrefix, overlapFilenameSignature)))
		indelRealignmentJobWalltime=300
		indelRealignmentJob = self.addGATKJob(executable=self.IndelRealignerJava, GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
					GATKAnalysisType='IndelRealigner',\
					inputFile=bamF, inputArgumentOption="-I", refFastaFList=passingData.refFastaFList, inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					interval=overlapInterval, outputFile=realignedBamFile, \
					parentJobLs=[realignerTargetIntervalJob]+parentJobLs, transferOutput=False, job_max_memory=7000,\
					frontArgumentList=None, extraArguments=None, \
					extraArgumentList=['-targetIntervals',realignerTargetIntervalJob.output], \
					extraOutputLs=None, \
					extraDependentInputLs=[realignerTargetIntervalJob.output, baiF], no_of_cpus=None, \
					walltime=indelRealignmentJobWalltime)
		
		# 2013.03.31 add the index job on bam file
		indexRealignedBamJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
										BuildBamIndexJar=self.BuildBamIndexJar, \
					inputBamF=indelRealignmentJob.output,\
					parentJobLs=[indelRealignmentJob], \
					transferOutput=transferOutput, job_max_memory=3000, \
					walltime=max(180, int(indelRealignmentJobWalltime/3)))
		if VCFFile:	#2013.03.31
			recalFile = File(os.path.join(topOutputDirJob.output, '%s_%s.recal_data.grp'%(bamFnamePrefix, overlapFilenameSignature)))
			countCovariatesJob = self.addGATKBaseRecalibratorJob(GenomeAnalysisTKJar=workflow.GenomeAnalysisTK2Jar, \
									inputFile=indelRealignmentJob.output, \
									VCFFile=VCFFile, interval=overlapInterval, outputFile=recalFile, \
									refFastaFList=passingData.refFastaFList, parentJobLs=[topOutputDirJob, \
																indelRealignmentJob, indexRealignedBamJob], 
									extraDependentInputLs=[indexRealignedBamJob.output, VCFFile.tbi_F], \
									transferOutput=True, \
									extraArguments=None, job_max_memory=4000, walltime=300)
			"""
			countCovariatesJob = mapEachChromosomeData.countCovariatesJob
			"""
			
			recalBAMFile = File(os.path.join(topOutputDirJob.output, '%s_%s.recal_data.bam'%(bamFnamePrefix, overlapFilenameSignature)))
			selectAlignmentParentJob, selectAlignmentParentBamIndexJob = self.addGATKPrintRecalibratedReadsJob(GenomeAnalysisTKJar=workflow.GenomeAnalysisTK2Jar, \
								inputFile=bamF, \
								recalFile=countCovariatesJob.recalFile, interval=overlapInterval, outputFile=recalBAMFile, \
								refFastaFList=passingData.refFastaFList, parentJobLs=[countCovariatesJob, indexRealignedBamJob], \
								extraDependentInputLs=[indexRealignedBamJob.output, VCFFile.tbi_F], transferOutput=False, \
								extraArguments=None, job_max_memory=3000, needBAMIndexJob=True, walltime=300)
		else:
			selectAlignmentParentJob = indelRealignmentJob
			selectAlignmentParentBamIndexJob = indexRealignedBamJob
		nonOverlapBamFile = File(os.path.join(topOutputDirJob.output, '%s_%s.bam'%(bamFnamePrefix, intervalFnameSignature)))
		selectAlignmentJob, bamIndexJob3 = self.addSelectAlignmentJob(executable=workflow.samtools, \
				inputFile=selectAlignmentParentJob.output, \
				outputFile=nonOverlapBamFile, region=mpileupInterval, parentJobLs=[selectAlignmentParentJob, \
														selectAlignmentParentBamIndexJob], \
				extraDependentInputLs=[selectAlignmentParentBamIndexJob.output], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, needBAMIndexJob=False)	#the next mergeBamJob doesn't need bai files.
		
		passingData.AlignmentJobAndOutputLs.append(PassingData(parentJobLs=[selectAlignmentJob,bamIndexJob3], \
															file=selectAlignmentJob.output))
		#add the sub-alignment to the alignment merge job
		self.no_of_jobs += 5
		return returnData
	
	def reduceAfterEachAlignment(self, workflow=None, passingData=None, transferOutput=False, data_dir=None, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		if workflow is None:
			workflow = self
		getattr(passingData, 'mapEach', [])
		AlignmentJobAndOutputLs = getattr(passingData, 'AlignmentJobAndOutputLs', [])
		bamFnamePrefix = passingData.bamFnamePrefix
		topOutputDirJob = passingData.topOutputDirJob
		individual_alignment = passingData.individual_alignment
		if len(AlignmentJobAndOutputLs)>0:	#2012.3.29	merge alignment output only when there is something to merge!
			mergedBamFile = File(os.path.join(topOutputDirJob.output, '%s_recal.bam'%(bamFnamePrefix)))
			alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=AlignmentJobAndOutputLs, \
					outputBamFile=mergedBamFile, \
					samtools=workflow.samtools, java=workflow.java, \
					MergeSamFilesJava=workflow.MergeSamFilesJava, MergeSamFilesJar=workflow.MergeSamFilesJar, \
					BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexJar=workflow.BuildBamIndexJar, \
					mv=workflow.mv, parentJobLs=[topOutputDirJob], \
					transferOutput=False)
			#2012.9.19 add/copy the alignment file to db-affliated storage
			#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
			logFile = File(os.path.join(topOutputDirJob.output, '%s_2db.log'%(bamFnamePrefix)))
			alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
								inputFile=alignmentMergeJob.output, otherInputFileList=[],\
								parent_individual_alignment_id=individual_alignment.id, \
								mask_genotype_method_id=self.mask_genotype_method_id,\
								logFile=logFile, data_dir=data_dir, \
								parentJobLs=[alignmentMergeJob, bamIndexJob], \
								extraDependentInputLs=[bamIndexJob.output], \
								extraArguments=None, transferOutput=transferOutput, \
								job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel, commit=True)
			self.no_of_jobs += 1
			returnData.jobDataLs.append(PassingData(jobLs=[alignment2DBJob], file=alignment2DBJob.logFile, \
											fileList=[alignment2DBJob.logFile]))
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
						"-knownSites", VCFFile,\
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
								recalFile=None, interval=None, outputFile=None, \
					refFastaFList=[], parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1, needBAMIndexJob=False, **keywords):
		"""
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
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputFile, recalFile, GenomeAnalysisTKJar] + refFastaFList)
		
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
		return job, bamIndexJob
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='BaseRecalibratorJava', clusterSizeMultipler=0.5)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='PrintRecalibratedReadsJava', clusterSizeMultipler=0.5)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='RealignerTargetCreatorJava', clusterSizeMultipler=0.5)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='IndelRealignerJava', clusterSizeMultipler=0.1)


if __name__ == '__main__':
	main_class = AlignmentReadBaseQualityRecalibrationWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()