#!/usr/bin/env python
"""
Examples:
	# 2011-8-30 workflow on condor, always commit (-c)
	%s -i 165-167 -o ShortRead2Alignment_isq_id_165_167_vs_9.xml -u yh -a 9 -l condorpool
		-n1 -z dl324b-1.cmb.usc.edu -c -H
	
	# 2011-8-30 a workflow with 454 long-read and short-read PE. need a ref index job (-n1). 
	%s -i 165-167 -o ShortRead2Alignment_isq_id_165_167_vs_9.xml -u yh -a 9
		-e /u/home/eeskin/polyacti -l hoffman2 -t /u/home/eeskin/polyacti/NetworkData/vervet/db -n1
		-z dl324b-1.cmb.usc.edu -c
		--tmpDir /work/ -H
	
	# 2011-8-30 output a workflow to run alignments on hoffman2's condor pool (-D changes local_data_dir. -t changes data_dir.)
	# 2012.3.20 use /work/ or /u/scratch/p/polyacti/tmp as TMP_DIR for MarkDuplicates.jar (/tmp is too small for 30X genome)
	# 2012.5.4 cluster 10 alignment jobs (before merging) as a unit (--cluster_size_for_aln_jobs 10), skip done alignment (-K)
	# 2012.9.21 add "-H" because AddAlignmentFile2DB need db conneciton
	# 2012.9.21 add "--alignmentPerLibrary" to also get alignment for each library within one individual_sequence
	%s  -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
		-l hcondor -j hcondor 
		-z localhost -u yh -c
		-i 631-700 -o workflow/ShortRead2Alignment_Isq_631-700_vs_524_hcondor.xml  -a 524 
		--tmpDir /work/ -e /u/home/eeskin/polyacti  --cluster_size_for_aln_jobs 10 -K  -H --alignmentPerLibrary
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files. need ref index job (--needRefIndexJob).
	# If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
	%s -i 176,178-183,207-211
		-o ShortRead2Alignment_8VWP_vs_9_condor_no_refIndex.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z dl324b-1.cmb.usc.edu -p secret  -c -H
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files.
	# If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
	%s -i 176,178-183,207-211
		-o ShortRead2Alignment_8VWP_vs_9_condor_no_refIndex.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z dl324b-1.cmb.usc.edu -p secret  -c -H
		
	# 2011-8-30 a workflow to run on uschpc, with ref index job. Note the site_handler and input_site_handler.
	# to enable replica-transfer.
	%s -i 391-397,456,473,493
		-o ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j local -l uschpc -n1 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10 -p secret  -c -H

	# 2011-8-30 a workflow to run on uschpc, Need ref index job (--needRefIndexJob), and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler and "-t ..." to enable symlink
	%s -i 391-397,456,473,493
		-o ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j uschpc -l uschpc --needRefIndexJob -p secret -c --no_of_aln_threads 4 -H
		-e /home/cmb-03/mn/yuhuang -z 10.8.0.10 
		-t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ -J /home/cmb-03/mn/yuhuang/bin/jdk/bin/java
	
	# 2011-11-16 a workflow to run on uschpc, Need ref index job (--needRefIndexJob), and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler. this will stage in all input and output (--notStageOutFinalOutput).
	%s -i 391-397,456,473,493
		-o workflow/ShortRead2Alignment/ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_local2usc.xml -u yh -a 9
		-j local -l uschpc --needRefIndexJob -p secret -c --no_of_aln_threads 4
		-e /home/cmb-03/mn/yuhuang -z 10.8.0.10 
		-J /home/cmb-03/mn/yuhuang/bin/jdk/bin/java  -H
	
	
	#2011-9-13 no ref index job, staging input files from localhost to uschpc, stage output files back to localhost
	# modify the refFastaFile's path in xml manually
	%s -i 1-3 -o ShortRead2Alignment_1_3_vs_524_local2uschpc.xml -u yh -a 524
		-j local -l uschpc --needRefIndexJob -p secret -c --no_of_aln_threads 4 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10
		-t /Network/Data/vervet/db/  -H
	
	# 2011-8-31 output the same workflow above but for condorpool
	%s -i 391-397,456,473,493, -o workflow/ShortRead2Alignment/ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_condorpool.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z 10.8.0.10  -p secret  -c --alignmentPerLibrary  -H
	
	# 2012-4-5 new alignment method, stampy (--alignment_method_name)
	%s -i 167,176,178,182,183,207-211,391-397,456,473,493
		-o workflow/ShortRead2Alignment/ShortRead2Alignment_10VWP_4DeepVRC_6LowCovVRC_392_397_vs_508_condorpool.xml
		-u yh -a 508 -j condorpool -l condorpool -n1 -z 10.8.0.10  -p secret  -c --alignment_method_name stampy  -H
	
Description:
	2012.11.14 a read-alignment workflow, extracted from vervet/src/ShortRead2AlignmentPipeline.py.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
				sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, yh_pegasus
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow


class ShortRead2AlignmentWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	alignment_option_dict = {
						('noCheckEmptyReadFile', 0, int):[0, 'R', 0, "toggle to not check whether each read file is empty (if empty, exclude it). \
							If IndividualSequenceFile.read_count is null, it'll try to count them on the fly and take a lot of time.\
							however, only toggle it if you know every input individual_sequence_file is not empty. empty read file fails alignment jobs."],\
						('additionalArguments', 0, ): ["-q 20", '', 1, 'a string of additional arguments passed to aln, not bwasw, add double quote if space'],\
						("bwa_path", 1, ): ["%s/bin/bwa", '', 1, 'bwa binary'],\
						("stampy_path", 1, ): ["%s/bin/stampy.py", '', 1, 'path to stampy.py'],\
						("alignment_method_name", 1, ): ["bwa-short-read", '', 1, 'alignment_method.short_name from db.\
								used only when unable to guess based on individual_sequence.sequencer and individual_sequence.sequence_type'],\
						("needRefIndexJob", 0, int): [0, '', 1, 'need to add a reference index job by bwa?'],\
						('no_of_aln_threads', 1, int): [1, '', 1, 'number of threads during alignment'],\
						('cluster_size_for_aln_jobs', 1, float): [0.01, '', 1, 'cluster size relative to self.clusters_size, \n\
	for bwa/PEAlignmentByBWA/LongSEAlignmentByBWA/addOrReplaceReadGroupsJava/SortSamFilesJava/samtools jobs'],\
						('notStageOutFinalOutput', 0, int):[0, '', 0, 'toggle to not stage out final output (bam + bam.bai)'],\
						("tmpDir", 1, ): ["/tmp/", '', 1, 'for MarkDuplicates.jar, default is /tmp/ but sometimes it is too small'],\
						
							}
	option_default_dict.update(alignment_option_dict.copy())
	option_default_dict.update({
						('refSequenceFname', 1, ): ["", '', 1, 'path to the reference file', ],\
						})

	def __init__(self,  **keywords):
		"""
		2012.3.29
			default to stage out final output.
			Argument stageOutFinalOutput morphs into notStageOutFinalOutput.
		2011-7-11
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		
		self.bwa_path =  self.insertHomePath(self.bwa_path, self.home_path)
		self.stampy_path =  self.insertHomePath(self.stampy_path, self.home_path)
		
		if self.notStageOutFinalOutput:
			self.stageOutFinalOutput = False
		else:
			self.stageOutFinalOutput = True
		
		if self.noCheckEmptyReadFile:
			self.ignoreEmptyReadFile = False
		else:
			self.ignoreEmptyReadFile = True
	
	def registerCustomJars(self, workflow=None):
		"""
		2012.1.9
		"""
		if workflow is None:
			workflow = self
		AbstractNGSWorkflow.registerCustomJars(self, workflow=workflow)
		
		site_handler = self.site_handler
		abs_path = os.path.join(self.picard_path, 'SortSam.jar')
		SortSamFilesJar = File(abs_path)
		SortSamFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(SortSamFilesJar)
		workflow.SortSamFilesJar = SortSamFilesJar
		
		abs_path = os.path.join(self.picard_path, 'SamFormatConverter.jar')
		SamFormatConverterJar = File(abs_path)
		SamFormatConverterJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(SamFormatConverterJar)
		workflow.SamFormatConverterJar = SamFormatConverterJar
		
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.5.3
			add clusters.size profile for alignment specific jobs (self.cluster_size_for_aln_jobs)
		2012.1.3
		"""
		if workflow is None:
			workflow = self
		AbstractNGSWorkflow.registerCustomExecutables(self, workflow=workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)		
		stampy = Executable(namespace=namespace, name="stampy", version=version, os=operatingSystem, \
						arch=architecture, installed=True)
		stampy.addPFN(PFN("file://" + self.stampy_path, site_handler))
		executableClusterSizeMultiplierList.append((stampy, self.cluster_size_for_aln_jobs))
		
		#workflow.BuildBamIndexFilesJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		executableClusterSizeMultiplierList.append((workflow.addOrReplaceReadGroupsJava, self.cluster_size_for_aln_jobs))
		executableClusterSizeMultiplierList.append((workflow.samtools, self.cluster_size_for_aln_jobs))
		
		
		bwa = Executable(namespace=namespace, name="bwa", version=version, os=operatingSystem, arch=architecture, installed=True)
		bwa.addPFN(PFN("file://" + self.bwa_path, site_handler))
		executableClusterSizeMultiplierList.append((bwa, self.cluster_size_for_aln_jobs))
		
		PEAlignmentByBWA = Executable(namespace=namespace, name="PEAlignmentByBWA", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		PEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/PEAlignmentByBWA.sh"), site_handler))
		executableClusterSizeMultiplierList.append((PEAlignmentByBWA, self.cluster_size_for_aln_jobs))
		
		ShortSEAlignmentByBWA = Executable(namespace=namespace, name="ShortSEAlignmentByBWA", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		ShortSEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/ShortSEAlignmentByBWA.sh"), site_handler))
		executableClusterSizeMultiplierList.append((ShortSEAlignmentByBWA, self.cluster_size_for_aln_jobs))
		
		LongSEAlignmentByBWA = Executable(namespace=namespace, name="LongSEAlignmentByBWA", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		LongSEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/LongSEAlignmentByBWA.sh"), site_handler))
		executableClusterSizeMultiplierList.append((LongSEAlignmentByBWA, self.cluster_size_for_aln_jobs))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)

	def addBWAReferenceIndexJob(self, workflow=None, refFastaFList=None, refSequenceBaseCount=3000000000, bwa=None,\
					bwaIndexFileSuffixLs = None,\
					transferOutput=True, job_max_memory=4000, max_walltime=200):
		"""
		2012.10.10
			renamed from addRefIndexJob to addBWAReferenceIndexJob()
			bwaIndexFileSuffixLs is controlled by AbstractNGSWorkflow.bwaIndexFileSuffixLs.
			and the suffices of index output has changed ('rbwt', 'rpac', 'rsa', are gone):
				524_superContigsMinSize2000.fasta.bwt
				524_superContigsMinSize2000.fasta.pac
				524_superContigsMinSize2000.fasta.ann
				524_superContigsMinSize2000.fasta.amb
				524_superContigsMinSize2000.fasta.sa
				524_superContigsMinSize2000.fasta.nhr
				524_superContigsMinSize2000.fasta.nin
				524_superContigsMinSize2000.fasta.nsq

		2011-8-28
		"""
		if workflow is None:
			workflow = self
		if bwaIndexFileSuffixLs is None:
			bwaIndexFileSuffixLs = self.bwaIndexFileSuffixLs	#self.bwaIndexFileSuffixLs is defined in AbstractNGSWorkflow
		
		if refSequenceBaseCount is None:	#default	#maybe add a base count job here
			index_algorithm = 'bwtsw'
		elif refSequenceBaseCount<500000000:	#500 million
			index_algorithm = "is"
		else:
			index_algorithm = "bwtsw"
		
		refFastaFile = refFastaFList[0]
		extraArgumentList = ["index", "-a", index_algorithm, refFastaFile]
		extraOutputLs = []
		for suffix in bwaIndexFileSuffixLs:
			file = File("%s.%s"%(refFastaFile.name, suffix))
			extraOutputLs.append(file)
		extraDependentInputLs = refFastaFList
		bwa_index_job = self.addGenericJob(executable=bwa, inputFile=None, outputFile=None, \
						parentJobLs=None, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArguments=None, extraArgumentList=extraArgumentList, \
						key2ObjectForJob=None,\
						job_max_memory=job_max_memory, max_walltime=max_walltime)
		return bwa_index_job
	
	def addStampyGenomeIndexHashJob(self, workflow=None, executable=None, refFastaFList=None, \
						parentJobLs=None, job_max_memory=100, job_max_walltime = 60, \
						extraDependentInputLs=None, \
						transferOutput=True, **keywords):
		"""
		2012.10.10 use addGenericJob()
		2012.2.23
		"""
		"""
		stampyGenomeIndexJob= Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		refFastaFile = refFastaFList[0]
		
		genomeIndexFile = File('%s.stidx'%(refFastaFile.name))
		stampyGenomeIndexJob.addArguments("-G ", refFastaFile, refFastaFile)
		stampyGenomeIndexJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		#genomeIndexFile is output of this stampyGenomeIndexJob but ignore it as output otherwise it'll get auto-cleaned by pegasus
		#stampyGenomeIndexJob.uses(genomeIndexFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(stampyGenomeIndexJob)
		for input in extraDependentInputLs:
			stampyGenomeIndexJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=stampyGenomeIndexJob)
		yh_pegasus.setJobProperRequirement(stampyGenomeIndexJob, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		
		stampyGenomeHashJob = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		stampyGenomeHashJob.addArguments("-g", refFastaFile, "-H", refFastaFile)
		
		stampyGenomeHashJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		genomeHashFile = File('%s.sthash'%(refFastaFile.name))
		#genomeIndexFile is input to this stampyGenomeHashJob but mark it as output otherwise it'll get auto-cleaned by pegasus
		stampyGenomeHashJob.uses(genomeIndexFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		stampyGenomeHashJob.uses(genomeHashFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		stampyGenomeHashJob.outputLs = [genomeIndexFile, genomeHashFile]
		
		workflow.addJob(stampyGenomeHashJob)
		workflow.depends(parent=stampyGenomeIndexJob, child=stampyGenomeHashJob)
		yh_pegasus.setJobProperRequirement(stampyGenomeHashJob, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		"""
		
		refFastaFile = refFastaFList[0]
		genomeIndexFile = File('%s.stidx'%(refFastaFile.name))
		extraOutputLs = [genomeIndexFile]
		extraArgumentList = [refFastaFile]
		stampyGenomeIndexJob = self.addGenericJob(executable=executable, inputFile=refFastaFile, \
						inputArgumentOption="-G", outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArguments=None, extraArgumentList=extraArgumentList, \
						key2ObjectForJob=None,\
						job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		
		genomeHashFile = File('%s.sthash'%(refFastaFile.name))
		extraOutputLs = [genomeHashFile, genomeIndexFile]
		extraArgumentList = ["-H", refFastaFile]
		stampyGenomeHashJob = self.addGenericJob(executable=executable, inputFile=refFastaFile, \
						inputArgumentOption="-g", outputFile=None, \
						parentJobLs=[stampyGenomeIndexJob], extraDependentInputLs=[genomeIndexFile], \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArguments=None, extraArgumentList=extraArgumentList, \
						key2ObjectForJob=None,\
						job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		return stampyGenomeHashJob
	
	def registerISQFileObjLsToWorkflow(self, fileObjectLs=None, workflow=None, data_dir=None):
		'''
		2012-2.24
			similar to registerFileToWorkflow but for a different input
		'''
		
		newFilePair = []
		for fileObject in fileObjectLs:
			relativePath = fileObject.db_entry.path
			fastqF = File(relativePath)
			fastqF.addPFN(PFN("file://" + fileObject.path, self.input_site_handler))
			workflow.addFile(fastqF)
			fileObject.fastqF = fastqF
			newFilePair.append(fileObject)
		return newFilePair


	def addRefIndexJobAndItsOutputAsParent(self, workflow=None, refIndexJob=None, childJob=None):
		"""
		2012.10.18
			updated by using addJobDependency() and self.addJobUse()
		2012.2.24
		"""
		self.addJobDependency(workflow=workflow, parentJob=refIndexJob, childJob=childJob)
		for output in refIndexJob.outputLs:
			self.addJobUse(childJob, file=output, transfer=True, register=True, link=Link.INPUT)
			#childJob.uses(output, transfer=True, register=True, link=Link.INPUT)
	
	def addStampyAlignmentJob(self, workflow=None, fileObjectLs=None, \
					refFastaFList=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, parentJobLs=None,\
					alignment_method=None, outputDir=None, \
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					java=None, SortSamFilesJava=None, SortSamFilesJar=None,\
					addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
					no_of_aln_threads=3, stampy=None, transferOutput=False, **keywords):
		"""
		2014.4.4
			update for stampy 1.0.17
				which added option --gatkcigarworkaround removing adjacent I/D events from CIGAR strings, which trips up GATK
		2012.3.5
			added "--overwrite" to stampy.py to overwrite partial output alignment file
			added two stampy options
				--baq                              (SAM format) Compute base-alignment quality (BAQ; BQ tag)
				--alignquals                       (SAM format) Compute posterior alignment probabilities (YQ tag)
		2012.2.26
			handle read files with quality_score_format="Illumina" (post 1.3 solexa).
			
		2012.2.24
			alignment job for stampy
			no_of_aln_threads is only for bwa.
		2011-9-13
			add argument java & SortSamFilesJar
		2011-9-9
			two steps:
				1. aln doesn't use pipe, outputs to sai files.
				2. sampe/samse, convert, sort => connected through pipe
		"""
		aln_job_max_memory = 6000	#in MB, bwa needs 3G. stampy needs 3G.
		#bwa: memory 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		
		bwasw_job_max_memory = 3800	#in MB, same ref, bwasw needs more memory
		samse_job_max_memory = 4500	#in MB
		sampe_job_max_memory = 6000	#in MB
		addRGJob_max_memory = 2500	#in MB
		aln_job_max_walltime= 4800	#80 hours, in minutes
		aln_job_max_walltime = 1320	#22 hours, because all reads are stored in chunks of 5-million-read files
		
		javaMemRequirement = "-Xms128m -Xmx%sm"%addRGJob_max_memory
		refFastaFile = refFastaFList[0]
		firstFileObject = fileObjectLs[0]
		fastqF = firstFileObject.fastqF
		relativePath = fastqF.name
		fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
		outputSamFile = File('%s.sam'%(os.path.join(outputDir, fname_prefix)))
		
		alignmentJob = Job(namespace=workflow.namespace, name=stampy.name, version=workflow.version)
		# make sure to use ', rather than ", to wrap the bwaoptions. double-quote(") would disappear during xml translation.
		alignmentJob.addArguments(" --bwa=%s "%(yh_pegasus.getAbsPathOutOfExecutable(bwa)), \
					"--bwaoptions='%s -t%s %s' "%(additionalArguments, no_of_aln_threads, refFastaFile.name),  \
					"-g", refFastaFile, "-h", refFastaFile, "-o", outputSamFile, '--overwrite',\
					'--baq', '--alignquals', '--gatkcigarworkaround')
		#Added option --gatkcigarworkaround removing adjacent I/D events from CIGAR strings, which trips up GATK
		if firstFileObject.db_entry.quality_score_format=='Illumina':
			alignmentJob.addArguments("--solexa")
		alignmentJob.uses(outputSamFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		alignmentJob.output = outputSamFile
		alignmentJob.fname_prefix = fname_prefix
		
		for refFastaFile in refFastaFList:
			alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		yh_pegasus.setJobProperRequirement(alignmentJob, job_max_memory=aln_job_max_memory, \
										no_of_cpus=no_of_aln_threads, max_walltime=aln_job_max_walltime)
		workflow.addJob(alignmentJob)
		#add fastq files after "-M"
		alignmentJob.addArguments(" -M ")	# -M FILE[,FILE] Map fastq file(s).  Use FILE.recaldata for recalibration if available
		for fileObject in fileObjectLs:
			fastqF = fileObject.fastqF
			relativePath = fastqF.name
			alignmentJob.addArguments(fastqF)
			alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
		if refIndexJob:
			self.addRefIndexJobAndItsOutputAsParent(workflow, refIndexJob, childJob=alignmentJob)
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					workflow.depends(parent=mkdirJob, child=alignmentJob)
		return alignmentJob
	
	def addBWAAlnJob(self, workflow=None, executable=None, bwaCommand='aln', fileObject=None, outputFile=None,\
					refFastaFList=None, no_of_aln_threads=3, \
					maxMissingAlignmentFraction=0.04, maxNoOfGaps=1, additionalArguments=None, \
					refIndexJob=None,\
					parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					extraArguments=None, extraArgumentList=None, job_max_memory=2000, \
					key2ObjectForJob=None, max_walltime=None, \
					**keywords):
		"""
		2012.10.10
		"""
		"""
		saiOutput = File(outputFname)
		alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
		alignmentJob.addArguments(alignment_method.command, additionalArguments,"-t %s"%no_of_aln_threads, \
								"-f", saiOutput)
		if fileObject.db_entry.quality_score_format=='Illumina':	#2012.4.5
			alignmentJob.addArguments("-I")
		alignmentJob.addArguments(refFastaFList[0], fastqF)
		alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
		for refFastaFile in refFastaFList:
			alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		alignmentJob.uses(saiOutput, transfer=False, register=True, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(alignmentJob, job_max_memory=aln_job_max_memory, \
										no_of_cpus=no_of_aln_threads, max_walltime=aln_job_max_walltime)
		
		workflow.addJob(alignmentJob)
		"""
		if extraArgumentList is None:
			extraArgumentList = []
		extraArgumentList.append(bwaCommand)
		if additionalArguments:
			extraArgumentList.append(additionalArguments)
		if fileObject.db_entry.quality_score_format=='Illumina':
			extraArgumentList.append("-I")
		if no_of_aln_threads:
			extraArgumentList.append("-t %s"%no_of_aln_threads)
		if maxMissingAlignmentFraction is not None:
			extraArgumentList.append('-n %s'%(maxMissingAlignmentFraction))
		if maxNoOfGaps is not None:
			extraArgumentList.append("-o %s"%(maxNoOfGaps))
		extraArgumentList.extend(["-f", outputFile])
		extraArgumentList.append(refFastaFList[0])
		extraArgumentList.append(fileObject.fastqF)
		
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraDependentInputLs.extend(refFastaFList[:])
		extraDependentInputLs.append(fileObject.fastqF)
		
		if extraOutputLs is None:
			extraOutputLs = []
		extraOutputLs = [outputFile]
		
		alignmentJob = self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
						key2ObjectForJob=None,\
						job_max_memory=job_max_memory, no_of_cpus=no_of_aln_threads, max_walltime=max_walltime)
		if refIndexJob:
			self.addRefIndexJobAndItsOutputAsParent(workflow, refIndexJob, childJob=alignmentJob)
		return alignmentJob
	
	def addBWAAlignmentJob(self, workflow=None, fileObjectLs=None, \
					refFastaFList=None, bwa=None, additionalArguments=None, \
					samtools=None, \
					refIndexJob=None, parentJobLs=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					java=None, SortSamFilesJava=None, SortSamFilesJar=None,\
					addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
					no_of_aln_threads=3, maxMissingAlignmentFraction=0.04, maxNoOfGaps=1, \
					transferOutput=False, **keywords):
		"""
		2012.10.10
		added argument maxMissingAlignmentFraction, maxNoOfGaps
		each fileObject in fileObjectLs should have 2 attributes
			fastqF: registered pegasus file
			db_entry: an individual_sequence_file db_entry or equivalent object. 2 attributes:
				quality_score_format: 'Standard', 'Illumina'
				individual_sequence: 2 attributes:
					sequencer:  'GA' , '454'
					sequence_type: 'PE', 'SR', 'genome'	#(not really used as 2012.10.10)
		alignment_method is also an object (supposedly from vervet db), two attributes.
			short_name (stampy, or bwa-short-read, or bwa-long-read)
			command (stampy.py , aln, bwasw)
		2012.4.5
			handle read files with quality_score_format="Illumina" (post 1.3 solexa).
			input fileObjectLs is different now.
		2012.2.23
			split out of addAlignmentJob()
		2011-9-13
			add argument java & SortSamFilesJar
		2011-9-9
			two steps:
				1. aln doesn't use pipe, outputs to sai files.
				2. sampe/samse, convert, sort => connected through pipe
		"""
		if workflow is None:
			workflow = self
		if namespace is None:
			namespace = workflow.namespace
		if version is None:
			version = workflow.version
		
		aln_job_max_memory = 2600	#in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		bwasw_job_max_memory = 3800	#in MB, same ref, bwasw needs more memory
		samse_job_max_memory = 4500	#in MB
		sampe_job_max_memory = 6000	#in MB
		addRGJob_max_memory = 2500	#in MB
		
		aln_job_max_walltime= 4800	#80 hours, in minutes
		aln_job_max_walltime = 1020	#22 hours, because all reads are stored in chunks of 5-million-read files
		
		memRequirementData = self.getJVMMemRequirment(job_max_memory=addRGJob_max_memory, minMemory=2000)
		job_max_memory = memRequirementData.memRequirement
		javaMemRequirement = memRequirementData.memRequirementInStr
		
		if len(fileObjectLs)==1:	#single end
			fileObject = fileObjectLs[0]
			fastqF = fileObject.fastqF
			relativePath = fastqF.name
			sequence_type = fileObject.db_entry.individual_sequence.sequence_type
			sequencer = fileObject.db_entry.individual_sequence.sequencer
			
			if alignment_method.command=='aln' and sequencer!='454':	#short single-end read
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sai'%fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob = self.addBWAAlnJob(workflow=workflow, executable=bwa, bwaCommand=alignment_method.command, \
									fileObject=fileObject, outputFile=saiOutput,\
					refFastaFList=refFastaFList, no_of_aln_threads=no_of_aln_threads, \
					maxMissingAlignmentFraction=maxMissingAlignmentFraction, maxNoOfGaps=maxNoOfGaps, \
					additionalArguments=additionalArguments, \
					refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
					extraDependentInputLs=None, extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=None, job_max_memory=aln_job_max_memory, \
					key2ObjectForJob=None, max_walltime=aln_job_max_walltime)
				
				alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fname_prefix)))
				sai2samJob = Job(namespace=namespace, name=ShortSEAlignmentByBWA.name, version=version)
				sai2samJob.addArguments(refFastaFList[0], saiOutput, fastqF, alignmentSamF)
				for refFastaFile in refFastaFList:
					sai2samJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				yh_pegasus.setJobProperRequirement(sai2samJob, job_max_memory=samse_job_max_memory)
				workflow.addJob(sai2samJob)
				self.addJobDependency(workflow=workflow, parentJob=alignmentJob, childJob=sai2samJob)
				workflow.no_of_jobs += 1
				
			elif alignment_method.command=='bwasw' or sequencer=='454':	#long single-end read
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fname_prefix)))
				alignmentJob = Job(namespace=namespace, name=LongSEAlignmentByBWA.name, version=version)
				
				alignmentJob.addArguments(refFastaFList[0], fastqF, alignmentSamF, repr(no_of_aln_threads))
				for refFastaFile in refFastaFList:
					alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				yh_pegasus.setJobProperRequirement(alignmentJob, job_max_memory=bwasw_job_max_memory, \
												no_of_cpus=no_of_aln_threads, max_walltime=aln_job_max_walltime)
				workflow.addJob(alignmentJob)
				#fake a sai2samJob
				sai2samJob = alignmentJob
				workflow.no_of_jobs += 1
			
		elif len(fileObjectLs)==2:	#paired end
			fileObject = fileObjectLs[0]
			fastqF1 = fileObject.fastqF
			relativePath = fastqF1.name
			sequence_type = fileObject.db_entry.individual_sequence.sequence_type
			sequencer = fileObject.db_entry.individual_sequence.sequencer
			
			#fastqF1, format, sequence_type = fileObjectLs[0][:3]
			#fastqF2, format, sequence_type = fileObjectLs[1][:3]
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			fname_prefix = fname_prefix[:-2]
			alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fname_prefix)))
			sai2samJob = Job(namespace=namespace, name=PEAlignmentByBWA.name, version=version)
			sai2samJob.addArguments(refFastaFList[0])
			yh_pegasus.setJobProperRequirement(sai2samJob, job_max_memory=sampe_job_max_memory)
			for refFastaFile in refFastaFList:
				sai2samJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
			workflow.addJob(sai2samJob)
			workflow.no_of_jobs += 1
			for fileObject in fileObjectLs:
				fastqF = fileObject.fastqF
				relativePath = fastqF.name
				sequence_type = fileObject.db_entry.individual_sequence.sequence_type
				sequencer = fileObject.db_entry.individual_sequence.sequencer
				
				#fastqF, format, sequence_type = fileObject[:3]
				
				#relativePath, format, sequence_type = fileObject[:3]
				tmp_fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sai'%tmp_fname_prefix)
				saiOutput = File(outputFname)
				
				alignmentJob = self.addBWAAlnJob(workflow=workflow, executable=bwa, bwaCommand=alignment_method.command, \
									fileObject=fileObject, outputFile=saiOutput,\
					refFastaFList=refFastaFList, no_of_aln_threads=no_of_aln_threads, \
					maxMissingAlignmentFraction=maxMissingAlignmentFraction, maxNoOfGaps=maxNoOfGaps, \
					additionalArguments=additionalArguments, \
					refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
					extraDependentInputLs=None, extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=None, job_max_memory=aln_job_max_memory, \
					key2ObjectForJob=None, max_walltime=aln_job_max_walltime)
				
				sai2samJob.addArguments(saiOutput)
				sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=alignmentJob, child=sai2samJob)
			
			#add a pair of fastq files to sampe in the end
			for fileObject in fileObjectLs:
				fastqF = fileObject.fastqF
				sai2samJob.addArguments(fastqF)
				sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
			sai2samJob.addArguments(alignmentSamF)
		sai2samJob.uses(alignmentSamF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		sai2samJob.output = alignmentSamF
		sai2samJob.outputLs = [alignmentSamF]
		sai2samJob.fname_prefix = fname_prefix
		if refIndexJob:
			self.addRefIndexJobAndItsOutputAsParent(workflow, refIndexJob, childJob=sai2samJob)
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					self.addJobDependency(workflow=workflow, parentJob=parentJob, childJob=sai2samJob)
		return sai2samJob
	
	def addAlignmentJob(self, workflow=None, fileObjectLs=None, individual_alignment=None, \
					refFastaFList=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, parentJobLs=None, \
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					java=None, SortSamFilesJava=None, SortSamFilesJar=None,\
					addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
					no_of_aln_threads=3, stampy=None, \
					maxMissingAlignmentFraction=0.04, maxNoOfGaps=1, addBamIndexJob=False,\
					transferOutput=False, **keywords):
		"""
		2012.10.18 add argument addBamIndexJob,
		2012.10.10
		not necessary to add refIndexJob into parentJobLs, because it'll be added as parent job.
		added argument maxMissingAlignmentFraction, maxNoOfGaps for bwa
		each fileObject in fileObjectLs should have 2 attributes
			fastqF: a registered pegasus file
			db_entry: an individual_sequence_file db_entry or equivalent object. 2 attributes:
				quality_score_format: 'Standard', 'Illumina'
				individual_sequence: 2 attributes:
					sequencer:  'GA' , '454'
					sequence_type: 'PE', 'SR', 'genome'	#(not really used as 2012.10.10)
		alignment_method is also an object (supposedly from vervet db), two attributes.
			short_name (stampy, or bwa-short-read, or bwa-long-read)
			command (stampy.py , aln, bwasw)
		2012.9.19
			modify it substantially
		2012.4.5
			choose which alignment program to use based on alignment_method.short_name
		2012.2.23
			split the BWA alignment part to addBWAAlignmentJob()
		2011-9-13
			add argument java & SortSamFilesJar
		2011-9-9
			two steps:
				1. aln doesn't use pipe, outputs to sai files.
				2. sampe/samse, convert, sort => connected through pipe
		"""
		aln_job_max_memory = 2600	#in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		bwasw_job_max_memory = 3800	#in MB, same ref, bwasw needs more memory
		samse_job_max_memory = 4500	#in MB
		sampe_job_max_memory = 6000	#in MB
		addRGJob_max_memory = 2500	#in MB
		
		aln_job_max_walltime= 4800	#80 hours, in minutes		
		
		if alignment_method.short_name=='stampy':
			alignmentJob = self.addStampyAlignmentJob(workflow=workflow, fileObjectLs=fileObjectLs,\
						refFastaFList=refFastaFList, bwa=bwa, \
						additionalArguments=additionalArguments, samtools=samtools, \
						refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
						alignment_method=alignment_method, \
						outputDir=outputDir, namespace=namespace, version=version,\
						PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
						LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
						java=java, SortSamFilesJava=SortSamFilesJava, SortSamFilesJar=SortSamFilesJar,\
						addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
						no_of_aln_threads=no_of_aln_threads, stampy=stampy,\
						transferOutput=transferOutput)
		elif alignment_method.short_name.find('bwa')==0:
			alignmentJob = self.addBWAAlignmentJob(workflow=workflow, fileObjectLs=fileObjectLs, \
						refFastaFList=refFastaFList, bwa=bwa, \
						additionalArguments=additionalArguments, samtools=samtools, \
						refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
						alignment_method=alignment_method, \
						outputDir=outputDir, namespace=namespace, version=version,\
						PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
						LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
						java=java, SortSamFilesJava=SortSamFilesJava, SortSamFilesJar=SortSamFilesJar,\
						addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
						no_of_aln_threads=no_of_aln_threads, transferOutput=transferOutput)
		else:
			sys.stderr.write("Alignment method %s is not supported.\n"%(alignment_method_name))
			sys.exit(3)
		fname_prefix = alignmentJob.fname_prefix
		
		## convert sam into bam	
		bamOutputF = File(os.path.join(outputDir, "%s.bam"%(fname_prefix)))
		sam_convert_job = self.addGenericJob(executable=samtools, inputFile=None,\
							outputFile=bamOutputF, outputArgumentOption="view -bSh -o", \
							parentJobLs=[alignmentJob], extraDependentInputLs=[alignmentJob.output], \
							extraOutputLs=[],\
							transferOutput=transferOutput, \
							extraArgumentList=[alignmentJob.output], \
							job_max_memory=2000)
		
		if individual_alignment:	#if not given then , no read-group addition job
			#2012.9.19 add a AddReadGroup job
			outputRGBAM = File(os.path.join(outputDir, "%s.RG.bam"%(fname_prefix)))
			addRGJob = self.addReadGroupInsertionJob(workflow=workflow, individual_alignment=individual_alignment, \
								inputBamFile=sam_convert_job.output, \
								outputBamFile=outputRGBAM,\
								addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, \
								addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
								parentJobLs=[sam_convert_job], extraDependentInputLs=None, \
								extraArguments=None, job_max_memory = 2500, transferOutput=transferOutput)
			sortAlnParentJob = addRGJob
		else:
			sortAlnParentJob = sam_convert_job
		
		"""
		# 2010-2-4
			sort it so that it could be used for merge
		"""
		bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
		sortBamF = File('%s.bam'%(bam_output_fname_prefix))
		sortAlignmentJob = self.addSortAlignmentJob(workflow=workflow, inputBamFile=sortAlnParentJob.output, \
					outputBamFile=sortBamF,\
					SortSamFilesJava=SortSamFilesJava, SortSamFilesJar=SortSamFilesJar,\
					parentJobLs=[sortAlnParentJob], extraDependentInputLs=None, \
					extraArguments=None, job_max_memory = 2500, transferOutput=transferOutput)
		if addBamIndexJob:
			# add the index job on the merged bam file
			bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
						BuildBamIndexFilesJar=self.BuildBamIndexFilesJar, \
						inputBamF=sortAlignmentJob.output, parentJobLs=[sortAlignmentJob], \
						transferOutput=transferOutput, javaMaxMemory=2000)
			returnJob = bamIndexJob	#bamIndexJob.parentJobLs[0] is sortAlignmentJob.
		else:
			returnJob = sortAlignmentJob
		return returnJob, returnJob.output
	
	def addReadGroupInsertionJob(self, workflow=None, individual_alignment=None, inputBamFile=None, \
								outputBamFile=None,\
								addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
								parentJobLs=None, extraDependentInputLs=None, \
								extraArguments=None, job_max_memory = 2500, transferOutput=False, **keywords):
		"""
		2012.9.19 split out of addAlignmentJob()
		"""
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementData.memRequirement
		javaMemRequirement = memRequirementData.memRequirementInStr
		
		# add RG to this bam
		sequencer = individual_alignment.individual_sequence.sequencer
		read_group = individual_alignment.getReadGroup()	#2012.9.19
		if sequencer=='454':
			platform_id = 'LS454'
		elif sequencer=='GA':
			platform_id = 'ILLUMINA'
		else:
			platform_id = 'ILLUMINA'
		# the add-read-group job
		
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', addOrReplaceReadGroupsJar,\
							"INPUT=", inputBamFile,\
							'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
							'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
							'OUTPUT=', outputBamFile, "VALIDATION_STRINGENCY=LENIENT"]
					#not including 'SORT_ORDER=coordinate'
					#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexFilesJar won't fail.)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputBamFile])
		
		job= self.addGenericJob(executable=addOrReplaceReadGroupsJava, inputFile=None,\
							outputFile=None, outputArgumentOption="-o", \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraOutputLs=[outputBamFile],\
							transferOutput=transferOutput, \
							extraArgumentList=extraArgumentList, \
							job_max_memory=memRequirementData.memRequirement, **keywords)
		return job
	
	def addSortAlignmentJob(self, workflow=None, inputBamFile=None, \
					outputBamFile=None,\
					SortSamFilesJava=None, SortSamFilesJar=None,\
					parentJobLs=None, extraDependentInputLs=None, \
					extraArguments=None, job_max_memory = 2500, transferOutput=False, **keywords):
		"""
		2012.9.19
		"""
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementData.memRequirement
		javaMemRequirement = memRequirementData.memRequirementInStr
		
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', SortSamFilesJar,\
							"SORT_ORDER=coordinate", "I=", inputBamFile, \
							"O=", outputBamFile, "VALIDATION_STRINGENCY=LENIENT"]
					#not including 'SORT_ORDER=coordinate'
					#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexFilesJar won't fail.)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputBamFile])
		
		job= self.addGenericJob(executable=SortSamFilesJava, inputFile=None,\
							outputFile=None, outputArgumentOption="-o", \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraOutputLs=[outputBamFile],\
							transferOutput=transferOutput, \
							extraArgumentList=extraArgumentList, \
							job_max_memory=memRequirementData.memRequirement, **keywords)
		return job
	
	def addAlignmentMergeBySAMtoolsJob(self, workflow, AlignmentJobAndOutputLs=[], outputBamFile=None,samtools=None,\
					java=None, mergeSamFilesJava=None, mergeSamFilesJar=None, \
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					mv=None, transferOutput=False, job_max_memory=2500, **keywords):
		"""
		2012-3.22
			samtools merge version of addAlignmentMergeJob(), which uses picard's MergeSamFiles.jar
			**untested**
		"""
		javaMaxMemory=2500
		if len(AlignmentJobAndOutputLs)>1:
			merge_sam_job = Job(namespace=workflow.namespace, name=samtools.name, version=workflow.version)
			merge_sam_job.addArguments('-f', outputBamFile)		#'overwrite the output BAM if exist'
			merge_sam_job.uses(outputBamFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			yh_pegasus.setJobProperRequirement(merge_sam_job, job_max_memory=job_max_memory)
			workflow.addJob(merge_sam_job)
			for AlignmentJobAndOutput in AlignmentJobAndOutputLs:
				alignmentJob, alignmentOutput = AlignmentJobAndOutput[:2]
				merge_sam_job.addArguments(alignmentOutput)
				merge_sam_job.uses(alignmentOutput, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=alignmentJob, child=merge_sam_job)
		else:	#one input file, no samtools merge. use "mv" to rename it instead
			alignmentJob, alignmentOutput = AlignmentJobAndOutputLs[0][:2]
			merge_sam_job = Job(namespace=workflow.namespace, name=mv.name, version=workflow.version)
			merge_sam_job.addArguments(alignmentOutput, outputBamFile)
			workflow.depends(parent=alignmentJob, child=merge_sam_job)
			merge_sam_job.uses(alignmentOutput, transfer=True, register=True, link=Link.INPUT)
			merge_sam_job.uses(outputBamFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			workflow.addJob(merge_sam_job)
		
		# add the index job on the merged bam file
		bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, \
						BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
						inputBamF=outputBamFile, parentJobLs=[merge_sam_job], \
						transferOutput=transferOutput, javaMaxMemory=javaMaxMemory)
		return merge_sam_job, bamIndexJob
	
	def addMarkDupJob(self, workflow, parentJobLs=[], inputBamF=None, inputBaiF=None, outputBamFile=None,\
					MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir="/Network/Data/vervet/vervetPipeline/tmp/",\
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					namespace='workflow', version='1.0', transferOutput=True, no_of_cpus=1):
		"""
		2012.3.21
			improve it
			set no_of_cpus=1 (was 2) to avoid thread problem in some linux kernels.
		#2011-11-10 duplicate-marking job
		"""
		MarkDupJobMaxMemory=4000
		MarkDupJob = Job(namespace=getattr(workflow, 'namespace', namespace), name=MarkDuplicatesJava.name, \
						version=getattr(workflow, 'version', version))
		bamFnamePrefix = os.path.splitext(outputBamFile.name)[0]
		MarkDupOutputF = outputBamFile
		MarkDupOutputMetricF = '%s.metric'%(bamFnamePrefix)
		
		memRequirementData = self.getJVMMemRequirment(job_max_memory=MarkDupJobMaxMemory, minMemory=2000)
		job_max_memory = memRequirementData.memRequirement
		javaMemRequirement = memRequirementData.memRequirementInStr
		
		MarkDupJob.addArguments(javaMemRequirement, '-jar', MarkDuplicatesJar, "MAX_FILE_HANDLES=500",\
			"VALIDATION_STRINGENCY=LENIENT", "ASSUME_SORTED=true", "INPUT=", inputBamF, \
			'OUTPUT=', MarkDupOutputF, "M=", MarkDupOutputMetricF, "MAX_RECORDS_IN_RAM=500000",\
			"TMP_DIR=%s"%tmpDir)
		MarkDupJob.uses(inputBaiF, transfer=True, register=True, link=Link.INPUT)
		MarkDupJob.uses(inputBamF, transfer=True, register=True, link=Link.INPUT)
		MarkDupJob.output = MarkDupOutputF
		MarkDupJob.MarkDupOutputMetricF = MarkDupOutputMetricF
		MarkDupJob.uses(MarkDupOutputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		MarkDupJob.uses(MarkDupOutputMetricF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		#pass	#don't register the files so leave them there
		workflow.addJob(MarkDupJob)
		yh_pegasus.setJobProperRequirement(MarkDupJob, job_max_memory=job_max_memory, no_of_cpus=no_of_cpus)
		for parentJob in parentJobLs:
			if parentJob:
				workflow.depends(parent=parentJob, child=MarkDupJob)
		
		
		# add the index job on the bam file
		bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, \
								BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
								inputBamF=MarkDupOutputF,\
								parentJobLs=[MarkDupJob], namespace=namespace, version=version,\
								transferOutput=transferOutput)
		return MarkDupJob, bamIndexJob
		
	@classmethod
	def addSAMtoolsCalmdJob(cls, workflow, samtoolsCalmd=None, inputBamF=None, \
					refFastaFList=None, outputBamF=None, \
					parentJob=None, \
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					namespace='workflow', version='1.0', transferOutput=True,\
					**keywords):
		"""
		2011-11-20
			run "samtools calmd" on the input bam and index the output bam
		"""
		job = Job(namespace=namespace, name=BuildBamIndexFilesJava.name, version=version)
		job.addArguments(inputBamF, refFastaFList[0], outputBamF)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=1000)
		job.uses(inputBamF, transfer=True, register=True, link=Link.INPUT)
		for refFastaFile in refFastaFList:
			job.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		if transferOutput:
			job.uses(outputBamF, transfer=True, register=True, link=Link.OUTPUT)
		else:
			pass	#don't register the files so leave them there
		workflow.addJob(job)
		if parentJob:
			workflow.depends(parent=parentJob, child=job)
			
		# add the index job on the bam
		return self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					inputBamF=outputBamF,\
					parentJobLs=[job], namespace=namespace, version=version,\
					transferOutput=transferOutput)