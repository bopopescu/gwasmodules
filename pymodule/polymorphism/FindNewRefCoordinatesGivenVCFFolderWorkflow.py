#!/usr/bin/env python
"""
Examples:
	#2012.10.09 on hoffman condor, no clustering (-C0), always need db connection on hcondor (-H)
	# add -U 0 -Z 3000 if u want to change the partitioning configuration (how many sites in one blast job)
	%s  -I ~/NetworkData/vervet/db/genotype_file/method_42/
		--oldRefFastaFname ~/NetworkData/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
		--newRefFastaFname ~/NetworkData/vervet/db/individual_sequence/3231_6483ContigsVervetRef1.0.3.fasta
		--formatdbPath ~/bin/blast/bin/formatdb --blastallPath ~/bin/blast/bin/blastall
		--maxNoOfMismatches 2 -H -C 0
		-j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
		-u yh -z localhost -o workflow/polymorphism/FindNewRefCoordinates_Method42_vs_3231_maxContigID2.xml -x 2
		#-U 0 -Z 3000
	
	#2012.10.18 use bwa to align
	%s -I ~/NetworkData/vervet/db/genotype_file/method_42/ 
		--oldRefFastaFname ~/NetworkData/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
		--newRefFastaFname ~/NetworkData/vervet/db/individual_sequence/3231_6483ContigsVervetRef1.0.3.fasta
		--maxNoOfMismatches 2 -H -C 4 --no_of_aln_threads 1
		-j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
		-u yh -z localhost -o workflow/polymorphism/FindNewRefCoordinates_Method42_vs_3231_BWA.xml  --alignmentMethodType 2
		--intervalSize 40000
		#--flankingLength 50

Description:
	2012-10-03
	#. preReduce: split the fasta into difference chromosomes (necessary?) 
    #. mapEachVCF:
    		split VCF into several small ones, N-SNP each
      #. mapEachInterval
         #. extract flanking sequences from the input VCF ():
         	ExtractFlankingSequenceForVCFLoci.py -i vcf -a reference fasta -o output.fasta -f flankingLength (=24)
         #. blast them
         #. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
             #. where hit length match query length, and no of mismatches <=2 => good => infer new coordinates
         #. output a mapping file between old SNP and new SNP coordinates.
         		#. reduce this thing by combining everything
    #. reduce()
         #. merge all output into one big one
         
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule.pegasus.BlastWorkflow import BlastWorkflow
from pymodule.pegasus.ShortRead2AlignmentWorkflow import ShortRead2AlignmentWorkflow

class FindNewRefCoordinatesGivenVCFFolderWorkflow(AbstractVCFWorkflow, BlastWorkflow, ShortRead2AlignmentWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update(ShortRead2AlignmentWorkflow.alignment_option_dict.copy())
	option_default_dict.update({
						('oldRefFastaFname', 1, ): ['', '', 1, 'path to the old reference sequence file (on which input VCF is based)', ],\
						("formatdbPath", 1, ): ["/usr/bin/formatdb", 'f', 1, 'path to formatdb, index fasta database file'],\
						("blastallPath", 1, ): ["/usr/bin/blastall", 's', 1, 'path to blastall'],\
						('newRefFastaFname', 1, ): ['', '', 1, 'path to the new reference sequence file (blast db)', ],\
						
						('minNoOfIdentities', 0, int): [None, '', 1, 'minimum number of identities between a query and target', ],\
						('maxNoOfMismatches', 0, int): [None, '', 1, 'minimum number of mismatches between a query and target', ],\
						('minIdentityFraction', 0, float): [None, '', 1, 'minimum percentage of identities between a query and target', ],\
						('flankingLength', 1, int): [24, '', 1, 'number of flanking bases on either side of the locus.\n\
	length of flanking = 2*flankingLength+locusLength', ],\
						('minAlignmentSpan', 1, int): [10, '', 1, 'minimum length of alignment in blast', ],\
						
						('alignmentMethodType', 1, int): [1, '', 1, 'which alignment program to use: 1 blast, 2 bwa', ],\
						('maxMissingAlignmentFraction', 1, float): [0.04, '', 1, ' max fraction of missing alignments given 2% uniform base error rate if FLOAT.', ],\
						('maxNoOfGaps', 1, int): [1, '', 1, 'Maximum number of gap opens', ],\
						
						})
	
	#2012.9.25 no overlap and make the interval a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 10000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
		ShortRead2AlignmentWorkflow.__init__(self, **keywords)
		
		self.oldRefFastaFname = os.path.abspath(self.oldRefFastaFname)
		self.newRefFastaFname = os.path.abspath(self.newRefFastaFname)
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = AbstractVCFWorkflow.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix,\
								passingData=passingData, transferOutput=transferOutput, **keywords)
		
		
		passingData.oldRefFastaFile = self.registerOneInputFile(workflow=workflow, inputFname=self.oldRefFastaFname, \
															folderName=self.pegasusFolderName)
		refIndexJob = None
		if self.alignmentMethodType==1:	#blast
			registerReferenceData = self.registerBlastNucleotideDatabaseFile(self.newRefFastaFname, \
						folderName=self.pegasusFolderName, input_site_handler=self.input_site_handler)
			passingData.newRefFastaFileList = registerReferenceData.refFastaFList
			if len(passingData.newRefFastaFileList)<4:	#some nt-database index file is missing
				sys.stderr.write("Adding blast-db-making job ...")
				refIndexJob = self.addMakeBlastDBJob(executable=self.formatdb,\
											inputFile=passingData.newRefFastaFileList[0], transferOutput=True)
				#add the index files
				passingData.newRefFastaFileList = [passingData.newRefFastaFileList[0]] + refIndexJob.outputList
				sys.stderr.write(".\n")
		else:	#bwa
			registerReferenceData = self.registerBWAIndexFile(refFastaFname=self.newRefFastaFname, \
												folderName=self.pegasusFolderName, input_site_handler=self.input_site_handler)
			passingData.newRefFastaFileList = registerReferenceData.refFastaFList
			if registerReferenceData.needBWARefIndexJob:	#2013.3.20
				sys.stderr.write("Adding bwa reference index job ...")
				refIndexJob = self.addBWAReferenceIndexJob(workflow=workflow, refFastaFList=passingData.newRefFastaFileList, \
					refSequenceBaseCount=3000000000, bwa=workflow.bwa,\
					transferOutput=True)
				#add the index files
				passingData.newRefFastaFileList = [passingData.newRefFastaFileList[0]] + refIndexJob.outputList
				sys.stderr.write(".\n")
				
		
		
		passingData.refIndexJob = refIndexJob
		return returnData
	
	def registerCustomJars(self, workflow=None):
		"""
		2012.10.18
		"""
		#super(FindNewRefCoordinatesGivenVCFFolderWorkflow, self).registerCustomJars(workflow=workflow)
		AbstractVCFWorkflow.registerCustomJars(self, workflow)
		ShortRead2AlignmentWorkflow.registerCustomJars(self, workflow)
		BlastWorkflow.registerCustomJars(self, workflow)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if workflow==None:
			workflow=self
		
		AbstractVCFWorkflow.registerCustomExecutables(self, workflow)
		ShortRead2AlignmentWorkflow.registerCustomExecutables(self, workflow)
		BlastWorkflow.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = workflow.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)		
		
		ExtractFlankingSequenceForVCFLoci = Executable(namespace=namespace, name="ExtractFlankingSequenceForVCFLoci", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		ExtractFlankingSequenceForVCFLoci.addPFN(PFN("file://" + os.path.join(workflow.pymodulePath, \
												"pegasus/mapper/extractor/ExtractFlankingSequenceForVCFLoci.py"), \
										site_handler))
		executableClusterSizeMultiplierList.append((ExtractFlankingSequenceForVCFLoci, 2))
		
		FindSNPPositionOnNewRefFromFlankingBlastOutput = Executable(namespace=namespace, name="FindSNPPositionOnNewRefFromFlankingBlastOutput", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		FindSNPPositionOnNewRefFromFlankingBlastOutput.addPFN(PFN("file://" + os.path.join(workflow.pymodulePath, \
														"polymorphism/FindSNPPositionOnNewRefFromFlankingBlastOutput.py"), \
										site_handler))
		executableClusterSizeMultiplierList.append((FindSNPPositionOnNewRefFromFlankingBlastOutput, 2))
		
		FindSNPPositionOnNewRefFromFlankingBWAOutput = Executable(namespace=namespace, name="FindSNPPositionOnNewRefFromFlankingBWAOutput", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		FindSNPPositionOnNewRefFromFlankingBWAOutput.addPFN(PFN("file://" + os.path.join(workflow.pymodulePath, \
														"polymorphism/FindSNPPositionOnNewRefFromFlankingBWAOutput.py"), \
										site_handler))
		executableClusterSizeMultiplierList.append((FindSNPPositionOnNewRefFromFlankingBWAOutput, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def addFindNewRefCoordinateJob(self, workflow=None, executable=None, inputFile=None, \
							outputFile=None, maxNoOfMismatches=None, minNoOfIdentities=None, minIdentityFraction=None,\
							minAlignmentSpan=None,\
							parentJobLs=None, transferOutput=False, job_max_memory=500,\
							extraArguments=None, extraArgumentList=None, extraDependentInputLs=None):
		"""
		2012.10.14
		
		"""
		if workflow is None:
			workflow = self
		# a FindSNPPositionOnNewRefFromFlankingBlastOutput job
		if extraArgumentList is None:
			extraArgumentList = []
		if minAlignmentSpan is not None:
			extraArgumentList.append('--minAlignmentSpan %s'%(minAlignmentSpan))
		if maxNoOfMismatches is not None:
			extraArgumentList.append('--maxNoOfMismatches %s'%(maxNoOfMismatches))
		if minNoOfIdentities is not None:
			extraArgumentList.append("--minNoOfIdentities %s"%(minNoOfIdentities))
		if minIdentityFraction is not None:
			extraArgumentList.append("--minIdentityFraction %s"%(minIdentityFraction))
			
		findNewRefCoordinateJob = self.addAbstractMapperLikeJob(workflow=workflow, \
					executable=executable, \
					inputF=inputFile, outputF=outputFile, \
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					extraDependentInputLs=extraDependentInputLs)
		return findNewRefCoordinateJob
	
	def mapEachInterval(self, workflow=None, \
					VCFFile=None, passingData=None, transferOutput=False, **keywords):
		"""
		2012.10.3
			#. extract flanking sequences from the input VCF (ref sequence file => contig ref sequence)
			#. blast them
			#. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
				#. where hit length match query length, and no of mismatches <=2 => good => infer new coordinates
			#. output a mapping file between old SNP and new SNP coordinates.
				#. reduce this thing by combining everything
			#. make a new VCF file based on the input split VCF file
				(replace contig ID , position with the new one's, remove the header part regarding chromosomes or replace it)

		"""
		if workflow is None:
			workflow = self
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []

		topOutputDirJob = passingData.topOutputDirJob
		mapDirJob = passingData.mapDirJob
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		intervalFnamePrefix = passingData.intervalFnamePrefix
		jobData = passingData.jobData
		
		
		splitVCFJob = passingData.mapEachVCFData.splitVCFJob
		chr = passingData.chromosome
		
		# a flanking sequence extraction job
		
		outputFnamePrefix = os.path.join(mapDirJob.output, '%s_flankSequence'%(intervalFnamePrefix)) 
		outputFile = File('%s.fasta'%(outputFnamePrefix))
		oldRefFastaFile = passingData.oldRefFastaFile
		extraArgumentList = ['--flankingLength %s'%(self.flankingLength), '--refFastaFname', oldRefFastaFile,\
							"--outputFormatType %s"%(self.alignmentMethodType)]
		# alignmentMethodType 1 => fasta
		# alignmentMethodType 2 => fastq
		extractFlankSeqJob = self.addAbstractMapperLikeJob(workflow=workflow, executable=self.ExtractFlankingSequenceForVCFLoci, \
					inputF=VCFFile, outputF=outputFile, \
					parentJobLs=[mapDirJob, splitVCFJob], transferOutput=transferOutput, job_max_memory=2000,\
					extraArguments=None, extraArgumentList=extraArgumentList, extraDependentInputLs=[oldRefFastaFile])
		
		if self.alignmentMethodType==1:	#blast
			# a blast job
			blastOutputFnamePrefix = '%s_blast'%(outputFnamePrefix)
			outputFile = File('%s.tsv'%(blastOutputFnamePrefix))
			newRefFastaFile = passingData.newRefFastaFileList[0]
			blastJob = self.addBlastWrapperJob(executable=self.BlastWrapper, inputFile=extractFlankSeqJob.output, \
											outputFile=outputFile, \
							outputFnamePrefix=blastOutputFnamePrefix, databaseFile=newRefFastaFile,\
							maxNoOfMismatches=self.maxNoOfMismatches, minNoOfIdentities=self.minNoOfIdentities, \
							minIdentityPercentage=self.minIdentityFraction, blastallPath=self.blastallPath, \
							parentJobLs=[mapDirJob, extractFlankSeqJob, passingData.refIndexJob], \
							extraDependentInputLs=passingData.newRefFastaFileList, \
							transferOutput=transferOutput, \
							extraArguments=None, job_max_memory=5000)
			
			# a FindSNPPositionOnNewRefFromFlankingBlastOutput job
			outputFile = File('%s_2NewRefCoordinateMap.tsv'%(outputFnamePrefix))
			findNewRefCoordinateJob = self.addFindNewRefCoordinateJob(workflow=workflow, executable=self.FindSNPPositionOnNewRefFromFlankingBlastOutput, \
						inputFile=blastJob.output, outputFile=outputFile, maxNoOfMismatches=self.maxNoOfMismatches, \
						minNoOfIdentities=self.minNoOfIdentities, minIdentityFraction=self.minIdentityFraction, \
						minAlignmentSpan=self.minAlignmentSpan, parentJobLs=[blastJob], \
						transferOutput=transferOutput, job_max_memory=500, extraArguments=None, \
						extraArgumentList=None)
			
		else:	#bwa
			#fake a fileObjectLs, which contains a single-end read fastq object
			# db_entry is supposedly a individual_sequence_file object 
			fileObjectLs = [PassingData(fastqF=extractFlankSeqJob.output, \
									db_entry=PassingData(quality_score_format='Standard', \
														individual_sequence=PassingData(sequencer='GA', sequence_type='SR')))]
			#2012.10.10 fake one
			alignment_method = PassingData(short_name='bwa-short-read', command='aln')
			
			#2012.10.10 individual_alignment is not passed so that ReadGroup addition job is not added in addAlignmentJob()
			bamIndexJob = self.addAlignmentJob(workflow=workflow, fileObjectLs=fileObjectLs, \
				refFastaFList=passingData.newRefFastaFileList, bwa=workflow.bwa, \
				additionalArguments=self.additionalArguments, samtools=workflow.samtools, \
				refIndexJob=passingData.refIndexJob, parentJobLs=[mapDirJob, extractFlankSeqJob], \
				alignment_method=alignment_method, \
				outputDir=mapDirJob.output,\
				PEAlignmentByBWA=workflow.PEAlignmentByBWA, ShortSEAlignmentByBWA=workflow.ShortSEAlignmentByBWA, \
				LongSEAlignmentByBWA=workflow.LongSEAlignmentByBWA,\
				java=workflow.java, SortSamFilesJava=workflow.SortSamFilesJava, SortSamJar=workflow.SortSamJar,\
				addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, \
				AddOrReplaceReadGroupsJar=workflow.AddOrReplaceReadGroupsJar,\
				no_of_aln_threads=self.no_of_aln_threads,\
				stampy=workflow.stampy, \
				maxMissingAlignmentFraction=self.maxMissingAlignmentFraction, maxNoOfGaps=self.maxNoOfGaps, \
				addBamIndexJob=True, transferOutput = False)[0]
			
			alignmentJob = bamIndexJob.parentJobLs[0]
			# a FindSNPPositionOnNewRefFromFlankingBlastOutput job
			outputFile = File('%s_2NewRefCoordinateMap.tsv'%(outputFnamePrefix))
			findNewRefCoordinateJob = self.addFindNewRefCoordinateJob(workflow=workflow, \
						executable=self.FindSNPPositionOnNewRefFromFlankingBWAOutput, \
						inputFile=alignmentJob.output, outputFile=outputFile, maxNoOfMismatches=self.maxNoOfMismatches, \
						minNoOfIdentities=self.minNoOfIdentities, minIdentityFraction=self.minIdentityFraction, \
						minAlignmentSpan=self.minAlignmentSpan, parentJobLs=[alignmentJob, bamIndexJob], \
						transferOutput=transferOutput, job_max_memory=500, extraArguments=None, \
						extraArgumentList=None, extraDependentInputLs=[bamIndexJob.output])
			

		returnData.findNewRefCoordinateJob = findNewRefCoordinateJob
		return returnData
	
	def reduce(self, workflow=None, passingData=None, reduceEachChromosomeDataLs=None, transferOutput=True, **keywords):
		"""
		2012.10.3
			#. merge all output of input jobs (passingData.mapEachIntervalDataLsLs) into one big one
		
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		outputFile = File(os.path.join(reduceOutputDirJob.output, 'SNPID2NewCoordinates.tsv'))
		reduceJob = self.addStatMergeJob(workflow, \
									statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
									outputF=outputFile, \
									parentJobLs=[reduceOutputDirJob],extraOutputLs=[], \
									extraDependentInputLs=[], transferOutput=transferOutput,)
		returnData.jobDataLs.append(PassingData(jobLs=[reduceJob], file=reduceJob.output, \
											fileList=[reduceJob.output]))
		
		for mapEachIntervalDataLs in passingData.mapEachIntervalDataLsLs:
			for mapEachIntervalData in mapEachIntervalDataLs:
				parentJob = mapEachIntervalData.findNewRefCoordinateJob
				self.addInputToStatMergeJob(workflow, statMergeJob=reduceJob, \
						parentJobLs=[parentJob])
		
		return returnData

if __name__ == '__main__':
	main_class = FindNewRefCoordinatesGivenVCFFolderWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()