#!/usr/bin/env python
"""
Examples:
	#
	%s  ...
	
	#2013.11.24
	%s --inputDir LiftPolymorphismCoordinates/FindNewRefCoordinates_Method109_vs_3488_BWA_F99.2013.Jul.11T191341/folderReduceLiftOverVCF/
		-H -C 10 -j hcondor -l hcondor -D /u/home/p/polyacti/NetworkData/vervet/db/ -t /u/home/p/polyacti/NetworkData/vervet/db/
		-o dags/SameSiteConcordance/Method109_vs_3488_BWA_F99.sameSiteConcordance.xml --notToUseDBToInferVCFNoOfLoci
		--db_user yh -z localhost

Description:
	2013.11.24 a generic workflow that map-reduces inputs of one or multiple genomic files (i.e. multi-chromosome, tabix-indexed )
		parent class is AbstractNGSWorkflow.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import Executable, File, PFN
from pymodule import ProcessOptions, PassingData, utils
from pymodule.pegasus import yh_pegasus
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow
from pymodule.yhio.FastaFile import FastaFile

parentClass = AbstractNGSWorkflow

class MapReduceGenomeFileWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						})
	
	#2013.11.22 default overlap is 5Mb, overlap is 500K
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 500000
	option_default_dict[('intervalSize', 1, int)][0] = 5000000
	option_default_dict[('max_walltime', 1, int)][0] = 1300	#under 23 hours
	option_default_dict[('inputSuffixList', 0, )][0] = None	#do not set input suffix
	
	def __init__(self,  **keywords):
		"""
		2013.11.24
		"""
		
		self.needSplitChrIntervalData = True
		parentClass.__init__(self, **keywords)
		self.needSplitChrIntervalData = True
		
		self.mapReduceType = 1	#2013.06.27 type 1: split VCF with fixed number of sites
		# type 2: SelectVariants from VCF with fixed-size windows
		# child classes could change its value in the end of their own __init__()
		
		#2013.07.18 offer child classes option to turn it off
		self.needGzipPreReduceReturnData = True
		self.needGzipReduceReturnData = True
	
	def connectDB(self):
		"""
		2012.9.24
			place holder.
		"""
		parentClass.connectDB(self)
		
		self.registerReferenceData = None
		self.refFastaFList= None
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2013.06.14
			move topOutputDirJob from addAllJobs to here. 
		2012.9.17
		"""
		if workflow is None:
			workflow = self
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		self.topOutputDirJob = self.addMkDirJob(outputDir="%sRun"%(outputDirPrefix))
		passingData.topOutputDirJob = self.topOutputDirJob
		
		mapDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
										outputDir="%sMap"%(outputDirPrefix))
		passingData.mapDirJob = mapDirJob
		returnData.mapDirJob = mapDirJob
		self.mapDirJob = mapDirJob
		
		reduceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
												outputDir="%sReduce"%(outputDirPrefix))
		passingData.reduceOutputDirJob = reduceOutputDirJob
		returnData.reduceOutputDirJob = reduceOutputDirJob
		
		self.plotDirJob = self.addMkDirJob(outputDir="%sPlot"%(outputDirPrefix))
		self.statDirJob = self.addMkDirJob(outputDir="%sStat"%(outputDirPrefix))
		self.reduceStatDirJob = self.addMkDirJob(outputDir="%sReduceStat"%(outputDirPrefix))
		self.reduceEachInputDirJob = self.addMkDirJob(outputDir="%sReduceEachInput"%(outputDirPrefix))
		self.reduceEachChromosomeDirJob = self.addMkDirJob(outputDir="%sReduceEachChromosome"%(outputDirPrefix))
		self.reduceOutputDirJob = reduceOutputDirJob
		return returnData
	
	def selectIntervalFromInputFile(self, jobData=None, chromosome=None,\
								intervalData=None, mapEachChromosomeData=None,\
								passingData=None, transferOutput=False,\
								**keywords):
		"""
		2013.11.24
		"""
		inputSuffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(jobData.file.name)[1]
		outputFile = File(os.path.join(self.mapDirJob.output, '%s_%s%s'%(passingData.fileBasenamePrefix, \
																intervalData.overlapInterval, inputSuffix)))
		tabixRetrieveJob = self.addTabixRetrieveJob(executable=self.tabixRetrieve, \
							tabixPath=self.tabixPath, \
							inputF=jobData.file, outputF=outputFile, \
							regionOfInterest=intervalData.overlapInterval, includeHeader=True,\
							parentJobLs=jobData.jobLs + [self.mapDirJob], job_max_memory=100, \
							extraDependentInputLs=jobData.fileLs[1:], \
							transferOutput=False)
		return self.constructJobDataFromJob(job=tabixRetrieveJob)
		

	def addAllJobs(self, workflow=None, inputData=None, chr2IntervalDataLs=None, \
				data_dir=None, \
				intervalSize=3000, intervalOverlapSize=0, \
				outputDirPrefix="", passingData=None, \
				transferOutput=True, job_max_memory=2000, **keywords):
		"""
		2013.06.14 bugfix regarding noOfUnits, which was all inferred from one file
		2012.7.26
			architect of the whole map-reduce framework
		"""
		sys.stderr.write("Adding jobs for %s input genome files... \n"%(len(inputData.jobDataLs)) )
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		#2012.9.22 
		# 	mapEachAlignmentDataLs is never reset.
		#	mapEachChromosomeDataLs is reset upon new alignment
		#	mapEachIntervalDataLs is reset upon each new chromosome
		#	all reduce lists never get reset.
		#	fileBasenamePrefix is the prefix of input file's basename, to be used for temporary output files in reduceEachInput()
		#		but not for output files in mapEachInterval()
		passingData = PassingData(\
					fileBasenamePrefix=None, \
					chromosome=None, \
					
					outputDirPrefix=outputDirPrefix, \
					intervalFileBasenamePrefix=None,\
					
					registerReferenceData=None, \
					refFastaFList=None, \
					refFastaF=None,\
					
					fastaDictJob = None,\
					refFastaDictF = None,\
					fastaIndexJob = None,\
					refFastaIndexF = None,\
					
					intervalOverlapSize =intervalOverlapSize, intervalSize=intervalSize,\
					jobData=None,\
					splitInputFile=None,\
					intervalDataLs=None,\
					preReduceReturnData=None,\
					
					mapEachIntervalData=None,\
					mapEachIntervalDataLs=None,\
					mapEachIntervalDataLsLs=[],\
					mapEachInputData=None,\
					mapEachInputDataLs=None,\
					mapEachInputDataLsLs=[],\
					mapEachChromosomeData=None, \
					mapEachChromosomeDataLs=[], \
					
					chromosome2mapEachIntervalDataLs = {},\
					chromosome2mapEachInputDataLs = {},\
					
					reduceEachInputData=None,\
					reduceEachChromosomeData=None,\
					reduceEachInputDataLs=None,\
					reduceEachInputDataLsLs=[],\
					reduceEachChromosomeDataLs=[],\
					)
		# mapEachIntervalDataLsLs is list of mapEachIntervalDataLs by each Input file.
		# mapEachInputDataLsLs is list of mapEachInputDataLs by each chromosome
		# reduceEachInputDataLsLs is list of reduceEachInputDataLs by each chromosome
		
		preReduceReturnData = self.preReduce(workflow=workflow, outputDirPrefix=outputDirPrefix, \
									passingData=passingData, transferOutput=True,\
									**keywords)
		passingData.preReduceReturnData = preReduceReturnData
		
		#gzip folder jobs (to avoid repeatedly creating the same folder
		gzipReduceEachInputFolderJob = None
		gzipReduceEachChromosomeFolderJob = None
		gzipReduceFolderJob = None
		gzipPreReduceFolderJob = None
		no_of_input_files = 0
		
		firstInterval = True
		
		for chromosome, intervalDataLs in chr2IntervalDataLs.iteritems():
			passingData.chromosome = chromosome
			mapEachChromosomeData = self.mapEachChromosome(workflow=workflow, chromosome=chromosome, \
										passingData=passingData, \
										transferOutput=False, **keywords)
			passingData.mapEachChromosomeData = mapEachChromosomeData
			passingData.mapEachChromosomeDataLs.append(mapEachChromosomeData)
			
			passingData.mapEachInputDataLsLs.append([])
			#the last one from the double list is the current one
			passingData.mapEachInputDataLs = passingData.mapEachInputDataLsLs[-1]
			passingData.mapEachIntervalDataLs = []
			passingData.chromosome2mapEachIntervalDataLs[chromosome] = []
			
			passingData.reduceEachInputDataLsLs.append([])
			passingData.reduceEachInputDataLs = passingData.reduceEachInputDataLsLs[-1]
			
			for i in xrange(len(inputData.jobDataLs)):
				jobData = inputData.jobDataLs[i]
				passingData.jobData = jobData
				passingData.inputJobData = jobData
				
				InputFile = jobData.file
				commonFileBasenamePrefix = utils.getFileBasenamePrefixFromPath(InputFile.name)
				passingData.fileBasenamePrefix = commonFileBasenamePrefix
				
				no_of_input_files += 1
				if no_of_input_files%10==0:
					sys.stderr.write("%s\t%s Inputs."%('\x08'*40, no_of_input_files))
				
				for intervalData in intervalDataLs:
					selectIntervalJobData = self.selectIntervalFromInputFile(jobData=jobData, chromosome=chromosome,\
													intervalData=intervalData, mapEachChromosomeData=mapEachChromosomeData,\
													passingData=passingData, transferOutput=firstInterval,\
													**keywords)
					mapEachIntervalData = self.mapEachInterval(workflow=workflow, inputJobData=jobData, \
															selectIntervalJobData=selectIntervalJobData, \
										chromosome=chromosome,intervalData=intervalData,\
										mapEachChromosomeData=mapEachChromosomeData, \
										passingData=passingData, transferOutput=firstInterval, **keywords)
					
					passingData.mapEachIntervalData = mapEachIntervalData
					passingData.mapEachIntervalDataLs.append(mapEachIntervalData)
					passingData.chromosome2mapEachIntervalDataLs[chromosome].append(mapEachIntervalData)
					
					linkMapToReduceData = self.linkMapToReduce(workflow=workflow, mapEachIntervalData=mapEachIntervalData, \
										preReduceReturnData=preReduceReturnData, \
										passingData=passingData, \
										**keywords)
					if firstInterval==True:
						firstInterval = False
				reduceEachInputData = self.reduceEachInput(workflow=workflow, chromosome=chromosome, passingData=passingData, \
								mapEachIntervalDataLs=passingData.mapEachIntervalDataLs,\
								transferOutput=False, data_dir=data_dir, \
								**keywords)
				passingData.reduceEachInputData = reduceEachInputData
				passingData.reduceEachInputDataLs.append(reduceEachInputData)
				
				gzipReduceEachInputData = self.addGzipSubWorkflow(workflow=workflow, \
					inputData=reduceEachInputData, transferOutput=transferOutput,\
					outputDirPrefix="%sReduceEachInput"%(outputDirPrefix), topOutputDirJob=gzipReduceEachInputFolderJob, \
					report=False)
				gzipReduceEachInputFolderJob = gzipReduceEachInputData.topOutputDirJob
			reduceEachChromosomeData = self.reduceEachChromosome(workflow=workflow, chromosome=chromosome, passingData=passingData, \
								mapEachInputDataLs=passingData.mapEachInputDataLs, \
								chromosome2mapEachIntervalDataLs=passingData.chromosome2mapEachIntervalDataLs,\
								reduceEachInputDataLs=passingData.reduceEachInputDataLs,\
								transferOutput=False, data_dir=data_dir, \
								**keywords)
			passingData.reduceEachChromosomeData = reduceEachChromosomeData
			passingData.reduceEachChromosomeDataLs.append(reduceEachChromosomeData)
			
			gzipReduceEachChromosomeData = self.addGzipSubWorkflow(workflow=workflow, \
					inputData=reduceEachChromosomeData, transferOutput=transferOutput,\
					outputDirPrefix="%sReduceEachChromosome"%(outputDirPrefix), \
					topOutputDirJob=gzipReduceEachChromosomeFolderJob, report=False)
			gzipReduceEachChromosomeFolderJob = gzipReduceEachChromosomeData.topOutputDirJob
			
		reduceReturnData = self.reduce(workflow=workflow, passingData=passingData, transferOutput=False, \
							mapEachChromosomeDataLs=passingData.mapEachInputDataLs,\
							reduceEachChromosomeDataLs=passingData.reduceEachChromosomeDataLs,\
							**keywords)
		passingData.reduceReturnData = reduceReturnData
		
		if self.needGzipPreReduceReturnData:
			gzipPreReduceReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=preReduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%sPreReduce"%(outputDirPrefix), \
						topOutputDirJob= gzipPreReduceFolderJob, report=False)
			gzipPreReduceFolderJob = gzipPreReduceReturnData.topOutputDirJob
		
		if self.needGzipReduceReturnData:
			gzipReduceReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=reduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%sReduce"%(outputDirPrefix), \
						topOutputDirJob=gzipReduceFolderJob, report=False)
			gzipReduceFolderJob = gzipReduceReturnData.topOutputDirJob
		
		sys.stderr.write("\n %s%s Input files.\n"%('\x08'*40, no_of_input_files))
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return reduceReturnData

	
	def mapEachChromosome(self, workflow=None, chromosome=None,\
				VCFJobData=None, jobData=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachInterval(self, workflow=None, inputJobData=None, selectIntervalJobData=None, \
					chromosome=None,intervalData=None,\
					mapEachChromosomeData=None, \
					passingData=None, transferOutput=False, **keywords):
		"""
		2013.04.08 use inputJobData
		2012.10.3
			#. extract flanking sequences from the input Input (ref sequence file => contig ref sequence)
			#. blast them
			#. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
				#. where hit length match query length, and no of mismatches <=2 => good => infer new coordinates
			#. output a mapping file between old SNP and new SNP coordinates.
				#. reduce this thing by combining everything
			#. make a new Input file based on the input split Input file
				(replace contig ID , position with the new one's, remove the header part regarding chromosomes or replace it)

		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		#passingData.intervalFileBasenamePrefix
		#passingData.splitInputFile
		#passingData.unitNumber
		"""
		## 2013.06.19 structures available from passingData, specific to the interval
		passingData.splitInputFile = splitInputFile
		passingData.unitNumber = unitNumber
		passingData.intervalFileBasenamePrefix = '%s_%s_splitInput_u%s'%(chromosome, commonPrefix, unitNumber)
		passingData.noOfIndividuals = jobData.file.noOfIndividuals
		passingData.span = self.intervalSize + self.intervalOverlapSize*2 	#2013.06.19 for memory/walltime gauging
		"""
		return returnData

	def linkMapToReduce(self, workflow=None, mapEachIntervalData=None, preReduceReturnData=None, passingData=None, transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def reduceEachChromosome(self, workflow=None, chromosome=None, passingData=None, mapEachInputDataLs=None, \
						chromosome2mapEachIntervalDataLs=None,\
						reduceEachInputDataLs=None,\
						transferOutput=True, \
						**keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachInputDataLs = mapEachInputDataLs
		returnData.reduceEachInputDataLs = reduceEachInputDataLs
		return returnData
	
	def reduceEachInput(self, workflow=None, chromosome=None, passingData=None, mapEachIntervalDataLs=None,\
					transferOutput=True, **keywords):
		"""
		2013.07.10
			#. concatenate all the sub-Inputs into one
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
		
		#intervalJobLs = [pdata for pdata in mapEachIntervalDataLs]
		
		"""
		realInputVolume = passingData.jobData.file.noOfIndividuals * passingData.jobData.file.noOfLoci
		baseInputVolume = 200*20000
		walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value
		job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=5000, \
							minJobPropertyValue=5000, maxJobPropertyValue=10000).value
		"""
		return returnData
	
	def reduce(self, workflow=None, reduceEachChromosomeDataLs=None, \
			mapEachChromosomeDataLs=None, passingData=None, transferOutput=True, \
			**keywords):
		"""
		2013.07.18 return each processed-Input job data so that followup workflows could carry out map-reduce
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachChromosomeDataLs = mapEachChromosomeDataLs
		returnData.reduceEachChromosomeDataLs = reduceEachChromosomeDataLs
		"""
		#2013.07.18 example to return each processed-Input job data so that followup workflows could carry out map-reduce
		for reduceEachInputDataLs in passingData.reduceEachInputDataLsLs:
			if reduceEachInputDataLs:
				for reduceEachInputData in reduceEachInputDataLs:
					if reduceEachInputData:
						returnData.jobDataLs.append(reduceEachInputData.WHATEVERJobData)
		"""
		return returnData
	
	def setup_run(self):
		"""
		2013.06.11 added firstInputJobData in return
			assign all returned data to self, rather than pdata (pdata has become self)
		2013.04.07 wrap all standard pre-run() related functions into this function.
			setting up for run(), called by run()
		"""
		pdata = parentClass.setup_run(self)
		
		#self.chr2size = {}
		#self.chr2size = set(['Contig149'])	#temporary when testing Contig149
		#self.chr2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		if self.needSplitChrIntervalData:	#2013.06.21 defined in parentClass.__init__()
			chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=self.chr2size, \
													intervalSize=self.intervalSize, \
													intervalOverlapSize=self.intervalOverlapSize)
		else:
			chr2IntervalDataLs = None
		inputData = None
		firstInputJobData = None
		if getattr(self, 'inputDir', None):
			inputData = self.registerAllInputFiles(inputDir=self.inputDir, input_site_handler=self.input_site_handler, \
											pegasusFolderName=self.pegasusFolderName,\
											inputSuffixSet=self.inputSuffixSet,\
											indexFileSuffixSet=set(['.tbi', '.fai']),\
											checkEmptyInputByReading=getattr(self, 'checkEmptyInputByReading', None),\
											maxContigID=self.maxContigID, \
											minContigID=self.minContigID,\
											db_vervet=getattr(self, 'db_vervet', None), \
											needToKnowNoOfLoci=getattr(self, 'needToKnowNoOfLoci', True),\
											minNoOfLociInInput=getattr(self, 'minNoOfLociInInput', 10))
			if inputData and inputData.jobDataLs:
				firstInputJobData = inputData.jobDataLs[0]
				#job=None, jobLs=[], vcfFile=inputF, tbi_F=tbi_F, file=inputF, fileLs=[inputF, tbi_F]
				firstInputFile = firstInputJobData.file
				sys.stderr.write("\t Input file %s is chosen as an example Input for any job that needs a random Input file.\n"%(firstInputFile))
		
		registerReferenceData = self.getReferenceSequence()
		
		self.inputData = inputData
		self.chr2IntervalDataLs = chr2IntervalDataLs
		self.registerReferenceData = registerReferenceData
		self.firstInputJobData = firstInputJobData
		#self.firstInputFile = firstInputFile
		return self
	
	def run(self):
		"""
		2012.9.24
		"""
		pdata = self.setup_run()
		workflow = pdata.workflow
		
		
		inputData=pdata.inputData
		
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("No VCF files in this folder , %s.\n"%self.inputDir)
			sys.exit(0)
				
		self.addAllJobs(inputData=inputData, \
				chr2IntervalDataLs=self.chr2IntervalDataLs, \
				data_dir=self.data_dir, \
				intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
				outputDirPrefix=self.pegasusFolderName, transferOutput=True,)
		
		self.end_run()
	

if __name__ == '__main__':
	main_class = MapReduceGenomeFileWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()