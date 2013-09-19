#!/usr/bin/env python
"""
2012.1.17
	a common class for pegasus workflows that work on VCF variant files
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, utils
from pymodule.yhio.MatrixFile import MatrixFile
from pymodule.yhio.VCFFile import VCFFile
from AbstractNGSWorkflow import AbstractNGSWorkflow

class AbstractVCFWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractNGSWorkflow.option_default_dict)
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('minDepth', 0, float): [0, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						
						("ligateVcfPerlPath", 1, ): ["%s/bin/umake/scripts/ligateVcf.pl", '', 1, 'path to ligateVcf.pl'],\
						("notToKnowNoOfLoci", 0, int): [0, '', 0, 'By default, this program infers the number of loci for each VCF file (for splitting, etc.) \n\
	either using first number of filename as db id or reading through the file. Toggle this to disable it.'],\
						("notToUseDBToInferVCFNoOfLoci", 0, int): [0, '', 0, 'assuming --notToKnowNoOfLoci is off, \n\
	toggle this to infer the number of loci by strictly reading through the file'],\
						('minNoOfLociInVCF', 0, int): [5, '', 1, 'minimum number of loci for an input VCF file to be included.', ],\
						})
	#update these two as they mean the number of loci now, not base pairs
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 100
	option_default_dict[('intervalSize', 1, int)][0] = 10000
	
	
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		
		self.mapReduceType = 1	#2013.06.27 type 1: split VCF with fixed number of sites
		# type 2: SelectVariants from VCF with fixed-size windows
		# child classes could change its value in the end of their own __init__()
		
		#2013.07.18 offer child classes option to turn it off
		self.needGzipPreReduceReturnData = True
		self.needGzipReduceReturnData = True
		
	def extra__init__(self):
		"""
		2013.2.14
		"""
		AbstractNGSWorkflow.extra__init__(self)
		if getattr(self, "inputDir", None):
			self.inputDir = os.path.abspath(self.inputDir)
		
		if hasattr(self, 'ligateVcfPerlPath'):
			self.ligateVcfPerlPath =  self.insertHomePath(self.ligateVcfPerlPath, self.home_path)
		
	
	def addAddVCFFile2DBJob(self, executable=None, inputFile=None, genotypeMethodShortName=None,\
						logFile=None, format=None, data_dir=None, checkEmptyVCFByReading=None, commit=False, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.12.12 use extended argument name
		2012.10.6 use addGenericDBJob() instead of addGenericJob()
		2012.8.30 moved from vervet/src/AddVCFFolder2DBWorkflow.py
		2012.6.27
		"""
		extraArgumentList = ['--format', format]
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
		if data_dir:
			extraArgumentList.extend(['--data_dir', data_dir])
		if checkEmptyVCFByReading:
			extraArgumentList.extend(['--checkEmptyVCFByReading'])
		if genotypeMethodShortName:
			extraArgumentList.extend(['--genotypeMethodShortName', genotypeMethodShortName, ])
		if commit:
			extraArgumentList.append('--commit')
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		job= self.addGenericDBJob(executable=executable, inputFile=inputFile, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[logFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		return job
	
	def registerExecutables(self, workflow=None):
		"""
		"""
		if not workflow:
			workflow = self
		AbstractNGSWorkflow.registerExecutables(self, workflow=workflow)
		
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		#2012.8.30 moved from vervet/src/AddVCFFolder2DBWorkflow.py
		AddVCFFile2DB = Executable(namespace=namespace, name="AddVCFFile2DB", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		AddVCFFile2DB.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/input/AddVCFFile2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((AddVCFFile2DB, 1))
		
		FilterVCFSNPCluster = Executable(namespace=namespace, name="FilterVCFSNPCluster", version=version, os=operatingSystem, arch=architecture, installed=True)
		FilterVCFSNPCluster.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "pegasus/mapper/filter/FilterVCFSNPCluster.py"), site_handler))
		executableClusterSizeMultiplierList.append((FilterVCFSNPCluster, 1))
		
		JuxtaposeAlleleFrequencyFromMultiVCFInput = Executable(namespace=namespace, name="JuxtaposeAlleleFrequencyFromMultiVCFInput", \
											version=version, os=operatingSystem, arch=architecture, installed=True)
		JuxtaposeAlleleFrequencyFromMultiVCFInput.addPFN(PFN("file://" + os.path.join(self.pymodulePath, \
											"pegasus/mapper/extractor/JuxtaposeAlleleFrequencyFromMultiVCFInput.py"), \
											site_handler))
		executableClusterSizeMultiplierList.append((JuxtaposeAlleleFrequencyFromMultiVCFInput, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		
		#2013.07.12
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='SelectVariantsJavaInReduce', \
															clusterSizeMultipler=0.001)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "polymorphism/qc/RemoveRedundantLociFromVCF.py"), \
									name='RemoveRedundantLociFromVCF_InReduce', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "polymorphism/qc/RemoveRedundantLociFromVCF.py"), \
									name='RemoveRedundantLociFromVCF', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "polymorphism/qc/ClearVCFBasedOnSwitchDensity.py"), \
									name='ClearVCFBasedOnSwitchDensity', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "polymorphism/qc/CalculateSameSiteConcordanceInVCF.py"), \
									name='CalculateSameSiteConcordanceInVCF', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/extractor/ExtractInfoFromVCF.py"), \
									name='ExtractInfoFromVCF', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "mapper/ExtractSamplesFromVCF.py"), \
									name='ExtractSamplesFromVCF', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "shell/ligateVcf.sh"), \
									name="ligateVcf", clusterSizeMultipler=1)
	
	def registerCommonExecutables(self, workflow=None):
		"""
		"""
		AbstractNGSWorkflow.registerCommonExecutables(self, workflow=workflow)
	
	def registerAllInputFiles(self, workflow=None, inputDir=None, input_site_handler=None, \
					checkEmptyVCFByReading=False, pegasusFolderName='',\
					maxContigID=None, minContigID=None, db_vervet=None, needToKnowNoOfLoci=False,
					minNoOfLociInVCF=None, includeIndelVCF=True, notToUseDBToInferVCFNoOfLoci=None):
		"""
		2013.07.17 bugfix
		2013.3.1 flip includeIndelVCF to true (now indel and SNP vcf files from AlignmentToCall workflows are in separate folders.
		2012.8.15 add argument db_vervet, needToKnowNoOfLoci, to get noOfLoci by parsing inputFname and find db-entry...
			argument minNoOfLociInVCF, only used when it's not None and needToKnowNoOfLoci is True
		2012.8.10 add maxContigID and minContigID to restrict input
		2012.7.27 add attribute file to each object in returnData.jobDataLs
		2012.5.9
			register the tbi file if it exists
		2012.3.1
			moved from CalculateTrioInconsistencyPipeline.py
		2012.1.9
			the returning data structure is changed to conform to some standard used across several functions
		2011-9-29
			vcf files only
		"""
		if workflow is None:
			workflow = self
		if input_site_handler is None:
			input_site_handler = workflow.input_site_handler
		if notToUseDBToInferVCFNoOfLoci:
			notToUseDBToInferVCFNoOfLoci = self.notToUseDBToInferVCFNoOfLoci
		
		sys.stderr.write("Registering input files from %s, maxContigID=%s, minContigID=%s, needToKnowNoOfLoci=%s, minNoOfLociInVCF=%s, includeIndelVCF=%s, notToUseDBToInferVCFNoOfLoci=%s ..."%\
						(inputDir, maxContigID, minContigID, needToKnowNoOfLoci, minNoOfLociInVCF, includeIndelVCF, notToUseDBToInferVCFNoOfLoci))
		returnData = PassingData(jobDataLs = [])
		counter = 0
		real_counter = 0
		if inputDir and os.path.isdir(inputDir):	#2013.04.07
			fnameLs = os.listdir(inputDir)
			previous_reported_real_counter = ''
			for fname in fnameLs:
				counter += 1
				inputFname = os.path.realpath(os.path.join(inputDir, fname))
				if (maxContigID is not None and maxContigID!=0) or (minContigID is not None and minContigID!=0):
					try:
						contigID = int(self.getContigIDFromFname(os.path.basename(fname)))
						if (maxContigID is not None and maxContigID!=0) and contigID>maxContigID:
							continue
						if (minContigID is not None and minContigID!=0) and contigID<minContigID:
							continue
					except:
						sys.stderr.write('Except type in handling file %s: %s\n'%(inputFname, repr(sys.exc_info())))
						import traceback
						traceback.print_exc()
						continue	#2013.08.12 ignore VCFs that contig/chromosome IDs could not be derived. 
				if NextGenSeq.isFileNameVCF(fname, includeIndelVCF=includeIndelVCF) and \
						not NextGenSeq.isVCFFileEmpty(inputFname, checkContent=checkEmptyVCFByReading):
					inputBaseFname = os.path.basename(inputFname)
					inputF = File(os.path.join(pegasusFolderName, inputBaseFname))
					inputF.addPFN(PFN("file://" + inputFname, input_site_handler))
					inputF.absPath = inputFname
					inputF.abspath = inputFname
					no_of_loci = None
					no_of_individuals = None
					if needToKnowNoOfLoci:
						if notToUseDBToInferVCFNoOfLoci:	#2013.06.21 do not use db to infer no of loci
							pass
						elif db_vervet:
							genotype_file = db_vervet.parseGenotypeFileGivenDBAffiliatedFilename(filename=inputFname)
							if genotype_file and inputFname.find(genotype_file.path)>=0:	#2012.9.6 make sure same file
								no_of_loci = genotype_file.no_of_loci
								no_of_individuals = genotype_file.no_of_individuals
						if no_of_loci is None:
							#do file parsing
							vcfFile = VCFFile(inputFname=inputFname, report=False)
							no_of_loci = vcfFile.getNoOfLoci()
							no_of_individuals = len(vcfFile.getSampleIDList())
							vcfFile.close()
					inputF.noOfLoci = no_of_loci
					inputF.no_of_loci = no_of_loci
					inputF.no_of_individuals = no_of_individuals
					inputF.noOfIndividuals = no_of_individuals
					
					if minNoOfLociInVCF is None or inputF.noOfLoci is None or \
							(minNoOfLociInVCF is not None and inputF.noOfLoci is not None and  inputF.noOfLoci >=minNoOfLociInVCF):
						workflow.addFile(inputF)
						tbi_F_absPath = "%s.tbi"%inputFname
						if os.path.isfile(tbi_F_absPath):	#it exists
							tbi_F = File(os.path.join(pegasusFolderName, "%s.tbi"%inputBaseFname))
							tbi_F.addPFN(PFN("file://" + tbi_F_absPath, input_site_handler))
							tbi_F.abspath = tbi_F_absPath
							workflow.addFile(tbi_F)
						else:
							tbi_F = None
						inputF.tbi_F = tbi_F
						returnData.jobDataLs.append(PassingData(job=None, jobLs=[], \
													vcfFile=inputF, tbi_F=tbi_F, file=inputF, fileLs=[inputF, tbi_F]))
						real_counter += 1
					if real_counter%10==0:
						sys.stderr.write("%s%s"%('\x08'*len(previous_reported_real_counter), real_counter))
						previous_reported_real_counter = repr(real_counter)
		sys.stderr.write("  %s non-empty VCF out of %s files, real_counter=%s.\n"%(len(returnData.jobDataLs), \
																				counter, real_counter))
		return returnData
	
	def addSplitVCFFileJob(self, workflow=None, executable=None, inputFile=None, outputFnamePrefix=None, \
					noOfOverlappingSites=1000, noOfSitesPerUnit=5000, noOfTotalSites=10000, \
					parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True, job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		2012.8.26
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		#turn them into nonnegative	
		noOfSitesPerUnit = abs(noOfSitesPerUnit)
		noOfOverlappingSites = abs(noOfOverlappingSites)
		if noOfTotalSites is not None:
			noOfTotalSites = abs(noOfTotalSites)
		else:
			noOfTotalSites = 10000000	#make it really big ,so no split
			noOfSitesPerUnit = noOfTotalSites
		
		if noOfSitesPerUnit>noOfTotalSites:
			noOfSitesPerUnit = noOfTotalSites
		if noOfOverlappingSites>noOfSitesPerUnit:
			noOfOverlappingSites = noOfSitesPerUnit
		key2ObjectForJob = {}
		extraArgumentList = ["-O %s"%outputFnamePrefix,]
		extraOutputLs = []
		if noOfOverlappingSites is not None:
			extraArgumentList.append('--noOfOverlappingSites %s'%(noOfOverlappingSites))
		if noOfSitesPerUnit is not None:
			extraArgumentList.append('--noOfSitesPerUnit %s'%(noOfSitesPerUnit))
		if noOfTotalSites is not None:
			extraArgumentList.append('--noOfTotalSites %s'%(noOfTotalSites))			
		noOfUnits = max(1, utils.getNoOfUnitsNeededToCoverN(N=noOfTotalSites, s=noOfSitesPerUnit, o=noOfOverlappingSites)-1)
		
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 
		for i in xrange(1, noOfUnits+1):
			suffixAndNameTupleList.append(['_unit%s.vcf'%(i), 'unit%s'%(i)])
		if extraArguments:
			extraArgumentList.append(extraArguments)
		self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, suffixAndNameTupleList=suffixAndNameTupleList, \
													extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
		
		job = self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		job.noOfSitesPerUnit = noOfSitesPerUnit	#2013.05.21 add this
		job.noOfUnits  = noOfUnits	#2013.06.19 added this for easy access
		return job
	
	
	def addLigateVcfJob(self, executable=None, ligateVcfExecutableFile=None, outputFile=None, \
						parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2013.07.04 moved from AlignmentToCallPipeline.py
		2012.06.26
			rename argument ligateVcfPerlPath to ligateVcfExecutableFile
		2012.6.1
			ligateVcf ligates overlapping VCF files.
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = [ligateVcfExecutableFile, outputFile]
		extraDependentInputLs.append(ligateVcfExecutableFile)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		#do not pass outputFile as argument to addGenericJob() because addGenericJob will add "-o" in front of it.
		return self.addGenericJob(executable=executable, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=[outputFile], \
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		
	
	def getChr2IntervalDataLsBySelectVCFFile(self, vcfFname=None, noOfSitesPerUnit=5000, noOfOverlappingSites=1000, \
											folderName=None, parentJobLs= None):
		"""
		2012.8.9 update it so that the interval encompassing all lines in one block/unit is known.
			good for mpileup to only work on that interval and then "bcftools view" select from sites from the block.
			TODO: offer partitioning by equal-chromosome span, rather than number of sites.
				Some sites could be in far from each other in one block, which could incur long-running mpileup. goal is to skip these deserts.
		2012.8.8 bugfix add -1 to the starting number below cuz otherwise it's included in the next block's start
				blockStopLineNumber = min(startLineNumber+(i+1)*noOfSitesPerUnit-1, stopLineNumber)	
		2012.8.14
			1.
			2. folderName is the relative path of the folder in the pegasus workflow, that holds vcfFname.
				it'll be created upon file stage-in. no mkdir job for it.
				
			get the number of lines in vcfFname.
			get chr2StartStopDataLsTuple
			for each chromosome, split its lines into units  that don't exceed noOfSitesPerUnit
				add the split job
				 
		"""
		sys.stderr.write("Splitting %s into blocks of %s lines ... "%(vcfFname, noOfSitesPerUnit))
		#from pymodule import utils
		#noOfLines = utils.getNoOfLinesInOneFileByWC(vcfFname)
		chr2StartStopDataLs = {}
		reader = MatrixFile(vcfFname)
		#csv.reader(inf, delimiter=figureOutDelimiter(vcfFname))
		lineNumber = 0
		previousChromosome = None
		previousLine = None
		chromosome = None
		for row in reader:
			lineNumber += 1
			chromosome, start, stop = row[:3]
			start = int(start)	#0-based, starting base
			stop = int(stop)	#0-based, stopping base but not inclusive, i.e. [start, stop)
			
			if previousLine is None or chromosome!=previousLine.chromosome:	#first line or different chromosome
				if previousLine is not None and previousLine.chromosome is not None:
					
					prevChrLastStartStopData = chr2StartStopDataLs[previousLine.chromosome][-1]
					if prevChrLastStartStopData.stopLineNumber is None:
						prevChrLastStartStopData.stopLineNumber = previousLine.lineNumber
						prevChrLastStartStopData.stopLineStart = previousLine.start
						prevChrLastStartStopData.stopLineStop = previousLine.stop
					
				if chromosome not in chr2StartStopDataLs:
					StartStopData = PassingData(startLineNumber=lineNumber, startLineStart=start, startLineStop=stop, \
											stopLineNumber=None, stopLineStart=None, stopLineStop=None)
					chr2StartStopDataLs[chromosome] = [StartStopData]
			else:	#same chromosome and not first line
				lastStartStopData = chr2StartStopDataLs[chromosome][-1]
				if lastStartStopData.stopLineNumber is None:	#last block hasn't been closed yet.
					noOfLinesInCurrentBlock = lineNumber - lastStartStopData.startLineNumber +1
					if noOfLinesInCurrentBlock>=noOfSitesPerUnit:	#time to close it
						lastStartStopData.stopLineNumber = lineNumber
						lastStartStopData.stopLineStart = start
						lastStartStopData.stopLineStop = stop
				else:	#generate a new block
					StartStopData = PassingData(startLineNumber=lineNumber, startLineStart=start, startLineStop=stop, \
											stopLineNumber=None, stopLineStart=None, stopLineStop=None)
					chr2StartStopDataLs[chromosome].append(StartStopData)
			previousLine = PassingData(chromosome = chromosome, start=start, stop=stop, lineNumber=lineNumber)
		#final closure
		if previousLine is not None:	#vcfFname is not empty
			lastStartStopData = chr2StartStopDataLs[previousLine.chromosome][-1]
			if lastStartStopData.stopLineNumber is None:	#last block hasn't been closed yet.
				#close it regardless of whether it has enough lines in it or not.
				lastStartStopData.stopLineNumber = previousLine.lineNumber
				lastStartStopData.stopLineStart = previousLine.start
				lastStartStopData.stopLineStop = previousLine.stop
		sys.stderr.write("%s chromosomes out of %s lines.\n"%(len(chr2StartStopDataLs), lineNumber))
		
		intervalFile = self.registerOneInputFile(inputFname=vcfFname, folderName=folderName)
		chr2IntervalDataLs = {}
		counter = 0
		for chromosome, startStopDataLs in chr2StartStopDataLs.iteritems():
			for startStopData in startStopDataLs:
				blockStartLineNumber = startStopData.startLineNumber
				blockStopLineNumber = startStopData.stopLineNumber
				# 2012.8.9 the large interval that encompasses all BED lines 
				interval = '%s:%s-%s'%(chromosome, startStopData.startLineStart, startStopData.stopLineStop)
				blockIntervalFile = File(os.path.join(folderName, '%s_line_%s_%s_bed.tsv'%(chromosome, \
												blockStartLineNumber, blockStopLineNumber)))
				blockIntervalJob = self.addSelectLineBlockFromFileJob(executable=self.SelectLineBlockFromFile, \
						inputFile=intervalFile, outputFile=blockIntervalFile,\
						startLineNumber=blockStartLineNumber, stopLineNumber=blockStopLineNumber, \
						parentJobLs=parentJobLs, extraDependentInputLs=None, \
						transferOutput=False, job_max_memory=500)
				intervalFileBasenameSignature = '%s_%s_%s'%(chromosome, blockStartLineNumber, blockStopLineNumber)
				if chromosome not in chr2IntervalDataLs:
					chr2IntervalDataLs[chromosome] = []
				intervalData = PassingData(file=blockIntervalFile, intervalFileBasenameSignature=intervalFileBasenameSignature, \
										interval=interval,\
										chr=chromosome, chromosome=chromosome, jobLs=[blockIntervalJob], job=blockIntervalJob)
				chr2IntervalDataLs[chromosome].append(intervalData)
				counter += 1
		sys.stderr.write("%s intervals and %s SelectLineBlockFromFile jobs.\n"%(counter, counter))
		return chr2IntervalDataLs
	
	
	def connectDB(self):
		"""
		2012.9.24
			place holder. AbstractVervetMapper.py will use it 
		"""
		self.registerReferenceData = None
		self.refFastaFList= None
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2013.06.14
			move topOutputDirJob from addAllJobs to here. 
		2012.9.17
		"""
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
		self.reduceOutputDirJob = reduceOutputDirJob
		
		return returnData
	
	def mapEachChromosome(self, workflow=None, chromosome=None,\
				VCFJobData=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachInterval(self, workflow=None, VCFJobData=None, chromosome=None,intervalData=None,\
					mapEachChromosomeData=None, passingData=None, transferOutput=False, \
					**keywords):
		"""
		2012.9.22
			argument VCFJobData looks like PassingData(file=splitVCFFile, vcfFile=splitVCFFile, fileLs=[splitVCFFile], \
																		job=splitVCFJob, jobLs=[splitVCFJob], tbi_F=None)
		"""
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		passingData.intervalFileBasenamePrefix
		passingData.splitVCFFile
		passingData.unitNumber
		"""
		## 2013.06.19 structures available from passingData, specific to the interval
		passingData.splitVCFFile = splitVCFFile
		passingData.unitNumber = unitNumber
		passingData.intervalFileBasenamePrefix = '%s_%s_splitVCF_u%s'%(chromosome, commonPrefix, unitNumber)
		passingData.noOfIndividuals = jobData.file.noOfIndividuals
		passingData.span = self.intervalSize + self.intervalOverlapSize*2 	#2013.06.19 for memory/walltime gauging
		"""
		return returnData
	
	def mapEachVCF(self, workflow=None, chromosome=None,VCFJobData=None, passingData=None, transferOutput=False, **keywords):
		"""
		2013.06.20 new data structure available used to split/select VCF
			intervalDataLs = passingData.intervalDataLs
		2013.04.07 use VCFJobData
		2012.9.22
			default is to split each VCF into intervals
		"""
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		fileBasenamePrefix = passingData.fileBasenamePrefix
		VCFJobData = passingData.VCFJobData
		VCFFile = VCFJobData.file	#2013.04.08
		jobData = passingData.jobData
		intervalOverlapSize = passingData.intervalOverlapSize
		intervalSize = passingData.intervalSize
		
		#intervalDataLs = passingData.intervalDataLs	#2013.06.20 new data structure available used to split/select VCF
		outputFnamePrefix = os.path.join(topOutputDirJob.output, '%s_splitVCF'%fileBasenamePrefix)
		
		splitVCFJob = self.addSplitVCFFileJob(executable=self.SplitVCFFile, inputFile=VCFFile, \
											outputFnamePrefix=outputFnamePrefix, \
				noOfOverlappingSites=intervalOverlapSize, noOfSitesPerUnit=intervalSize, \
				noOfTotalSites=getattr(VCFFile, 'noOfLoci', None), \
				parentJobLs=jobData.jobLs+[topOutputDirJob], \
				extraDependentInputLs=[jobData.tbi_F], \
				extraArguments=None, transferOutput=transferOutput, job_max_memory=2000)
		self.no_of_jobs +=1
		returnData.jobDataLs.append(PassingData(jobLs=[splitVCFJob], file=splitVCFJob.output, \
											fileLs=splitVCFJob.outputLs))
		returnData.splitVCFJob = splitVCFJob
		
		return returnData
	
	def linkMapToReduce(self, workflow=None, mapEachIntervalData=None, preReduceReturnData=None, passingData=None, transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def reduceEachChromosome(self, workflow=None, chromosome=None, passingData=None, mapEachVCFDataLs=None, reduceEachVCFDataLs=None,\
						transferOutput=True, \
						**keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachVCFDataLs = mapEachVCFDataLs
		returnData.reduceEachVCFDataLs = reduceEachVCFDataLs
		return returnData
	
	def concatenateOverlapIntervalsIntoOneVCFSubWorkflow(self, fileBasenamePrefix=None,\
					passingData=None, intervalJobLs=None,\
					outputDirJob=None,
					transferOutput=True, job_max_memory=None, walltime=None, needBGzipAndTabixJob=True, **keywords):
		"""
		Examples:
			realInputVolume = passingData.jobData.file.noOfIndividuals * passingData.jobData.file.noOfLoci
			baseInputVolume = 200*2000000
			#base is 4X coverage in 20Mb region => 120 minutes
			walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
								baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
								minJobPropertyValue=60, maxJobPropertyValue=500).value
			#base is 4X, => 5000M
			job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
								baseInputVolume=baseInputVolume, baseJobPropertyValue=2000, \
								minJobPropertyValue=2000, maxJobPropertyValue=8000).value
			self.concatenateOverlapIntervalsIntoOneVCFSubWorkflow(passingData=passingData, \
						intervalJobLs=[pdata.beagleJob for pdata in mapEachIntervalDataLs],\
						outputDirJob=self.beagleReduceDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		
		2013.07.09 added argument needBGzipAndTabixJob
		2013.06.14
			#. concatenate overlapping (share some loci) VCFs into one, used in reduceEachVCF
		"""
		#ligate vcf job (different segments of a chromosome into one chromosome) for replicate VCFs
		if fileBasenamePrefix is None:
			fileBasenamePrefix = getattr(passingData, 'fileBasenamePrefix', None)
		concatVCFFilename = os.path.join(outputDirJob.folder, '%s.vcf'%(fileBasenamePrefix))
		concatVCFFile = File(concatVCFFilename)
		if needBGzipAndTabixJob:
			transferConcatOutput=False
		else:
			transferConcatOutput = transferOutput
		concatJob = self.addLigateVcfJob(executable=self.ligateVcf, \
									ligateVcfExecutableFile=self.ligateVcfExecutableFile, \
									outputFile=concatVCFFile, \
									parentJobLs=[outputDirJob], \
									extraDependentInputLs=None, transferOutput=transferConcatOutput, \
									extraArguments=None, job_max_memory=job_max_memory, walltime=walltime/2)
		
		for intervalJob in intervalJobLs:
			#add this output to the union job
			# 2012.6.1 done it through addInputToStatMergeJob()
			self.addInputToStatMergeJob(statMergeJob=concatJob, inputF=intervalJob.output, \
							parentJobLs=[intervalJob], extraDependentInputLs=intervalJob.outputLs[1:])
		
		if needBGzipAndTabixJob:
			#bgzip and tabix the trio caller output
			gzipVCFFile = File("%s.gz"%concatVCFFilename)
			bgzip_tabix_job = self.addBGZIP_tabix_Job(bgzip_tabix=self.bgzip_tabix_in_reduce, \
									parentJob=concatJob, \
									inputF=concatJob.output, outputF=gzipVCFFile, transferOutput=transferOutput,\
									job_max_memory=job_max_memory/4, walltime=walltime/4)
			
			returnData = PassingData(file=bgzip_tabix_job.output, vcfFile=bgzip_tabix_job.output, \
					fileLs=[bgzip_tabix_job.output, bgzip_tabix_job.tbi_F], \
					job=bgzip_tabix_job, jobLs=[bgzip_tabix_job], tbi_F=bgzip_tabix_job.tbi_F)
		else:
			returnData = PassingData(file=concatJob.output, vcfFile=concatJob.output, \
					fileLs=concatJob.outputLs, \
					job=concatJob, jobLs=[concatJob], tbi_F=None)
		return returnData
	
	def concatenateIntervalsIntoOneVCFSubWorkflow(self, executable=None, refFastaFList=None, fileBasenamePrefix=None, \
					passingData=None, intervalJobLs=None,\
					outputDirJob=None,needBGzipAndTabixJob=True,\
					transferOutput=True, job_max_memory=None, walltime=None, \
					**keywords):
		"""
		2013.07.08
			#. concatenate non-overlapping VCFs into one, used in reduceEachVCF or reduce()
			assuming all VCFs with same set of samples
		"""
		
		if executable is None:
			executable = self.CombineVariantsJavaInReduce
		#ligate vcf job (different segments of a chromosome into one chromosome) for replicate VCFs
		if fileBasenamePrefix is None:
			fileBasenamePrefix = getattr(passingData, 'fileBasenamePrefix', None)
		
		if fileBasenamePrefix is None:
			message = "ERRROR in AbstractVCFWorkflow.concatenateIntervalsIntoOneVCFSubWorkflow(): %s is None."%\
							(fileBasenamePrefix)
			utils.pauseForUserInput(message=message, continueAnswerSet=None, exitType=1)
		
		concatVCFFilename = os.path.join(outputDirJob.folder, '%s.vcf'%(fileBasenamePrefix))
		concatVCFFile = File(concatVCFFilename)
		if needBGzipAndTabixJob:
			transferConcatOutput=False
		else:
			transferConcatOutput = transferOutput
		concatJob = self.addGATKCombineVariantsJob(executable=executable, GenomeAnalysisTKJar=None, \
							refFastaFList=refFastaFList, inputFileList=None, argumentForEachFileInInputFileList="--variant",\
							outputFile=concatVCFFile, genotypeMergeOptions='UNSORTED', \
					parentJobLs=[outputDirJob], transferOutput=transferConcatOutput, job_max_memory=job_max_memory,\
					walltime=walltime,\
					extraArguments=None, extraArgumentList=['--assumeIdenticalSamples'], extraDependentInputLs=None)
		
		for intervalJob in intervalJobLs:
			self.addInputToStatMergeJob(statMergeJob=concatJob, inputF=intervalJob.output, inputArgumentOption="--variant",\
							parentJobLs=[intervalJob], extraDependentInputLs=intervalJob.outputLs[1:])
		
		
		if needBGzipAndTabixJob:
			
			#bgzip and tabix the trio caller output
			gzipVCFFile = File("%s.gz"%concatVCFFilename)
			bgzip_tabix_job = self.addBGZIP_tabix_Job(bgzip_tabix=self.bgzip_tabix_in_reduce, \
									parentJob=concatJob, \
									inputF=concatJob.output, outputF=gzipVCFFile, transferOutput=transferOutput,\
									job_max_memory=job_max_memory/4, walltime=walltime/4)
			
			returnData = self.constructJobDataFromJob(bgzip_tabix_job)
		else:
			returnData = self.constructJobDataFromJob(concatJob)
		return returnData
	
	def reduceManySmallIntervalVCFIntoBigIntervalVCF(self, executable=None, registerReferenceData=None, fileBasenamePrefix=None, \
					intervalJobLs=None, outputDirJob=None, bigIntervalSize=10000000,
					transferOutput=False, job_max_memory=None, walltime=None, needBGzipAndTabixJob=False,\
					**keywords):
		"""
		
			#. each intervalJob should have an intervalData attached to it , detailing how big an interval it is covering
			#. assuming all jobs in intervalJobLs are in chromosomal order
		
		2013.09.05
			
		"""
		sys.stderr.write("Reducing %s interval VCFs into VCFs spanning %s bp, current job count=%s ..."%\
						(len(intervalJobLs), bigIntervalSize, self.no_of_jobs))
		returnData = PassingData(jobDataLs=[])
		_intervalJobLs = []
		startIntervalData = None	#an instance of IntervalData in AbstractNGSWorkflow()
		for intervalJob in intervalJobLs:
			intervalData = intervalJob.intervalData
			if startIntervalData is None:
				_intervalJobLs.append(intervalJob)
				startIntervalData = copy.deepcopy(intervalData)	#not interfering original data
			elif intervalData.chromosome!=startIntervalData.chromosome or (intervalData.span+startIntervalData.span)>bigIntervalSize:
				_filenamePartLs = []
				if fileBasenamePrefix:
					_filenamePartLs.append(fileBasenamePrefix)
				_filenamePartLs.extend([startIntervalData.chromosome, \
										startIntervalData.overlapStart, startIntervalData.overlapStop])
				_filenamePartLs = map(str, _filenamePartLs)
				_fileBasenamePrefix = "_".join(_filenamePartLs)
				concateJobData = self.concatenateIntervalsIntoOneVCFSubWorkflow(executable=executable, \
						refFastaFList=registerReferenceData.refFastaFList, \
						fileBasenamePrefix=_fileBasenamePrefix, \
						passingData=None, intervalJobLs=_intervalJobLs,\
						outputDirJob=outputDirJob,
						needBGzipAndTabixJob=needBGzipAndTabixJob,\
						transferOutput=transferOutput, \
						job_max_memory=job_max_memory, walltime=walltime)
				concateJobData.intervalData = startIntervalData
				returnData.jobDataLs.append(concateJobData)
				
				#reset
				_intervalJobLs = [intervalJob]
				startIntervalData = copy.deepcopy(intervalData)	#not interfering original data
			else:
				_intervalJobLs.append(intervalJob)
				startIntervalData.unionOneIntervalData(intervalData)
		if _intervalJobLs:
			_fileBasenamePrefix = "%s_%s_%s_%s"%(fileBasenamePrefix, startIntervalData.chromosome, \
													startIntervalData.overlapStart, startIntervalData.overlapStop)
			concateJobData = self.concatenateIntervalsIntoOneVCFSubWorkflow(executable=executable, \
									refFastaFList=registerReferenceData.refFastaFList, \
					fileBasenamePrefix=_fileBasenamePrefix, \
					passingData=None, intervalJobLs=_intervalJobLs,\
					outputDirJob=outputDirJob,
					needBGzipAndTabixJob=needBGzipAndTabixJob,\
					transferOutput=transferOutput, job_max_memory=job_max_memory, walltime=walltime)
			concateJobData.intervalData = startIntervalData
			returnData.jobDataLs.append(concateJobData)
		sys.stderr.write(" reduced to %s jobs, current job count=%s.\n"%\
						(len(returnData.jobDataLs), self.no_of_jobs))
		return returnData
	
	def reduceEachVCF(self, workflow=None, chromosome=None, passingData=None, mapEachIntervalDataLs=None,\
					transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
		"""
		realInputVolume = passingData.jobData.file.noOfIndividuals * passingData.jobData.file.noOfLoci
		baseInputVolume = 200*2000000
		#base is 4X coverage in 20Mb region => 120 minutes
		walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value
		#base is 4X, => 5000M
		job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=2000, \
							minJobPropertyValue=2000, maxJobPropertyValue=8000).value
		concatenateIntoOneVCFJob = self.concatenateOverlapIntervalsIntoOneVCFSubWorkflow(chromosome=chromosome, passingData=passingData, \
						intervalJobLs=[pdata.lastVCFJob for pdata in mapEachIntervalDataLs],\
						outputDirJob=self.reduceOutputDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		"""
		return returnData
	
	def reduce(self, workflow=None, reduceEachChromosomeDataLs=None, \
			mapEachChromosomeDataLs=None, passingData=None, transferOutput=True, \
			**keywords):
		"""
		2013.07.18 return each processed-VCF job data so that followup workflows could carry out map-reduce
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachChromosomeDataLs = mapEachChromosomeDataLs
		returnData.reduceEachChromosomeDataLs = reduceEachChromosomeDataLs
		"""
		#2013.07.18 example to return each processed-VCF job data so that followup workflows could carry out map-reduce
		for reduceEachVCFDataLs in passingData.reduceEachVCFDataLsLs:
			if reduceEachVCFDataLs:
				for reduceEachVCFData in reduceEachVCFDataLs:
					if reduceEachVCFData:
						returnData.jobDataLs.append(reduceEachVCFData.WHATEVERJobData)
		"""
		return returnData
	
	def addAllJobs(self, workflow=None, inputVCFData=None, chr2IntervalDataLs=None, \
				GenomeAnalysisTKJar=None, samtools=None, \
				CreateSequenceDictionaryJava=None, CreateSequenceDictionaryJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexJar=None,\
				mv=None, \
				registerReferenceData=None, \
				needFastaIndexJob=False, needFastaDictJob=False, \
				data_dir=None, no_of_gatk_threads = 1, \
				intervalSize=3000, intervalOverlapSize=0, \
				outputDirPrefix="", passingData=None, \
				transferOutput=True, job_max_memory=2000, **keywords):
		"""
		2013.06.14 bugfix regarding noOfUnits, which was all inferred from one file
		2012.7.26
			architect of the whole map-reduce framework
		"""
		chr2jobDataLs = {}
		for jobData in inputVCFData.jobDataLs:
			chromosome = self.getChrFromFname(os.path.basename(jobData.file.name))
			if chromosome not in chr2jobDataLs:
				chr2jobDataLs[chromosome] = []
			chr2jobDataLs[chromosome].append(jobData)
		
		sys.stderr.write("Adding jobs for %s chromosomes/contigs of %s VCF files... \n"%(len(chr2jobDataLs), len(inputVCFData.jobDataLs)))
		if getattr(registerReferenceData, 'refFastaFList', None):
			refFastaFList = registerReferenceData.refFastaFList
		else:
			refFastaFList = None
		if refFastaFList:
			refFastaF = refFastaFList[0]
		else:
			refFastaF = None
		
		if needFastaDictJob or getattr(registerReferenceData, 'needPicardFastaDictJob', None):
			fastaDictJob = self.addRefFastaDictJob(workflow, CreateSequenceDictionaryJava=CreateSequenceDictionaryJava, \
												refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
		else:
			fastaDictJob = None
			refFastaDictF = getattr(registerReferenceData, 'refPicardFastaDictF', None)
		
		if needFastaIndexJob or getattr(registerReferenceData, 'needSAMtoolsFastaIndexJob', None):
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
		else:
			fastaIndexJob = None
			refFastaIndexF = getattr(registerReferenceData, 'refSAMtoolsFastaIndexF', None)
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		#2012.9.22 
		# 	mapEachAlignmentDataLs is never reset.
		#	mapEachChromosomeDataLs is reset upon new alignment
		#	mapEachIntervalDataLs is reset upon each new chromosome
		#	all reduce lists never get reset.
		#	fileBasenamePrefix is the prefix of input file's basename, to be used for temporary output files in reduceEachVCF()
		#		but not for output files in mapEachInterval()
		passingData = PassingData(\
					fileBasenamePrefix=None, \
					chromosome=None, \
					
					outputDirPrefix=outputDirPrefix, \
					intervalFileBasenamePrefix=None,\
					
					refFastaFList=refFastaFList, \
					registerReferenceData=registerReferenceData, \
					refFastaF=refFastaF,\
					
					fastaDictJob = fastaDictJob,\
					refFastaDictF = refFastaDictF,\
					fastaIndexJob = fastaIndexJob,\
					refFastaIndexF = refFastaIndexF,\
					
					intervalOverlapSize =intervalOverlapSize, intervalSize=intervalSize,\
					jobData=None,\
					VCFJobData=None,\
					splitVCFFile=None,\
					intervalDataLs=None,\
					preReduceReturnData=None,\
					
					mapEachIntervalData=None,\
					mapEachVCFData=None,\
					mapEachChromosomeData=None, \
					
					mapEachIntervalDataLs=None,\
					mapEachVCFDataLs=None,\
					
					mapEachIntervalDataLsLs=[],\
					mapEachVCFDataLsLs=[],\
					
					mapEachChromosomeDataLs=[], \
					
					reduceEachVCFData=None,\
					reduceEachChromosomeData=None,\
					
					reduceEachVCFDataLs=None,\
					
					reduceEachVCFDataLsLs=[],\
					
					reduceEachChromosomeDataLs=[],\
					
					chr2jobDataLs = chr2jobDataLs,\
					)
		# mapEachIntervalDataLsLs is list of mapEachIntervalDataLs by each VCF file.
		# mapEachVCFDataLsLs is list of mapEachVCFDataLs by each chromosome
		# reduceEachVCFDataLsLs is list of reduceEachVCFDataLs by each chromosome
		
		preReduceReturnData = self.preReduce(workflow=workflow, outputDirPrefix=outputDirPrefix, \
									passingData=passingData, transferOutput=True,\
									**keywords)
		passingData.preReduceReturnData = preReduceReturnData
		
		#gzip folder jobs (to avoid repeatedly creating the same folder
		gzipReduceEachVCFFolderJob = None
		gzipReduceEachChromosomeFolderJob = None
		gzipReduceFolderJob = None
		gzipPreReduceFolderJob = None
		no_of_vcf_files = 0
		for chromosome, jobDataLs in chr2jobDataLs.iteritems():
			passingData.chromosome = chromosome
			mapEachChromosomeData = self.mapEachChromosome(workflow=workflow, chromosome=chromosome, \
										passingData=passingData, \
										transferOutput=False, **keywords)
			passingData.mapEachChromosomeData = mapEachChromosomeData
			passingData.mapEachChromosomeDataLs.append(mapEachChromosomeData)
			
			passingData.mapEachVCFDataLsLs.append([])
			#the last one from the double list is the current one
			passingData.mapEachVCFDataLs = passingData.mapEachVCFDataLsLs[-1]
			
			passingData.reduceEachVCFDataLsLs.append([])
			passingData.reduceEachVCFDataLs = passingData.reduceEachVCFDataLsLs[-1]
			
			#2013.06.20
			if chr2IntervalDataLs:
				passingData.intervalDataLs = chr2IntervalDataLs.get(chromosome)
			else:
				pass
			for i in xrange(len(jobDataLs)):
				jobData = jobDataLs[i]
				passingData.jobData = jobData
				passingData.VCFJobData = jobData
				
				VCFFile = jobData.vcfFile
				inputFBaseName = os.path.basename(VCFFile.name)
				commonPrefix = inputFBaseName.split('.')[0]
				
				passingData.fileBasenamePrefix = commonPrefix
				
				no_of_vcf_files += 1
				if no_of_vcf_files%10==0:
					sys.stderr.write("%s\t%s VCFs."%('\x08'*40, no_of_vcf_files))
				if self.mapReduceType==1:
					mapEachVCFData = self.mapEachVCF(workflow=workflow,\
													VCFJobData=jobData, passingData=passingData, \
													transferOutput=False, **keywords)
					passingData.mapEachVCFData = mapEachVCFData
					passingData.mapEachVCFDataLs.append(mapEachVCFData)
					
					passingData.mapEachIntervalDataLsLs.append([])
					passingData.mapEachIntervalDataLs = passingData.mapEachIntervalDataLsLs[-1]
					
					splitVCFJob = mapEachVCFData.splitVCFJob
					#noOfUnits = max(1, utils.getNoOfUnitsNeededToCoverN(N=jobData.file.noOfLoci, s=intervalSize, o=intervalOverlapSize)-1)
					noOfUnits = splitVCFJob.noOfUnits
					jobData.file.noOfUnits = noOfUnits	#2013.06.14
					for unitNumber in xrange(1, noOfUnits+1):
						splitVCFFile = getattr(splitVCFJob, 'unit%sFile'%(unitNumber), None)
						if splitVCFFile is not None:
							passingData.splitVCFFile = splitVCFFile
							passingData.unitNumber = unitNumber
							passingData.intervalFileBasenamePrefix = '%s_%s_splitVCF_u%s'%(chromosome, commonPrefix, unitNumber)
							passingData.noOfIndividuals = getattr(jobData.file, 'noOfIndividuals', None)
							passingData.span = self.intervalSize + self.intervalOverlapSize*2 	#2013.06.19 for memory/walltime gauging
							mapEachIntervalData = self.mapEachInterval(workflow=workflow, \
													VCFJobData=PassingData(file=splitVCFFile, vcfFile=splitVCFFile, fileLs=[splitVCFFile], \
																		job=splitVCFJob, jobLs=[splitVCFJob], tbi_F=None), \
													chromosome=chromosome,intervalData=None,\
													mapEachChromosomeData=mapEachChromosomeData, \
													passingData=passingData, transferOutput=False, \
													**keywords)
							passingData.mapEachIntervalData = mapEachIntervalData
							passingData.mapEachIntervalDataLs.append(mapEachIntervalData)
							
							linkMapToReduceData = self.linkMapToReduce(workflow=workflow, mapEachIntervalData=mapEachIntervalData, \
												preReduceReturnData=preReduceReturnData, \
												passingData=passingData, \
												**keywords)
				else:
					if not chr2IntervalDataLs:
						sys.stderr.write("Error: mapReduceType=%s, but chr2IntervalDataLs (%s) is nothing.\n"%\
										(self.mapReduceType, repr(chr2IntervalDataLs)))
						sys.exit(2)
					if not chr2IntervalDataLs.get(chromosome):
						sys.stderr.write("Error: chromosome %s is not in chr2IntervalDataLs (mapReduceType=%s) while VCFFile %s is from this chromosome.\n"%\
										(chromosome, self.mapReduceType, jobData.file.name))
						sys.exit(3)
					for intervalData in chr2IntervalDataLs.get(chromosome):
						mapEachIntervalData = self.mapEachInterval(workflow=workflow, \
														VCFJobData=jobData, \
														chromosome=chromosome,intervalData=intervalData,\
														mapEachChromosomeData=mapEachChromosomeData, \
														passingData=passingData, transferOutput=False, \
														**keywords)
						passingData.mapEachIntervalData = mapEachIntervalData
						passingData.mapEachIntervalDataLs.append(mapEachIntervalData)
						
						linkMapToReduceData = self.linkMapToReduce(workflow=workflow, mapEachIntervalData=mapEachIntervalData, \
											preReduceReturnData=preReduceReturnData, \
											passingData=passingData, \
											**keywords)
				reduceEachVCFData = self.reduceEachVCF(workflow=workflow, chromosome=chromosome, passingData=passingData, \
								mapEachIntervalDataLs=passingData.mapEachIntervalDataLs,\
								transferOutput=False, data_dir=data_dir, \
								**keywords)
				passingData.reduceEachVCFData = reduceEachVCFData
				passingData.reduceEachVCFDataLs.append(reduceEachVCFData)
				
				gzipReduceEachVCFData = self.addGzipSubWorkflow(workflow=workflow, \
					inputData=reduceEachVCFData, transferOutput=transferOutput,\
					outputDirPrefix="%sReduceEachVCF"%(outputDirPrefix), topOutputDirJob=gzipReduceEachVCFFolderJob, \
					report=False)
				gzipReduceEachVCFFolderJob = gzipReduceEachVCFData.topOutputDirJob
			reduceEachChromosomeData = self.reduceEachChromosome(workflow=workflow, chromosome=chromosome, passingData=passingData, \
								mapEachVCFDataLs=passingData.mapEachVCFDataLs,\
								reduceEachVCFDataLs=passingData.reduceEachVCFDataLs,\
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
							mapEachChromosomeDataLs=passingData.mapEachVCFDataLs,\
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
		
		sys.stderr.write("\n %s%s VCF files.\n"%('\x08'*40, no_of_vcf_files))
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return reduceReturnData
	
	def setup_run(self):
		"""
		2013.06.11 added firstVCFJobData in return
			assign all returned data to self, rather than pdata (pdata has become self)
		2013.04.07 wrap all standard pre-run() related functions into this function.
			setting up for run(), called by run()
		"""
		pdata = AbstractNGSWorkflow.setup_run(self)
		
		#self.chr2size = {}
		#self.chr2size = set(['Contig149'])	#temporary when testing Contig149
		#self.chr2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		if self.needSplitChrIntervalData:	#2013.06.21 defined in AbstractNGSWorkflow.__init__()
			chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=self.chr2size, \
													intervalSize=self.intervalSize, \
													intervalOverlapSize=self.intervalOverlapSize)
		else:
			chr2IntervalDataLs = None
		inputData = None
		firstVCFJobData = None
		if getattr(self, 'inputDir', None):	#2013.05.20 bugfix
			inputData = self.registerAllInputFiles(inputDir=self.inputDir, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName=self.pegasusFolderName,\
											maxContigID=self.maxContigID, \
											minContigID=self.minContigID,\
											db_vervet=getattr(self, 'db_vervet', None), \
											needToKnowNoOfLoci=getattr(self, 'needToKnowNoOfLoci', True),\
											minNoOfLociInVCF=getattr(self, 'minNoOfLociInVCF', 10))
			if inputData and inputData.jobDataLs:
				firstVCFJobData = inputData.jobDataLs[0]
				#job=None, jobLs=[], vcfFile=inputF, tbi_F=tbi_F, file=inputF, fileLs=[inputF, tbi_F]
				firstVCFFile = firstVCFJobData.file
				sys.stderr.write("\t VCF file %s is chosen as an example VCF for any job that needs a random VCF file.\n"%(firstVCFFile))
		
		registerReferenceData = self.getReferenceSequence()
		
		self.inputData = inputData
		self.chr2IntervalDataLs = chr2IntervalDataLs
		self.registerReferenceData = registerReferenceData
		self.firstVCFJobData = firstVCFJobData
		#self.firstVCFFile = firstVCFFile
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
				
		self.addAllJobs(workflow=workflow, inputVCFData=inputData, \
				chr2IntervalDataLs=self.chr2IntervalDataLs, samtools=workflow.samtools, \
				GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
				CreateSequenceDictionaryJava=workflow.CreateSequenceDictionaryJava, \
				CreateSequenceDictionaryJar=workflow.CreateSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar,\
				mv=workflow.mv, \
				registerReferenceData=pdata.registerReferenceData,\
				needFastaIndexJob=getattr(self, 'needFastaIndexJob',False), \
				needFastaDictJob=getattr(self, 'needFastaDictJob', False), \
				data_dir=self.data_dir, no_of_gatk_threads = 1,\
				intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
				outputDirPrefix=self.pegasusFolderName, transferOutput=True,)
		
		self.end_run()
