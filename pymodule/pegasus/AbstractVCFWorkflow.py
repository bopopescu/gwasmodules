#!/usr/bin/env python
"""
2012.1.17
	a common class for pegasus workflows that work on VCF variant files
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, utils
from Pegasus.DAX3 import *
from AbstractNGSWorkflow import AbstractNGSWorkflow

class AbstractVCFWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('minDepth', 0, float): [0, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						})
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		if getattr(self, "inputDir", None):
			self.inputDir = os.path.abspath(self.inputDir)
	
	def addAddVCFFile2DBJob(self, executable=None, inputFile=None, genotypeMethodShortName=None,\
						logFile=None, format=None, dataDir=None, checkEmptyVCFByReading=None, commit=False, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.8.30 moved from vervet/src/AddVCFFolder2DBWorkflow.py
		2012.6.27
		"""
		extraArgumentList = ['-f', format]
		if logFile:
			extraArgumentList.extend(["-l", logFile])
		if dataDir:
			extraArgumentList.extend(['-t', dataDir])
		if checkEmptyVCFByReading:
			extraArgumentList.extend(['-E'])
		if genotypeMethodShortName:
			extraArgumentList.extend(['-s', genotypeMethodShortName, ])
		if commit:
			extraArgumentList.append('-c')
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[logFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		return job
	
	def registerExecutables(self, workflow=None):
		"""
		"""
		AbstractNGSWorkflow.registerExecutables(self, workflow)
		
		if not workflow:
			workflow = self
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		PlotVCFtoolsStat = Executable(namespace=namespace, name="PlotVCFtoolsStat", version=version, os=operatingSystem, arch=architecture, installed=True)
		PlotVCFtoolsStat.addPFN(PFN("file://" +  os.path.join(self.vervetSrcPath, "plot/PlotVCFtoolsStat.py"), site_handler))
		executableClusterSizeMultiplierList.append((PlotVCFtoolsStat, 0))
		
		SplitVCFFile = Executable(namespace=namespace, name="SplitVCFFile", version=version, os=operatingSystem, arch=architecture, installed=True)
		SplitVCFFile.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "pegasus/mapper/SplitVCFFile.py"), site_handler))
		executableClusterSizeMultiplierList.append((SplitVCFFile, 1))
		
		#2012.8.30 moved from vervet/src/AddVCFFolder2DBWorkflow.py
		AddVCFFile2DB = Executable(namespace=namespace, name="AddVCFFile2DB", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		AddVCFFile2DB.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/AddVCFFile2DB.py"), site_handler))
		executableClusterSizeMultiplierList.append((AddVCFFile2DB, 1))
		
		FilterVCFSNPCluster = Executable(namespace=namespace, name="FilterVCFSNPCluster", version=version, os=operatingSystem, arch=architecture, installed=True)
		FilterVCFSNPCluster.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "pegasus/mapper/FilterVCFSNPCluster.py"), site_handler))
		executableClusterSizeMultiplierList.append((FilterVCFSNPCluster, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
	def registerCommonExecutables(self, workflow=None):
		"""
		"""
		AbstractNGSWorkflow.registerCommonExecutables(self, workflow)
	
	def registerAllInputFiles(self, workflow=None, inputDir=None, input_site_handler=None, \
					checkEmptyVCFByReading=False, pegasusFolderName='',\
					maxContigID=None, minContigID=None, db_vervet=None, needToKnowNoOfLoci=False,
					minNoOfLoci=None):
		"""
		2012.8.15 add argument db_vervet, needToKnowNoOfLoci, to get no_of_loci by parsing inputFname and find db-entry...
			argument minNoOfLoci, only used when it's not None and needToKnowNoOfLoci is True
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
		sys.stderr.write("Registering input files from %s ..."%(inputDir))
		returnData = PassingData(jobDataLs = [])
		fnameLs = os.listdir(inputDir)
		counter = 0
		real_counter = 0
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
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			if NextGenSeq.isFileNameVCF(fname, includeIndelVCF=False) and \
					not NextGenSeq.isVCFFileEmpty(inputFname, checkContent=checkEmptyVCFByReading):
				real_counter += 1
				inputBaseFname = os.path.basename(inputFname)
				inputF = File(os.path.join(pegasusFolderName, inputBaseFname))
				inputF.addPFN(PFN("file://" + inputFname, input_site_handler))
				inputF.absPath = inputFname
				inputF.abspath = inputFname
				no_of_loci = None
				no_of_individuals = None
				if needToKnowNoOfLoci:
					if db_vervet:
						genotype_file = db_vervet.parseGenotypeFileGivenDBAffiliatedFilename(filename=inputFname)
						if genotype_file and inputFname.find(genotype_file.path)>=0:	#2012.9.6 make sure same file
							no_of_loci = genotype_file.no_of_loci
							no_of_individuals = genotype_file.no_of_individuals
					if no_of_loci is None:
						#do file parsing
						from pymodule.VCFFile import VCFFile
						vcfFile = VCFFile(inputFname=inputFname)
						counter = 0
						no_of_loci = 0
						for vcfRecord in vcfFile:
							no_of_loci += 1
						no_of_individuals = len(vcfFile.getSampleIDList())
						vcfFile.close()
				inputF.noOfLoci = no_of_loci
				inputF.no_of_loci = no_of_loci
				inputF.no_of_individuals = no_of_individuals
				inputF.noOfIndividuals = no_of_individuals
				
				if minNoOfLoci is None or (minNoOfLoci and inputF.no_of_loci and  inputF.no_of_loci >minNoOfLoci):
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
					returnData.jobDataLs.append(PassingData(vcfFile=inputF, jobLs=[], tbi_F=tbi_F, file=inputF, fileLs=[]))
				if real_counter%200==0:
					sys.stderr.write("%s%s"%('\x08'*len(previous_reported_real_counter), real_counter))
					previous_reported_real_counter = repr(real_counter)
		sys.stderr.write("  %s non-empty VCF out of %s files.\n"%(len(returnData.jobDataLs), counter))
		return returnData

	def addPlotVCFtoolsStatJob(self, workflow=None, executable=None, inputFileList=None, outputFnamePrefix=None, \
							whichColumn=None, whichColumnLabel=None, whichColumnPlotLabel=None, need_svg=False, \
							logWhichColumn=True, positiveLog=False, valueForNonPositiveYValue=-1, \
							posColumnPlotLabel=None, chrLengthColumnLabel=None, chrColumnLabel=None, \
							minChrLength=1000000, posColumnLabel=None, minNoOfTotal=100,\
							figureDPI=300, ylim_type=2, samplingRate=0.0001, logCount=False,\
							parentJobLs=None, \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		2012.8.31 add argument positiveLog and valueForNonPositiveYValue
		# whichColumnPlotLabel and posColumnPlotLabel should not contain spaces or ( or ). because they will disrupt shell commandline
		
		2012.8.2 moved from vervet/src/CalculateVCFStatPipeline.py
		2012.8.1
			
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnLabel', 0, ): ["", 'W', 1, 'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take log of whichColumn'],\
			('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
			('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, 'plot label for data of the whichColumn', ],\
			('posColumnPlotLabel', 1, ): ['position', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
			('chrLengthColumnLabel', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
			('chrColumnLabel', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
			('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
			('posColumnLabel', 1, ): ['BIN_START', 'l', 1, 'label of the position column, BIN_START for binned vcftools output. POS for others.', ],\
			('outputFnamePrefix', 0, ): [None, 'O', 1, 'output filename prefix (optional).'],\
			
				('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
				('title', 1, ): [None, 't', 1, 'title for the figure.'],\
				('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
				('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
				('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
				('samplingRate', 1, float): [0.001, 's', 1, 'how often you include the data'],\
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if inputFileList:
			extraDependentInputLs.extend(inputFileList)
		extraArgumentList = ["-O %s"%outputFnamePrefix, '-i %s'%(minNoOfTotal), \
							'-f %s'%(figureDPI), '-y %s'%(ylim_type), '-s %s'%(samplingRate), '-l %s'%(posColumnLabel)]
		extraOutputLs = [File('%s.png'%(outputFnamePrefix)), File('%s_hist.png'%(outputFnamePrefix))]
		if need_svg:
			extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
		key2ObjectForJob = {}
		if minChrLength is not None:
			extraArgumentList.append('-m %s'%(minChrLength))
		if whichColumnLabel:
			extraArgumentList.append("--whichColumnLabel %s"%(whichColumnLabel))
		if whichColumn:
			extraArgumentList.append("--whichColumn %s"%(whichColumn))
		if logWhichColumn:
			extraArgumentList.append('--logWhichColumn')
			if positiveLog:
				extraArgumentList.append('--positiveLog')
		if whichColumnPlotLabel:
			extraArgumentList.append("--whichColumnPlotLabel %s"%(whichColumnPlotLabel))
		if posColumnPlotLabel:
			extraArgumentList.append("--posColumnPlotLabel %s"%(posColumnPlotLabel))
		if chrLengthColumnLabel:
			extraArgumentList.append("--chrLengthColumnLabel %s"%(chrLengthColumnLabel))
		if chrColumnLabel:
			extraArgumentList.append("--chrColumnLabel %s"%(chrColumnLabel))
		if logCount:
			extraArgumentList.append("--logCount")
		if valueForNonPositiveYValue:
			extraArgumentList.append("--valueForNonPositiveYValue %s"%(valueForNonPositiveYValue))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		#add all input files to the last (after db arguments,) otherwise, it'll mask others (cuz these don't have options).
		if inputFileList:
			for inputFile in inputFileList:
				if inputFile:
					job.addArguments(inputFile)
		return job
	
	
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
		noOfOverlappingSites = abs(noOfOverlappingSites)
		noOfSitesPerUnit = abs(noOfSitesPerUnit)
		noOfTotalSites = abs(noOfTotalSites)
		
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
		
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		return job
	
	def getChr2IntervalDataLsBySelectVCFFile(self, vcfFolder=None, noOfSitesPerUnit=5000, noOfOverlappingSites=1000, \
											folderName=None, parentJobLs= None):
		"""
		2012.8.9 update it so that the interval encompassing all lines in one block/unit is known.
			good for mpileup to only work on that interval and then "bcftools view" select from sites from the block.
			TODO: offer partitioning by equal-chromosome span, rather than number of sites.
				Some sites could be in far from each other in one block, which could incur long-running mpileup. goal is to skip these deserts.
		2012.8.8 bugfix add -1 to the starting number below cuz otherwise it's included in the next block's start
				blockStopLineNumber = min(startLineNumber+(i+1)*noOfLinesPerUnit-1, stopLineNumber)	
		2012.8.14
			1.
			2. folderName is the relative path of the folder in the pegasus workflow, that holds intervalFname.
				it'll be created upon file stage-in. no mkdir job for it.
				
			get the number of lines in intervalFname.
			get chr2StartStopDataLsTuple
			for each chr, split its lines into units  that don't exceed noOfLinesPerUnit
				add the split job
				 
		"""
		sys.stderr.write("Splitting %s into blocks of %s lines ... "%(intervalFname, noOfLinesPerUnit))
		#from pymodule import utils
		#noOfLines = utils.getNoOfLinesInOneFileByWC(intervalFname)
		chr2StartStopDataLs = {}
		import csv
		from pymodule import figureOutDelimiter
		inf = open(intervalFname)
		reader = csv.reader(inf, delimiter=figureOutDelimiter(intervalFname))
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
					if noOfLinesInCurrentBlock>=noOfLinesPerUnit:	#time to close it
						lastStartStopData.stopLineNumber = lineNumber
						lastStartStopData.stopLineStart = start
						lastStartStopData.stopLineStop = stop
				else:	#generate a new block
					StartStopData = PassingData(startLineNumber=lineNumber, startLineStart=start, startLineStop=stop, \
											stopLineNumber=None, stopLineStart=None, stopLineStop=None)
					chr2StartStopDataLs[chromosome].append(StartStopData)
			previousLine = PassingData(chromosome = chromosome, start=start, stop=stop, lineNumber=lineNumber)
		#final closure
		if previousLine is not None:	#intervalFname is not empty
			lastStartStopData = chr2StartStopDataLs[previousLine.chromosome][-1]
			if lastStartStopData.stopLineNumber is None:	#last block hasn't been closed yet.
				#close it regardless of whether it has enough lines in it or not.
				lastStartStopData.stopLineNumber = previousLine.lineNumber
				lastStartStopData.stopLineStart = previousLine.start
				lastStartStopData.stopLineStop = previousLine.stop
		sys.stderr.write("%s chromosomes out of %s lines.\n"%(len(chr2StartStopDataLs), lineNumber))
		
		intervalFile = self.registerOneInputFile(inputFname=intervalFname, folderName=folderName)
		chr2IntervalDataLs = {}
		counter = 0
		for chr, startStopDataLs in chr2StartStopDataLs.iteritems():
			for startStopData in startStopDataLs:
				blockStartLineNumber = startStopData.startLineNumber
				blockStopLineNumber = startStopData.stopLineNumber
				# 2012.8.9 the large interval that encompasses all BED lines 
				interval = '%s:%s-%s'%(chr, startStopData.startLineStart, startStopData.stopLineStop)
				blockIntervalFile = File(os.path.join(folderName, '%s_line_%s_%s_bed.tsv'%(chr, \
												blockStartLineNumber, blockStopLineNumber)))
				blockIntervalJob = self.addSelectLineBlockFromFileJob(executable=self.SelectLineBlockFromFile, \
						inputFile=intervalFile, outputFile=blockIntervalFile,\
						startLineNumber=blockStartLineNumber, stopLineNumber=blockStopLineNumber, \
						parentJobLs=parentJobLs, extraDependentInputLs=None, \
						transferOutput=False, job_max_memory=500)
				intervalFnameSignature = '%s_%s_%s'%(chr, blockStartLineNumber, blockStopLineNumber)
				if chr not in chr2IntervalDataLs:
					chr2IntervalDataLs[chr] = []
				intervalData = PassingData(file=blockIntervalFile, intervalFnameSignature=intervalFnameSignature, interval=interval,\
										chr=chr, jobLs=[blockIntervalJob], job=blockIntervalJob)
				chr2IntervalDataLs[chr].append(intervalData)
				counter += 1
		sys.stderr.write("%s intervals and %s SelectLineBlockFromFile jobs.\n"%(counter, counter))
		return chr2IntervalDataLs